function sim_fun_OFDM(conn,num_frames,M_ary,EbN0,M,T,Fc,vel)

% Make parameter structure
systemName = 'OFDM';
params_inst = struct( ...
    'M_ary', M_ary, ...
    'EbN0', EbN0, ...
    'M', M, ...
    'T', T, ...
    'Fc', Fc, ...
    'vel', vel, ...
    'shape', shape, ...
    'alpha', alpha, ...
    'Q', Q...
    );

% Serialize to JSON for DB
paramsJSON  = jsonencode(params_inst);
paramHash = string(DataHash(paramsJSON,'SHA-256'));

% Load data from DB
sqlquery_results = sprintf("SELECT * FROM sim_results WHERE param_hash = '%s' AND system_name = '%s'", ...
    paramHash,systemName);
sim_result = fetch(conn, sqlquery_results);

if ~isempty(sim_result)
    % Find new frame count to simulate
    if sim_result.frames_simulated < num_frames
        new_frames = num_frames - sim_result.frames_simulated;
        run_flag = true;
    else
        run_flag = false;
    end
else
    % Simulate given frame count
    new_frames = num_frames;
    run_flag = true;
end

% % Simulation Settings
% data_folder = 'Data/OFDM';
% if ~exist(data_folder, 'dir')
%     mkdir(data_folder);
% end
%
% % Make sim_result file and filename
% sim_result.parameters = struct( ...
%     'M_ary', M_ary, ...
%     'EbN0', EbN0, ...
%     'M', M, ...
%     'T', T, ...
%     'Fc', Fc, ...
%     'vel', vel...
%     );
% param_fields = fieldnames(orderfields(sim_result.parameters));
% param_vals = struct2cell(orderfields(sim_result.parameters));
% key_parts = strcat(param_fields, '_', string(param_vals));
% hash_key = strjoin(key_parts, '__');
% data_filename = fullfile(data_folder,char('sim_OFDM_' + hash_key + '.mat'));
% index_filename = fullfile(data_folder,'results_index_OFDM.mat');
%
% % Add data to sim file
% if isfile(data_filename)
%
%     % Load file
%     load(data_filename, 'sim_result');
%
%     % Find new frame count to simulate
%     if sim_result.frames_simulated < num_frames
%         new_frames = num_frames - sim_result.frames_simulated;
%         run_flag = true;
%     else
%         run_flag = false;
%     end
%
% else
%
%     % Simulate given frame count
%     new_frames = num_frames;
%     run_flag = true;
%
% end

% Run if new frames are needed
if run_flag

    % Input conversion
    obj = comms_obj;
    obj.Eb_N0_db = EbN0;
    obj.M = M_ary;
    obj.N_subcarriers = M;
    obj.select_mod = "MPSK";
    F_d = (vel * (1000/3600))*Fc / physconst('LightSpeed');
    obj.FdT0 = F_d * T;
    obj.u = 1;

    % Convert to local variables
    N0 = obj.N0;
    S = obj.S;
    N = obj.N_subcarriers;
    L_efct = obj.L_efct;
    u = obj.u;
    FdT0 = obj.FdT0;
    F_d = obj.Fd;
    R_h_half = obj.R_h_half;
    R_p_half = obj.R_p_half;
    R_p = obj.R_p;
    T2 = obj.T2;
    syms_per_f = N;

    % Model F_N and F_uN
    F_N = generate_DFT(1,N);
    F_uN = generate_DFT(u,N);

    % Create pseudo-inverse Phi and noise-whitening matrix D if not
    % quasi-static
    if FdT0 ~= 0
        % Model pseudo-inverse of R_w, Phi
        [V_p,Omega_p,~] = svd(R_p);

        % Model "half" of Phi, D
        D = sqrt(1/u) * sqrt(inv(Omega_p)) * V_p' * F_uN';
    end

    % Create progress bar
    update_vals = floor(new_frames*(linspace(0.01,1,34)));
    fprintf("|           Progress:            |\n")

    % Reset bit errors for each SNR
    bit_errors = zeros(new_frames,1);
    sym_errors = zeros(new_frames,1);
    frm_errors = zeros(new_frames,1);

    % Initialize vectors
    for frame = 1:new_frames

        % Increment progress bar
        leq_locs = find(update_vals <= frame);
        inst_locs = find(update_vals <= frame & update_vals > frame - 1);
        if frame == 1
            for i = 1:length(leq_locs)
                fprintf("x")
            end
        else
            for i = 1:length(inst_locs)
                fprintf("x")
            end
        end

        % TX
        [TX_data,s] = generate_data(S,N);

        % Channel
        h_T = generate_fading(N*u,L_efct,F_d,T2);
        g_T = R_h_half * h_T;

        % Noise
        n1 = sqrt(1 / 2) * (randn(N*u,1) + 1j*randn(N*u,1));
        z_T = sqrt(N0) * R_p_half * n1;

        % Set up channel matrices
        G = zeros(u*N);
        G(1:size(g_T,1),1:size(g_T,2)) = g_T;
        for k = 1:u*N
            G(:,k) = circshift(G(:,k),k-1);
        end
        G_T = zeros(u*N,N);
        for k = 1:N
            G_T(:,k) = G(:,1+(k-1)*u);
        end

        % Make frequency domain components
        G_F = F_uN * G_T * F_N';
        z_F = F_uN * z_T;

        % Detection
        s_hat = zeros(N,1);
        RX_data = zeros(N,1);
        if FdT0 ~= 0 % Non-quasi-static fading
            % Create Block variables
            H_bar = D * G_F;
            w_bar = D * z_F;
            r_bar = H_bar * s + w_bar;
            s_hat = bdfe(r_bar.',H_bar,N0,S);
            RX_data = demodulator(s_hat,S);

        else % Quasi-static fading
            % Set up g_T for equalizer
            g_T_eq = zeros(u*N,1);
            g_T_eq(1:L_efct) = g_T(:,1);

            for k = 0:N-1
                % Find u-rows of normalized uN-point DFT
                F_uN_k = generate_DFT_k(u,N,k);

                % Find transformations unique for each k
                h_k = sqrt(N) * F_uN_k * g_T_eq;
                w_k = zeros(u,1);
                for m1 = 0:u-1
                    w_k(m1+1) = z_F(m1*N + k + 1);
                end
                r_k = h_k * s(k+1) + w_k;
                R_w = N0 * F_uN_k * R_p  * F_uN_k';

                % Model pseudoinverse matrix Phi_k and noise whitening matrix D
                [V_k,Omega_k] = eig(R_w / N0);
                Phi_k = V_k * (Omega_k \ V_k');

                % Create equalizer variables
                q_k = h_k' * Phi_k * h_k;
                phi_k = h_k' * Phi_k * r_k;

                % Decide most likely symbol
                arg_result = abs(phi_k - q_k*S).^2;
                [~,RX_data(k+1)] = min(arg_result);
                s_hat(k+1) = S(RX_data(k+1));
            end
        end

        % Compute number of errors in this frame and add to stack
        error_vec = RX_data ~= TX_data;
        bit_errors(frame) = log2(M_ary) * sum(error_vec(:));
        sym_errors(frame) = sum(error_vec(:));
        if sum(sym_errors(frame)) > 0
            frm_errors(frame) = 1;
        end

    end

    % Calculate BER, SER and FER
    BER = sum(bit_errors,"all") / (new_frames*syms_per_f*log2(M_ary));
    SER = sum(sym_errors,"all") / (new_frames*syms_per_f);
    FER = sum(frm_errors,"all") / (new_frames);

    % Set commands
    flag_sel = 0;
    sqlquery_flag = sprintf("SELECT * FROM system_flags WHERE id = '%d'", ...
        flag_sel);
    sqlquery_flagset = sprintf("UPDATE system_flags SET flag_value=%d WHERE id=%d", ...
        true, flag_sel);
    sqlquery_flagunset = sprintf("UPDATE system_flags SET flag_value=%d WHERE id=%d", ...
        false, flag_sel);

    % Write to database
    need_to_write = true;
    while need_to_write

        % Check system usage flag
        flag_row = fetch(conn, sqlquery_flag);
        flag_val = flag_row.flag_value;

        if ~flag_val

            % Set system usage flag 0
            exec(conn, sqlquery_flagset);

            % Load from DB again
            sim_result = fetch(conn, sqlquery_results);

            if ~isempty(sim_result) % Overwrite row in DB

                % Get new frame count
                N_total = sim_result.frames_simulated + new_frames;

                % Add new data to stack
                old_metrics = jsondecode(sim_result.metrics{1});
                metrics.BER = ...
                    (old_metrics.BER * sim_result.frames_simulated + BER * new_frames) / N_total;
                metrics.SER = ...
                    (old_metrics.SER * sim_result.frames_simulated + SER * new_frames) / N_total;
                metrics.FER = ...
                    (old_metrics.FER * sim_result.frames_simulated + FER * new_frames) / N_total;
                metricsJSON = jsonencode(metrics);

                % Format and execute SQL update string
                sqlupdate = sprintf("UPDATE sim_results SET metrics = '%s', frames_simulated = %d WHERE system_name = '%s' AND param_hash = '%s'", ...
                    metricsJSON, ...
                    N_total, ...
                    systemName, ...
                    paramHash);

                exec(conn, sqlupdate);

            else % Make new row in DB

                % Get new frame count
                N_total = new_frames;

                metrics.BER = BER;
                metrics.SER = SER;
                metrics.FER = FER;
                metricsJSON = jsonencode(metrics);

                sim_result_new = table( ...
                    string(systemName), ...
                    string(paramHash), ...
                    string(paramsJSON), ...
                    string(metricsJSON), ...
                    N_total, ...
                    'VariableNames', {'system_name', 'param_hash', 'parameters', 'metrics', 'frames_simulated'} );

                sqlwrite(conn,'sim_results',sim_result_new);

            end

            % Unset system usage flag 0
            exec(conn, sqlquery_flagunset);

            % No longer need to write to database
            need_to_write = false;

        else

            % Wait a random time between 1 and 5 seconds
            waitTime = 1 + (5 - 1) * rand();
            pause(waitTime);

        end

    end

%     % Add data to sim file
%     sim_result.timestamp = datetime('now');
%     if isfile(data_filename)
% 
%         % Load file
%         load(data_filename, 'sim_result');
% 
%         % Get new frame count
%         N_total = sim_result.frames_simulated + new_frames;
% 
%         % Add new data to stack
%         sim_result.metrics.BER = ...
%             (sim_result.metrics.BER * sim_result.frames_simulated + BER * new_frames) / N_total;
%         sim_result.metrics.SER = ...
%             (sim_result.metrics.SER * sim_result.frames_simulated + SER * new_frames) / N_total;
%         sim_result.metrics.FER = ...
%             (sim_result.metrics.FER * sim_result.frames_simulated + FER * new_frames) / N_total;
% 
%         % Update number of frames
%         sim_result.frames_simulated = N_total;
% 
%     else
% 
%         % Add new data to stack
%         sim_result.metrics = struct( ...
%             'BER', BER, ...
%             'SER', SER, ...
%             'FER', FER ...
%             );
% 
%         % Update number of frames
%         sim_result.frames_simulated = num_frames;
% 
%     end
% 
%     % Save data file
%     save(data_filename,"sim_result")
% 
%     % Load indexing file
%     if isfile(index_filename)
%         load(index_filename,'results_index');
%     else
%         results_index = table( ...
%             'Size',[0 6], ...
%             'VariableTypes',{'double','double','double','double','string','datetime'}, ...
%             'VariableNames',{'FramesSimulated','BER','SER','FER','HashKey','LastUpdated'});
%     end
% 
%     % Insert or update the row in the index
%     row_idx = find(results_index.HashKey == hash_key,1,'first');
%     if isempty(row_idx)
%         newRow = {sim_result.frames_simulated, ...
%             sim_result.metrics.BER, ...
%             sim_result.metrics.SER, ...
%             sim_result.metrics.FER, ...
%             hash_key, datetime('now')};
%         results_index = [results_index; newRow];
%     else
%         results_index{row_idx,{'FramesSimulated','BER','SER','FER'}} = ...
%             [sim_result.frames_simulated, ...
%             sim_result.metrics.BER, ...
%             sim_result.metrics.SER, ...
%             sim_result.metrics.FER];
%         results_index.LastUpdated(row_idx) = datetime('now');
%     end
% 
%     % Save indexing file
%     save(index_filename,'results_index');

end

end

%% Functions needed inside this function
function fft_matrix = generate_DFT(u,N)
% This is a function for generating the normalized N-point Discrete Fourier
% Transform matrix
omega = exp(-1j * 2*pi / (u*N));
fft_matrix = zeros(u*N);
for m1 = 1:1:u*N
    for n1 = 1:1:u*N
        fft_matrix(m1,n1) = omega^((m1-1) * (n1-1));
    end
end

fft_matrix = fft_matrix / sqrt(u*N);

end

function F_N_k = generate_DFT_k(u,N,k)
% This is a function for generating u-rows from the normalized uN-point
% Discrete Fourier Transform matrix
F_N_k = zeros(u,u*N);
for m2 = 1:1:u
    for n2 = 1:1:u*N
        F_N_k(m2,n2) = (1 / sqrt(u*N)) * exp(-1j * 2*pi * ((m2-1)*N + k) * (n2-1) / (u*N));
    end
end
end