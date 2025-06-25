function sim_fun_OTFS_v2(conn,num_frames,parameters)

% Make parameter structure
sim_type = 'OTFS';
fields = fieldnames(parameters);
for i = 1:numel(fields)
    eval([fields{i} ' = parameters.(fields{i});']);
end

% Serialize to JSON for DB
paramsJSON  = jsonencode(parameters);
paramHash = string(DataHash(paramsJSON,'SHA-256'));

% Load data from DB
sqlquery_results = sprintf("SELECT * FROM sim_results WHERE param_hash = '%s' AND system_name = '%s'", ...
    paramHash,sim_type);
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

% Run if new frames are needed
if run_flag

    % Input conversion
    sys = obj_comms_OTFS;
    sys.EbN0_db = EbN0;
    sys.M_ary = M_ary;
    sys.N_tsyms = N;
    sys.M_sbcars= M;
    sys.filter = shape;
    sys.q = Q;
    sys.sbcar_spacing = 1 / T;
    sys.Fc = Fc;
    sys.v_vel = vel;
    sys.rolloff = alpha;

    % Equalizer settings
    N_iters = 8;
    ambig_res = 1001;
    % rng('default');

    % Inputs from system object
    Es = sys.Es;
    S = sys.S;
    M_ary = sys.M_ary;
    N0 = sys.N0;
    M = sys.M_sbcars;
    N = sys.N_tsyms;
    Fc = sys.Fc;
    vel = sys.v_vel;
    if M_ary == 2
        bit_order = [0;1];
    elseif M_ary == 4
        bit_order = [0,0;0,1;1,0;1,1];
    end
    shape = sys.filter;
    alpha = sys.rolloff;
    q = sys.q;
    T = sys.T;
    Ts = sys.Ts;
    Lp = sys.Lp;
    Ln = sys.Ln;
    if shape == "rect"
        q = 1;
        alpha = 1;
    elseif shape == "sinc"
        alpha = 1;
    end

    load_test = sys.ambig_vals;
    if isempty(load_test)
        % % Set discrete ambiguity lookup table if not done already
        % fprintf("Checking ambiguity table and noise covariance matrix...\n")
        if shape ~= "rect"
            if ~exist("Pre-rendered Lookup Tables\\OTFS DD Cross-Ambiguity Tables", 'dir')
                mkdir("Pre-rendered Lookup Tables\\OTFS DD Cross-Ambiguity Tables")
            end
            filename = sprintf("Pre-rendered Lookup Tables\\OTFS DD Cross-Ambiguity Tables\\ambig_discrete_M%d_N%d_T%d_Fc%d_vel%d_%s_alpha%.1f_q%d.mat",M,N,T,Fc,vel,shape,alpha,q);
            if isfile(filename)
                loaded_file = load(filename);
                ambig_t_range = loaded_file.ambig_t_range;
                ambig_f_range = loaded_file.ambig_f_range;
                ambig_vals = loaded_file.ambig_vals;
            else
                ambig_t_lim = Ts*(q + ceil((2510*10^(-9))/Ts));
                ambig_f_lim = ((vel * (1000/3600))*Fc) / (physconst('LightSpeed'));
                ambig_t_range = linspace(-ambig_t_lim,ambig_t_lim,ambig_res);
                ambig_f_range = linspace(-ambig_f_lim,ambig_f_lim,ambig_res);
                ambig_vals = zeros(ambig_res);
                update_vals = floor(length(ambig_t_range)*(linspace(.01,.99,34)));
                % fprintf("Generating discrete ambiguity table (one time process)...\n")
                % fprintf("|           Progress:            |\n")
                for k = 1:length(ambig_t_range)
                    % if ismember(k,update_vals)
                    %     fprintf("x")
                    % end

                    for l = 1:length(ambig_f_range)
                        ambig_vals(k,l) = ambig_direct(ambig_t_range(k),ambig_f_range(l),Ts,shape,alpha,q,ambig_res);
                    end
                end
                fprintf("\n")

                % Save to file
                save(filename,"ambig_t_range","ambig_f_range","ambig_vals");
            end

            % fprintf("Complete!\n\n")
        else
            ambig_t_range = [];
            ambig_f_range = [];
            ambig_vals = [];
        end
        sys.ambig_t_range = ambig_t_range;
        sys.ambig_f_range = ambig_f_range;
        sys.ambig_vals = ambig_vals;
    end

    load_test = sys.R;
    if isempty(load_test)
        % Set up noise covariance matrix if not done already
        if shape ~= "rect"
            if ~exist("Pre-rendered Lookup Tables\\OTFS Noise Covariance Matrices", 'dir')
                mkdir("Pre-rendered Lookup Tables\\OTFS Noise Covariance Matrices")
            end
            filename = sprintf("Pre-rendered Lookup Tables\\OTFS Noise Covariance Matrices\\Rzddt_T%d_N%d_M%d_q%d_%s_alpha%.1f.mat",T,N,M,q,shape,alpha);
            if isfile(filename)
                loaded_file = load(filename);
                R_half = loaded_file.R_half;
                R = loaded_file.R;
            else
                % fprintf("Generating noise covariance matrix (one time process)...\n")
                [R_half,R] = gen_Rzddt(T,N,M,q,ambig_res,shape,alpha);
                % fprintf("\n")

                % Save to file
                save(filename,"R_half","R");
            end
        else
            R_half = eye(N*M);
            R = R_half;
        end
        sys.R_half = R_half;
        sys.R = R;
    end

    % Needed variables and matrices setup
    syms_per_f = M*N;
    Gamma_MN = gen_Gamma_MN(M,N);
    F_N = gen_DFT(N);

    % Bring in previously loaded data
    ambig_t_range = sys.ambig_t_range;
    ambig_f_range = sys.ambig_f_range;
    ambig_vals = sys.ambig_vals;
    R = sys.R;
    R_half = sys.R_half;
    if shape ~= "rect"
        res_chn_tau = (ambig_t_range(2)-ambig_t_range(1));
        res_chn_v = ambig_f_range(2)-ambig_f_range(1);
    end

    % % Create progress bar
    % update_vals = floor(new_frames*(linspace(0.01,1,34)));
    % fprintf("|           Progress:            |\n")

    % Reset bit errors for each SNR
    bit_errors = zeros(new_frames,syms_per_f*log2(M_ary));
    sym_errors = zeros(new_frames,syms_per_f);
    frm_errors = zeros(new_frames,1);

    % Simulation loop
    for frame = 1:new_frames

        % % Increment progress bar
        % leq_locs = find(update_vals <= frame);
        % inst_locs = find(update_vals <= frame & update_vals > frame - 1);
        % if frame == 1
        %     for i = 1:length(leq_locs)
        %         fprintf("x")
        %     end
        % else
        %     for i = 1:length(inst_locs)
        %         fprintf("x")
        %     end
        % end

        % Generate data
        [TX_1,TX_2,x_DD] = gen_data(bit_order,S,syms_per_f);
        TX_bit = Gamma_MN' * TX_1;
        TX_sym = Gamma_MN' * TX_2;
        x_tilde = Gamma_MN' * x_DD;

        % Generate channel
        [chn_g,chn_tau,chn_v] = channel_generation(Fc,vel);
        if shape == "rect" % rectangular ambiguity is closed form
            % Create H Matrix
            H = gen_H(T,N,M,Lp,Ln,chn_g,chn_tau,chn_v,shape,alpha);
        else
            % Normalize tau and v to cohere with discrete ambig values
            chn_tau = round(chn_tau/res_chn_tau)*res_chn_tau;
            chn_v = round(chn_v/res_chn_v)*res_chn_v;

            % Find direct tap indices and tap values
            l = (Ln:Lp).';
            tap_t_range = (l*Ts - chn_tau) .* ones(Lp-Ln+1,length(chn_g),N*M);
            tap_f_range = (ones(Lp-Ln+1,1) .* chn_v) .* ones(Lp-Ln+1,length(chn_g),N*M);
            tap_t_range = round(tap_t_range ./ res_chn_tau) + ceil(length(ambig_t_range)/2);
            tap_f_range = round(tap_f_range ./ res_chn_v) + ceil(length(ambig_f_range)/2);
            tap_t_range(tap_t_range < 1) = 1;
            tap_t_range(tap_t_range > 1001) = 1001;
            tap_f_range(tap_f_range < 1) = 1;
            tap_f_range(tap_f_range > 1001) = 1001;

            % Create H Matrix
            H = gen_H_direct(T,N,M,Lp,Ln,chn_g,chn_v,ambig_vals,tap_t_range,tap_f_range);
        end
        H_tilde = gen_H_tilde(H,M,N,Lp,Ln,F_N);

        % Generate noise
        z_tilde = sqrt(N0/2) * R_half * (randn(syms_per_f,1) + 1j*randn(syms_per_f,1));

        % Create receive vector
        y_tilde = H_tilde * x_tilde + z_tilde;

        % Iterative Detector - JRW
        x_hat = OTFS_pulse_equalizer_v2(y_tilde,H_tilde,N,M,Lp,Ln,Es,N0,S,N_iters,R);

        % Hard detection for final x_hat
        dist = abs(x_hat.' - S).^2;
        [~,min_index] = min(dist);
        RX_sym = min_index.' - 1;

        % Convert final RX_sym to RX_bit
        RX_bit = bit_order(RX_sym+1,:);

        % Error calculation
        bit_error_vec = TX_bit ~= RX_bit;
        sym_error_vec = TX_sym ~= RX_sym;
        bit_errors(frame) = sum(bit_error_vec,"all");
        sym_errors(frame) = sum(sym_error_vec,"all");
        if sum(bit_error_vec,"all") > 0
            frm_errors(frame) = 1;
        else
            frm_errors(frame) = 0;
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
                    sim_type, ...
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
                    string(sim_type), ...
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

end
