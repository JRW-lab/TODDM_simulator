function sim_fun_TODDM_v2(conn,num_frames,parameters)

% Make parameter structure
system_name = 'TODDM';
table_name = "sim_results";
fields = fieldnames(parameters);
for i = 1:numel(fields)
    eval([fields{i} ' = parameters.(fields{i});']);
end

% Serialize to JSON for DB
paramsJSON  = jsonencode(parameters);
paramHash = string(DataHash(paramsJSON,'SHA-256'));

% Load data from DB
sqlquery_results = sprintf("SELECT * FROM %s WHERE param_hash = '%s' AND system_name = '%s'", ...
    table_name,paramHash,system_name);
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

    % Cause error if settings are inconsistent with theory
    Ts = T / M;
    if 2*Q + floor(2510*10^(-9) / Ts) >= M
        error("Settings can not satisfy the ambiguity assumption.")
    end

    % Define parameters
    res = 10;
    Es = 1;
    Eb = Es / log2(M_ary);
    N_iters = 10;
    syms_per_f = M*N*U;
    N0 = Eb / (10^(EbN0 / 10)) * ((N+2)/N);

    % Add redundancy for rectangular and sinc pulses
    if shape == "rect"
        Q = 1;
        alpha = 1;
    elseif shape == "sinc"
        alpha = 1;
    end

    % Data setup
    if M_ary == 2
        % Define bit order
        bit_order = [0;1];

        % Define alphabet
        alphabet_set = linspace(1,M_ary,M_ary)';
        S = sqrt(Es) .* exp(-1j * 2*pi .* (alphabet_set) ./ M_ary);
    elseif M_ary == 4
        % Define bit order
        bit_order = [0,0;0,1;1,0;1,1];

        % Define alphabet
        S = zeros(4,1);
        S(1) = (sqrt(2)/2) + (1j*sqrt(2)/2);
        S(2) = (sqrt(2)/2) - (1j*sqrt(2)/2);
        S(3) = -(sqrt(2)/2) + (1j*sqrt(2)/2);
        S(4) = -(sqrt(2)/2) - (1j*sqrt(2)/2);
        S = sqrt(Es) .* S;
    end

    % Render ambiguity tables
    [Ambig_Table.vals,Ambig_Table.t_range,Ambig_Table.f_range] = gen_DF_ambig_table(U,N,M,T,Fc,vel,shape,alpha,Q,res,true);

    % % Create progress bar
    % update_vals = floor(new_frames*(linspace(0.01,1,34)));
    % fprintf("|           Progress:            |\n")

    % Reset bit errors for each SNR
    bit_errors = zeros(new_frames,syms_per_f*log2(M_ary));
    sym_errors = zeros(new_frames,syms_per_f);
    frm_errors = zeros(new_frames,1);
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
        [TX_bit,TX_sym,xDD] = gen_data(bit_order,S,syms_per_f);

        % Generate discrete channel matrix
        [HDD,L1,L2] = gen_HDD_triorth(U,N,M,T,Q,Fc,vel,Ambig_Table);

        % Generate noise
        zDD = sqrt(N0/2) * (randn(syms_per_f,1) + 1j*randn(syms_per_f,1));

        % Construct received signal
        yDD = HDD * xDD + zDD;

        % Equalize received signal
        x_hat = OTFS_pulse_equalizer_AWGN(yDD,HDD,N*U,M,L2,-L1,Es,N0,S,N_iters);

        % Hard detection for final x_hat
        dist = abs(x_hat.' - S).^2;
        [~,min_index] = min(dist);
        RX_sym = min_index.' - 1;

        % Convert final RX_sym to RX_bit
        RX_bit = bit_order(RX_sym+1,:);

        % Error calculation
        diff_bit = TX_bit ~= RX_bit;
        diff_sym = TX_sym ~= RX_sym;
        bit_errors(frame,:) = diff_bit(:).';
        sym_errors(frame,:) = diff_sym(:).';
        if sum(sym_errors(frame,:)) > 0
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
                sqlupdate = sprintf("UPDATE %s SET metrics = '%s', frames_simulated = %d WHERE system_name = '%s' AND param_hash = '%s'", ...
                    table_name, ...
                    metricsJSON, ...
                    N_total, ...
                    system_name, ...
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
                    string(system_name), ...
                    string(paramHash), ...
                    string(paramsJSON), ...
                    string(metricsJSON), ...
                    N_total, ...
                    'VariableNames', {'system_name', 'param_hash', 'parameters', 'metrics', 'frames_simulated'} );

                sqlwrite(conn,table_name,sim_result_new);

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