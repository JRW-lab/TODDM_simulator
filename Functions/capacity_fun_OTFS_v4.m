function capacity_fun_OTFS_v4(conn,EbN0_range,num_frames,parameters)

% Make parameter structure
systemName = 'OTFS';
tableName = "cap_results";
fields = fieldnames(parameters);
for i = 1:numel(fields)
    if fields{i} ~= "alpha"
        eval([fields{i} ' = parameters.(fields{i});']);
    else
        eval([fields{i} '_inst = parameters.(fields{i});']);
    end
end

% Find all run flags
new_frames = zeros(length(EbN0_range),1);
run_flags = zeros(length(EbN0_range),1);
sqlquery_results = strings(length(EbN0_range),1);
paramHash = strings(length(EbN0_range),1);
paramsJSON = cell(length(EbN0_range),1);
for i = 1:length(EbN0_range)

    % Select EbN0
    EbN0 = EbN0_range(i);
    parameters.EbN0 = EbN0;

    % Serialize to JSON for DB
    paramsJSON{i}  = jsonencode(parameters);
    paramHash(i) = string(DataHash(paramsJSON{i},'SHA-256'));

    % Load data from DB
    sqlquery_results(i) = sprintf("SELECT * FROM %s WHERE param_hash = '%s' AND system_name = '%s'", ...
        tableName,paramHash(i),systemName);
    sim_result = fetch(conn, sqlquery_results(i));

    if ~isempty(sim_result)
        % Find new frame count to simulate
        if sim_result.frames_simulated < num_frames
            new_frames(i) = num_frames - sim_result.frames_simulated;
            run_flags(i) = true;
        else
            run_flags(i) = false;
        end
    else
        % Simulate given frame count
        new_frames(i) = num_frames;
        run_flags(i) = true;
    end

end

% Run if new frames are needed
if sum(run_flags) > 0

    % Function settings
    EsN0_range = 10.^(log2(M_ary) .* EbN0_range ./ 10);

    % Input conversion
    sys = obj_comms_OTFS;
    sys.M_ary = M_ary;
    sys.N_tsyms = N;
    sys.M_sbcars= M;
    sys.filter = shape;
    sys.q = Q;
    sys.sbcar_spacing = 1 / T;
    sys.Fc = Fc;
    sys.v_vel = vel;
    sys.rolloff = alpha_inst;

    % Equalizer settings
    ambig_res = 1001;
    % rng('default');

    % Inputs from system object
    q = sys.q;
    Ts = sys.Ts;
    Lp = sys.Lp;
    Ln = sys.Ln;
    if shape == "rect"
        q = 1;
        alpha_inst = 1;
    elseif shape == "sinc"
        alpha_inst = 1;
    end

    load_test = sys.ambig_vals;
    if isempty(load_test)
        % Set discrete ambiguity lookup table if not done already
        fprintf("Checking ambiguity table and noise covariance matrix...\n")
        if shape ~= "rect"
            if ~exist("Pre-rendered Lookup Tables\\OTFS DD Cross-Ambiguity Tables", 'dir')
                mkdir("Pre-rendered Lookup Tables\\OTFS DD Cross-Ambiguity Tables")
            end
            filename = sprintf("Pre-rendered Lookup Tables\\OTFS DD Cross-Ambiguity Tables\\ambig_discrete_M%d_N%d_T%d_Fc%d_vel%d_%s_alpha%.1f_q%d.mat",M,N,T,Fc,vel,shape,alpha_inst,q);
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
                fprintf("Generating discrete ambiguity table (one time process)...\n")
                fprintf("|           Progress:            |\n")
                for k = 1:length(ambig_t_range)
                    if ismember(k,update_vals)
                        fprintf("x")
                    end

                    for l = 1:length(ambig_f_range)
                        ambig_vals(k,l) = ambig_direct(ambig_t_range(k),ambig_f_range(l),Ts,shape,alpha_inst,q,ambig_res);
                    end
                end
                fprintf("\n")

                % Save to file
                save(filename,"ambig_t_range","ambig_f_range","ambig_vals");
            end

            fprintf("Complete!\n\n")
        else
            ambig_t_range = [];
            ambig_f_range = [];
            ambig_vals = [];
        end
        sys.ambig_t_range = ambig_t_range;
        sys.ambig_f_range = ambig_f_range;
        sys.ambig_vals = ambig_vals;
    end

    % Bring in previously loaded data
    ambig_t_range = sys.ambig_t_range;
    ambig_f_range = sys.ambig_f_range;
    ambig_vals = sys.ambig_vals;
    if shape ~= "rect"
        res_chn_tau = (ambig_t_range(2)-ambig_t_range(1));
        res_chn_v = ambig_f_range(2)-ambig_f_range(1);
    end

    % Create progress bar
    fprintf("Generating H matrices...\n")
    update_vals = floor(num_frames*(linspace(0.01,1,34)));
    fprintf("|           Progress:            |\n")

    C_cell = cellfun(@(x) zeros(x,1), num2cell(new_frames),"UniformOutput",false);
    for frame = 1:max(new_frames)

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

        % Generate channel information
        [chn_g,chn_tau,chn_v] = channel_generation(Fc,vel);

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

        % Using SVD method
        s = svd(H);
        s_squared = s.^2;

        % Loop through SNR values
        for j = 1:length(EbN0_range)
            if frame <= new_frames(j)

                % Calculate capacity
                % SVD method
                C_cell{j}(frame) = sum(log2(1 + EsN0_range(j) * s_squared));

            end
        end


    end

    % Get average of all iterations
    C = cellfun(@(x) mean(x), C_cell);

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

            for i = 1:length(EbN0_range)

                if new_frames(i) > 0

                    % Load from DB again
                    sim_result = fetch(conn, sqlquery_results(i));

                    if ~isempty(sim_result) % Overwrite row in DB

                        % Get new frame count
                        N_total = sim_result.frames_simulated + new_frames(i);

                        % Add new data to stack
                        old_metrics = jsondecode(sim_result.metrics{1});
                        metrics.C = ...
                            (old_metrics.C * sim_result.frames_simulated + C(i) * new_frames(i)) / N_total;
                        metricsJSON = jsonencode(metrics);

                        % Format and execute SQL update string
                        sqlupdate = sprintf("UPDATE %s SET metrics = '%s', frames_simulated = %d WHERE system_name = '%s' AND param_hash = '%s'", ...
                            tableName, ...
                            metricsJSON, ...
                            N_total, ...
                            systemName, ...
                            paramHash(i));

                        % Update data row
                        exec(conn, sqlupdate);

                    else % Make new row in DB

                        % Get new frame count
                        N_total = new_frames(i);

                        % Create new metrics
                        metrics.C = C(i);
                        metricsJSON = jsonencode(metrics);

                        % Create data row
                        sim_result_new = table( ...
                            string(systemName), ...
                            string(paramHash(i)), ...
                            string(paramsJSON{i}), ...
                            string(metricsJSON), ...
                            N_total, ...
                            'VariableNames', {'system_name', 'param_hash', 'parameters', 'metrics', 'frames_simulated'} );

                        % Add new data row
                        sqlwrite(conn,tableName,sim_result_new);

                    end

                end

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