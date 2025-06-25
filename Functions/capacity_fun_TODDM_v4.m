function capacity_fun_TODDM_v4(conn,EbN0_range,num_frames,parameters)

% Make parameter structure
systemName = 'TODDM';
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
    res = 10;
    EsN0_range = 10.^(log2(M_ary) * EbN0_range ./ 10);

    % Render ambiguity tables
    [Ambig_Table.vals,Ambig_Table.t_range,Ambig_Table.f_range] = gen_DF_ambig_table(U,N,M,T,Fc,vel,shape,alpha_inst,Q,res,true);

    % Create progress bar
    fprintf("Generating H matrices...\n")
    update_vals = floor(max(new_frames)*(linspace(0.01,1,34)));
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

        % Generate H matrix and channel information
        H = gen_HDD_triorth(U,N,M,T,Q,Fc,vel,Ambig_Table);

        % Using SVD method
        s = svd(H);
        s_squared = abs(s).^2;

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