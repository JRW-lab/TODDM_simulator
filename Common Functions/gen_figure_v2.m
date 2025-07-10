function gen_figure_v2(conn,default_parameters,system_names,configs,figure_data)

% Import settings
data_type = figure_data.data_type;
primary_var = figure_data.primary_var;
primary_vals = figure_data.primary_vals;
title_vars = figure_data.title_vars;
legend_vec = figure_data.legend_vec;
line_styles = figure_data.line_styles;
line_colors = figure_data.line_colors;
save_sel = figure_data.save_sel;

% Figure settings
render_title = false;
figures_folder = 'Figures';
line_val = 2;
mark_val = 10;
font_val = 16;

% Select MySQL table
switch conn.DataSource
    case "comm_database"
        identifier = "system_name";
        table_name = "sim_results_v2";
        var_start = 1;
    case "med_database"
        identifier = "project_name";
        table_name = "lrm_results_v2";
        var_start = 1;
end

% SIM LOOP
var_names = fieldnames(default_parameters);
results_mat = zeros(length(primary_vals),length(configs));
for primvar_sel = 1:length(primary_vals)

    % Set primary variable
    primval_sel = primary_vals(primvar_sel);

    % Go through each settings profile
    for sel = 1:length(configs)

        % Define label of result
        if identifier == "system_name"
            result_label = system_names{primvar_sel,sel};
        else
            result_label = "BPDB";
        end

        % Populate parameters from configs or default_parameters
        params_inst = struct();
        for var_sel = var_start:length(var_names)
            var_name = var_names{var_sel};
            if isfield(configs{sel}, var_name)
                value = configs{sel}.(var_name);
            elseif ~strcmp(var_name, primary_var)
                value = default_parameters.(var_name);
            else
                value = primval_sel;
            end
            params_inst.(var_name) = value;
        end

        % Edit parameters file if comm sim
        if identifier == "system_name"
            switch result_label
                case "TODDM"

                    % Make number of symbols per frame
                    syms_per_f = params_inst.M*params_inst.N*params_inst.U;
                    U = params_inst.U;

                case "ODDM"

                    % Correct parameters list, make number of symbols per frame
                    params_inst = rmfield(params_inst, 'U');
                    syms_per_f = params_inst.M*params_inst.N;
                    U = 1;

                case "OTFS"

                    % Correct parameters list, make number of symbols per frame
                    params_inst = rmfield(params_inst, 'U');
                    syms_per_f = params_inst.M*params_inst.N;
                    U = 1;

                case "OFDM"

                    % Correct parameters list, make number of symbols per frame
                    params_inst = rmfield(params_inst, 'U');
                    params_inst = rmfield(params_inst, 'N');
                    syms_per_f = params_inst.M;
                    U = 1;

            end
        end

        % Serialize to JSON for DB
        [~,paramHash] = jsonencode_sorted(params_inst);

        % Load data from DB
        switch conn.DataSource
            case "comm_database"
                DB_data = mysql_load(conn,table_name,paramHash);
            case "med_database"
                sqlquery = sprintf("SELECT * FROM %s WHERE param_hash = '%s' AND %s = '%s'", ...
                    table_name, paramHash, identifier, result_label);
                DB_data = fetch(conn, sqlquery);
        end

        % Select data to extract
        results_inst = jsondecode(DB_data.metrics{1});
        switch data_type
            case "BER"
                results_val = results_inst.BER;
            case "SER"
                results_val = results_inst.SER;
            case "FER"
                results_val = results_inst.FER;
            case "Thr"
                frame_duration = params_inst.N * params_inst.T;
                bandwidth_hz = params_inst.M * U / params_inst.T;
                results_val = (log2(params_inst.M_ary) * syms_per_f * (1 - results_inst.FER)) / (frame_duration * bandwidth_hz);
            case "Cap"
                results_val = results_inst.C / 1000;
            case "accy"
                level_view = figure_data.level_view;
                results_val = results_inst.(level_view).(data_type);
        end

        % Add result to stack
        if isempty(results_val)
            results_mat(primvar_sel,sel) = NaN;
        else
            results_mat(primvar_sel,sel) = results_val;
        end

    end
end

% Create folders if they don't exist
subfolder = fullfile(figures_folder, ['/' char(data_type)]);
subsubfolder = fullfile(subfolder, ['/' char(primary_var)]);
if ~exist(figures_folder, 'dir')
    mkdir(figures_folder);
end
if ~exist(subfolder, 'dir')
    mkdir(subfolder);
end
if ~exist(subsubfolder, 'dir')
    mkdir(subsubfolder);
end

% Change label depending on range parameter
switch primary_var
    case "EbN0"
        xlabel_name = "E_b/N_0";
    case "vel"
        xlabel_name = "Vehicular velocity (km/hr)";
    case "T"
        primary_vals = primary_vals * 10^6;
        xlabel_name = "T (us)";
    case "frequency_limit"
        xlabel_name = "Cutoff Frequency (Hz)";
    otherwise
        xlabel_name = primary_var;
end

% Display figure
figure(1)
hold on
for i = 1:size(results_mat,2)
    plot(primary_vals,results_mat(:,i), ...
        line_styles{i}, ...
        Color=line_colors{i}, ...
        LineWidth=line_val, ...
        MarkerSize=mark_val)
end
if identifier == "system_name"
    if data_type == "Cap"
        ylim_vec = [0 ceil(max(results_mat,[],"all")/10)*10];
        % set(gca, 'YScale', 'log');
        ylabel("Ergodic Capacity (kbps/Hz)")
        loc = "northwest";
        % loc = "southeast";
    elseif data_type == "Thr"
        ylim_vec = [0 ceil(max(results_mat,[],"all"))];
        ylabel("Throughput (bps/Hz)")
        loc = "southeast";
    else
        ylim_vec = [1e-6 1e-1];
        set(gca, 'YScale', 'log');
        ylabel(data_type)
        loc = "northeast";
    end
else
    if default_parameters.dataset == "Human"
        ylim_vec = [0.9 1];
        % ylim_vec = [0.7 1];
    else
        ylim_vec = [0.6 1];
    end
    ylabel("Model Accuracy")
    % loc = "southwest";
    loc = "southeast";
end
grid on
ylim(ylim_vec)
xlabel(xlabel_name)
xlim([min(primary_vals) max(primary_vals)])
xticks(primary_vals)
legend(legend_vec,Location=loc);
set(gca, 'FontSize', font_val);

if render_title
    % Make title combo
    title_vec = "";
    for i = 1:length(title_vars)
        if i > 1
            title_vec = sprintf("%s, ",title_vec);
        end

        switch title_vars(i)
            case "T"
                title_val = default_parameters.T / 1e-6;
                title_vec = sprintf("%s%s = %.2f",title_vec,title_vars(i),title_val);
            otherwise
                title_val = eval("default_parameters." + title_vars(i));
                title_vec = sprintf("%s%s = %d",title_vec,title_vars(i),title_val);
        end

        if title_vars(i) == "EbN0"
            title_vec = sprintf("%sdB",title_vec);
        elseif title_vars(i) == "vel"
            title_vec = sprintf("%s km/hr",title_vec);
        elseif title_vars(i) == "T"
            title_vec = sprintf("%s us",title_vec);
        end

    end
    title(title_vec)
end

% Save figure
timestamp = datetime('now', 'Format', 'yyyyMMdd_HHmmss');
timestamp_str = char(timestamp);
figure_filename = fullfile(subsubfolder, "Figure_" + timestamp_str + ".png");
if save_sel
    saveas(figure(1), figure_filename);
end