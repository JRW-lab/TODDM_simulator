function gen_figure_v2(save_data,conn,table_name,hash_cell,configs,figure_data)

% Figure settings
figures_folder = 'Figures';
line_val = 2;
mark_val = 10;
font_val = 16;

% Load data from DB and set new frame count
switch save_data.priority
    case "mysql"
        if save_data.save_mysql
            T = mysql_load(conn,table_name,"*");
        elseif save_data.save_excel
            try
                T = readtable(save_data.excel_path, 'TextType', 'string');
            catch
                T = table;
            end
        end
    case "local"
        if save_data.save_excel
            try
                T = readtable(save_data.excel_path, 'TextType', 'string');
            catch
                T = table;
            end
        elseif save_data.save_mysql
            T = mysql_load(conn,table_name,"*");
        end
end

% Import settings
data_type = figure_data.data_type;
primary_var = figure_data.primary_var;
primary_vals = figure_data.primary_vals;
legend_vec = figure_data.legend_vec;
line_styles = figure_data.line_styles;
line_colors = figure_data.line_colors;
save_sel = figure_data.save_sel;

% Load loop
results_mat = zeros(length(primary_vals),length(configs));
for primvar_sel = 1:size(hash_cell,1)

    % Go through each settings profile
    for sel = 1:size(hash_cell,2)

        % Load data from DB
        paramHash = hash_cell{primvar_sel,sel};
        sim_result = T(string(T.param_hash) == paramHash, :);

        % Select data to extract
        metrics_loaded = jsondecode(sim_result.metrics{1});
        results_val = metrics_loaded.(data_type);

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
if data_type == "Thr"
    ylim_vec = [0 ceil(max(results_mat,[],"all"))];
    ylabel("Throughput (bps/Hz)")
    loc = "southeast";
else
    ylim_vec = [1e-6 1e-1];
    set(gca, 'YScale', 'log');
    ylabel(data_type)
    loc = "northeast";
end
grid on
ylim(ylim_vec)
xlabel(xlabel_name)
xlim([min(primary_vals) max(primary_vals)])
xticks(primary_vals)
legend(legend_vec,Location=loc);
set(gca, 'FontSize', font_val);

% Save figure
timestamp = datetime('now', 'Format', 'yyyyMMdd_HHmmss');
timestamp_str = char(timestamp);
figure_filename = fullfile(subsubfolder, "Figure_" + timestamp_str + ".png");
if save_sel
    saveas(figure(1), figure_filename);
end