%% MAIN_simulator_v2.m
% This file tests the BER/SER/FER for a few wireless communications
% systems (supported: OTFS, ODDM, TODDM), with settings specified in each
% profile. Data is saved in a MySQL server so a password is required.
%
% Coded 6/9/2025, JRW
clc; clear;

% Settings
frames_per_iter = 10;
create_database_tables = false;
save_data.priority = "local";
save_data.save_mysql = false;
save_data.save_excel = true;

% Set paths and data
addpath(fullfile(pwd, '..\Common Functions'));
addpath(fullfile(pwd, 'Functions'));
javaaddpath('..\mysql-connector-j-8.4.0.jar');
dbname     = 'comm_database';
table_name = "sim_results_v2";
save_data.excel_folder = 'Data';
save_data.excel_name = table_name;
save_data.excel_path = fullfile(save_data.excel_folder,save_data.excel_name + ".xlsx");

% Profile names (for profile selector)
profile_names = {
    "BER - System Comparison (max(M) = 128, N = 16)"
    "BER - M/U Comparison (max(M) = 128, N = 16)"
    "Thr - Subcarrier Spacing Comparison (15-480kHz)"
    "BER - Hexgrid (2048 symbols, Eb/N0=24dB)"
    "BER - Hexgrid (4096 symbols, Eb/N0=24dB)"
    };
[profile_sel,num_frames] = profile_select(profile_names,true);

% Set number of frames per iteration
if num_frames <= 0
    % Settings
    render_figure = true;
    save_sel = true;
    render_delay = 5;
    skip_simulations = true;
else
    % Settings
    skip_simulations = false;

    % Figure settings
    [render_figure,save_sel,render_delay] = figure_settings();
end

% Data select
switch profile_sel
    case 1
        % PROFILE 1
        vis_type = "figure";
        data_type = "BER";
        primary_var = "EbN0";
        primary_vals = 9:3:24;
        title_vars = ["T","N","vel"];
        default_parameters = struct(...
            'system_name', "TODDM",...
            'M_ary', 4, ...
            'EbN0', 18, ...
            'M', 16, ...
            'N', 16, ...
            'U', 8, ...
            'T', 1 / 15000, ...
            'Fc', 4e9, ...
            'vel', 500, ...
            'shape', "rrc", ...
            'alpha', 0.4, ...
            'Q', 8);
        configs = {
            struct('system_name', "OTFS", 'M', 128, 'U', 1)
            struct('system_name', "ODDM", 'M', 128, 'U', 1)
            struct('system_name', "TODDM", 'M', 128, 'U', 1)
            struct('system_name', "TODDM", 'M', 64, 'U', 2)
            struct('system_name', "TODDM", 'M', 32, 'U', 4)
            };
        legend_vec = {
            "OTFS, M = 128"
            "ODDM, M = 128"
            "TODDM, M = 128, U = 1"
            "TODDM, M = 64, U = 2"
            "TODDM, M = 32, U = 4"
            };
        line_styles = {
            "-*"
            "-o"
            "-^"
            "--^"
            "-.^"
            };
        line_colors = {...
            "#FF0000"
            "#0F62FE"
            "#000000"
            "#000000"
            "#000000"
            };

    case 2
        % PROFILE 2
        vis_type = "figure";
        data_type = "BER";
        primary_var = "EbN0";
        primary_vals = 9:3:24;
        title_vars = ["N","alpha","Q"];
        default_parameters = struct(...
            'system_name', "TODDM",...
            'M_ary', 4, ...
            'EbN0', 18, ...
            'M', 1024, ...
            'N', 16, ...
            'U', 1, ...
            'T', 1 / 15000, ...
            'Fc', 4e9, ...
            'vel', 500, ...
            'shape', "rrc", ...
            'alpha', 1, ...
            'Q', 2);
        configs = {
            struct('M', 128, 'U', 1)
            struct('M', 64, 'U', 2)
            struct('M', 32, 'U', 4)
            struct('M', 16, 'U', 8)
            };
        legend_vec = {
            "M=128, U=1"
            "M=64, U=2"
            "M=32, U=4"
            "M=16, U=8"
            };
        line_styles = {
            "-o"
            "-+"
            "-*"
            "-x"
            };
        line_colors = {
            "#FF0000"
            "#0F62FE"
            "#24A249"
            "#4C00C6"
            };

    case 3
        % PROFILE 3
        vis_type = "figure";
        data_type = "Thr";
        primary_var = "EbN0";
        primary_vals = 6:3:21;
        title_vars = ["M","N","U","alpha"];
        default_parameters = struct(...
            'system_name', "TODDM",...
            'M_ary', 4, ...
            'EbN0', 18, ...
            'M', 64, ...
            'N', 16, ...
            'U', 2, ...
            'T', 1 / 15000, ...
            'Fc', 4e9, ...
            'vel', 500, ...
            'shape', "rrc", ...
            'alpha', 0.4, ...
            'Q', 8);
        configs = {
            struct('T',1/15000)
            struct('T',1/30000)
            struct('T',1/60000)
            struct('T',1/120000)
            struct('T',1/240000)
            };
        legend_vec = {
            "\Deltaf = 15kHz"
            "\Deltaf = 30kHz"
            "\Deltaf = 60kHz"
            "\Deltaf = 120kHz"
            "\Deltaf = 240kHz"
            };
        line_styles = {
            "-o"
            "-+"
            "-*"
            "-x"
            "-^"
            };
        line_colors = {
            "#FF0000"
            "#0F62FE"
            "#24A249"
            "#4C00C6"
            "#00CAFF"
            };

    case 4
        % PROFILE 4
        vis_type = "hexgrid";
        data_type = "BER";
        primary_var = "EbN0";
        primary_vals = 24;
        title_vars = ["T","N","vel"];
        default_parameters = struct(...
            'system_name', "TODDM",...
            'M_ary', 4, ...
            'EbN0', 18, ...
            'M', 1024, ...
            'N', 32, ...
            'U', 1, ...
            'T', 1 / 15000, ...
            'Fc', 4e9, ...
            'vel', 500, ...
            'shape', "rrc", ...
            'alpha', 0.4, ...
            'Q', 8);
        % configs = gen_hexgrid_configs(2048, 32, 2, 1);
        configs = gen_hexgrid_configs(1024, 32, 2, 1);
        legend_vec = {};
        line_styles = {};
        line_colors = {};

    case 5
        % PROFILE 5
        vis_type = "hexgrid";
        data_type = "BER";
        primary_var = "EbN0";
        primary_vals = 24;
        title_vars = ["T","N","vel"];
        default_parameters = struct(...
            'system_name', "TODDM",...
            'M_ary', 4, ...
            'EbN0', 18, ...
            'M', 1024, ...
            'N', 32, ...
            'U', 1, ...
            'T', 1 / 15000, ...
            'Fc', 4e9, ...
            'vel', 500, ...
            'shape', "rrc", ...
            'alpha', 0.4, ...
            'Q', 8);
        configs = gen_hexgrid_configs(4096, 32, 2, 1);
        legend_vec = {};
        line_styles = {};
        line_colors = {};
end



% Find function files, get parameter list, modify sim data as needed
var_names = fieldnames(default_parameters);
prvr_len = length(primary_vals);
conf_len = length(configs);

% Progress tracking setup
num_iters = ceil(num_frames / frames_per_iter);
dq = parallel.pool.DataQueue;
completed_tasks = 0;
total_tasks = prvr_len*conf_len*num_iters;

% Callback to update progress
afterEach(dq, @updateProgressBar);

% Set up connection to MySQL server
if save_data.save_excel
    conn_local = mysql_login(dbname);
else
    conn_local = [];
end
if create_database_tables
    exec(conn, fileread('comm_database_sim_results.sql'));
    exec(conn, fileread('comm_database_system_flags.sql'));
end

% Ensure the folder exists
if ~isfolder(save_data.excel_folder)
    mkdir(save_data.excel_folder);
end

% Make parameters
system_names = cell(prvr_len,conf_len);
params_cell = cell(prvr_len,conf_len);
prior_frames = zeros(length(primary_vals),length(configs));
for primvar_sel = 1:prvr_len

    % Set primary variable
    primvar_val = primary_vals(primvar_sel);

    % Go through each settings profile
    for sel = 1:conf_len

        % Increment count variable
        iter_idx = (primvar_sel-1)*prvr_len + sel;

        % set all other variables
        var_vals = cell(length(var_names),1);
        for var_sel = 1:length(var_names)
            if isfield(configs{sel}, var_names{var_sel}) % Set secondary variables
                var_vals{var_sel} = configs{sel}.(cell2mat(var_names(var_sel)));
            elseif string(var_names{var_sel}) ~= primary_var % Set generic variables
                var_vals{var_sel} = default_parameters.(cell2mat(var_names(var_sel)));
            elseif string(var_names{var_sel}) == primary_var % Set primary variable variables
                var_vals{var_sel} = primvar_val;
            end
        end

        % Make parameters
        fields = fieldnames(default_parameters);
        parameters = cell2struct(var_vals, fields);
        system_names{primvar_sel,sel} = var_vals{1};

        % Remove unnecessary variables to get correct hash
        if system_names{primvar_sel,sel} == "ODDM"
            parameters = rmfield(parameters, 'U');
        elseif system_names{primvar_sel,sel} == "OTFS"
            parameters = rmfield(parameters, 'U');
        end

        % Add parameters to stack
        params_cell{primvar_sel,sel} = parameters;

        % Load data from DB
        [~,paramHash] = jsonencode_sorted(parameters);
        switch save_data.priority
            case "mysql"
                if save_data.save_mysql
                    sim_result = mysql_load(conn_local,table_name,paramHash);
                elseif save_data.save_excel
                    try
                        T = readtable(save_data.excel_path, 'TextType', 'string');
                        sim_result = T(T.param_hash == paramHash, :);
                    catch
                        T = table;
                        sim_result = [];
                    end
                end
            case "local"
                if save_data.save_excel
                    try
                        T = readtable(save_data.excel_path, 'TextType', 'string');
                        sim_result = T(T.param_hash == paramHash, :);
                    catch
                        T = table;
                        sim_result = [];
                    end
                elseif save_data.save_mysql
                    sim_result = mysql_load(conn_local,table_name,paramHash);
                end
        end

        % Set prior frames
        if ~isempty(sim_result)
            prior_frames(primvar_sel,sel) = sim_result.frames_simulated;
        else
            prior_frames(primvar_sel,sel) = 0;
        end

    end
end

% Set up connection to MySQL server
conn_thrall = conn_local; % UNCOMMENT TO USE LINEAR SOLVER
% if isempty(gcp('nocreate')) % UNCOMMENT TO USE PARALLELIZATION
%     poolCluster = parcluster('local');
%     maxCores = poolCluster.NumWorkers;  % Get the max number of workers available
%     parpool(poolCluster, maxCores);     % Start a parallel pool with all available workers
% end
% parfevalOnAll(@() javaaddpath('..\..\mysql-connector-j-8.4.0.jar'), 0);

% SIM LOOP
if ~skip_simulations
    for iter = 1:num_iters

        % Set current frame goal
        if iter < num_iters
            current_frames = iter*frames_per_iter;
        else
            current_frames = num_frames;
        end

        % Go through each settings profile
        for primvar_sel = 1:prvr_len % UNCOMMENT TO USE LINEAR SOLVER
        % parfor primvar_sel = 1:prvr_len % UNCOMMENT TO USE PARALLELIZATION

            for sel = 1:conf_len

                % Select parameters
                parameters = params_cell{primvar_sel,sel};

                % Continue to simulate if need more frames
                if current_frames > prior_frames(primvar_sel,sel)

                    % % Set up connection to MySQL server
                    % conn_thrall = mysql_login(dbname); % UNCOMMENT TO USE PARALLELIZATION

                    % Notify main thread of progress
                    progress_bar_data = parameters;
                    progress_bar_data.system_name = system_names{primvar_sel,sel};
                    progress_bar_data.num_iters = num_iters;
                    progress_bar_data.iter = iter;
                    progress_bar_data.primvar_sel = primvar_sel;
                    progress_bar_data.sel = sel;
                    progress_bar_data.prvr_len = prvr_len;
                    progress_bar_data.conf_len = conf_len;
                    progress_bar_data.current_frames = current_frames;
                    progress_bar_data.num_frames = num_frames;
                    send(dq, progress_bar_data);

                    % Simulate under current settings
                    sim_save(save_data,conn_thrall,table_name,current_frames,parameters);

                    % % Close connection instance
                    % close(conn_thrall) % UNCOMMENT TO USE PARALLELIZATION

                end

            end
        
        end
    end
end

% Set up figure data
figure_data.data_type = data_type;
figure_data.primary_var = primary_var;
figure_data.primary_vals = primary_vals;
figure_data.title_vars = title_vars;
figure_data.legend_vec = legend_vec;
figure_data.line_styles = line_styles;
figure_data.line_colors = line_colors;
figure_data.save_sel = false;

% Generate figure
clc;
fprintf("Displaying results for profile %d:\n",profile_sel)
if render_figure
    figure_data.save_sel = save_sel;
    clf
    switch vis_type
        case "figure"
            gen_figure_v3(save_data,conn_local,table_name,default_parameters,system_names,configs,figure_data);
        case "hexgrid"
            gen_hex_layout(conn_local,default_parameters,system_names,configs,figure_data);
    end
end

% Close connection with database
close(conn_local);