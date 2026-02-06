% Main driver script for plotting based on analysis of ICESat-2 data.
% This script sets up default paths and generates figures for the manuscript.

clear
close all

%% Setup Paths
% Get the location of this script to determine relative paths
[current_path, ~, ~] = fileparts(mfilename('fullpath'));
% Assuming structure:
%   repo_root/
%     Figures/Drive_Figures.m
%     Processing/
%     Data/
root_path = fileparts(current_path);

% Define Data Input Path
% Pointing to the output of drive_MIZ_width.m
data_input_folder = fullfile(root_path, 'Data', 'Output');
load_str = fullfile(data_input_folder, 'AT_stats_SH_v6_all.mat');

% Display the configuration settings.
disp('-----')
disp('Loading MIZ-referenced IS2 data')
fprintf('Using IS2 data located at: %s \n', load_str);
disp('-----')

% Check if data file exists
if exist(load_str, 'file')
    load(load_str);
else
    warning('Data file not found: %s. Please run Processing/drive_MIZ_width.m first.', load_str);
    % In a real scenario we might stop here, but for now we continue setup
end

%% Setup analysis options and paths.
% This will overwrite the pre-existing OPTS from the loaded data (if any).
% Reset these with local paths.

% Location of the analysis code.
OPTS.code_folder = root_path;

% Location of figure scripts and utilities.
OPTS.plot_folder = fullfile(root_path, 'Figures');
OPTS.plot_utils_folder = fullfile(root_path, 'Utilities', 'Figures');

% Add necessary paths
addpath(OPTS.plot_folder);
% Also add subfolders for figures
addpath(fullfile(OPTS.plot_folder, 'Figure_1_PM_Comp'));
addpath(fullfile(OPTS.plot_folder, 'Figures_IS2'));
addpath(fullfile(OPTS.plot_folder, 'Figure_SAR_Comp'));

% Add Utilities path if needed (might be needed for some utility functions)
addpath(fullfile(root_path, 'Utilities', 'Processing'));
addpath(fullfile(root_path, 'Utilities', 'Analysis'));

% Output folder for figures
OPTS.plot_save_str = fullfile(root_path, 'Figures', 'Output');
if ~exist(OPTS.plot_save_str, 'dir')
    mkdir(OPTS.plot_save_str);
end
fprintf('Saving figures to: %s \n', OPTS.plot_save_str);

%% Now start with figures

if exist('MIZ_DATA', 'var') && exist('IS2_DATA', 'var')

    %% ----------- Passive Microwave Related Figures ---------------

    disp('Calculating Statistics of PM-SIC Products');
    if exist('calc_PM_stats', 'file')
        calc_PM_stats;
    else
        warning('calc_PM_stats.m not found on path.');
    end

    %%
    % Figure 1 is the base comparison of PM data
    if exist('plot_PM_global_stats', 'file')
        plot_PM_global_stats;
    end

    %%

    if exist('plot_PM_bias_map', 'file')
        plot_PM_bias_map;
    end

    %% ----------- IS2 Related Figures ---------------

    if exist('load_MIZ_waves', 'file')
        load_MIZ_waves;
    end

    %% Figure on bias between AMSR2-NT2 and CDR

    if exist('create_composite_figure', 'file')
        create_composite_figure;
    end

    %% Examine differences in MIZ width

    if exist('create_MIZ_width_panel', 'file')
        create_MIZ_width_panel;
    end

    %% Parameteric figure of LIF offset, AMSR2 offset

    if exist('create_parametric_plots', 'file')
        create_parametric_plots;
    end

else
    disp('Skipping figure generation steps because data was not loaded.');
end

%% ----------- SAR Related Figures ---------------

%% The localized SAR comparison with a single IS2 track
% This seems to be standalone or requires separate data?
% Keeping the path add for now.
