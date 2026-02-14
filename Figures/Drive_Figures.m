% Main driver script for plotting based on analysis of ICESat-2 data.
% This script sets up default paths and generates figures for the manuscript.

clear
close all

%% Setup Paths
% Get the location of this script to determine relative paths

% Assuming structure:
%   repo_root/
%     Figures/Drive_Figures.m
%     Processing/
%     Data/
OPTS.code_folder = '/Users/chorvat/Code/compare_MIZ_width';

% Define Data Input Path
% Pointing to the output of drive_MIZ_width.m
IS2_data_folder = fullfile(OPTS.code_folder,'Data','IS2_Data');
SIC_data_folder = fullfile(OPTS.code_folder,'Data','SIC_Data');

IS2_load_str = fullfile(IS2_data_folder, 'AT_stats_SH_v6_all.mat');

%% Setup analysis options and paths.
% This will overwrite the pre-existing OPTS from the loaded data (if any).
% Reset these with local paths.

% Location of the analysis code.
OPTS.code_folder = OPTS.code_folder;

% Location of figure scripts and utilities.
OPTS.plot_folder = fullfile(OPTS.code_folder, 'Figures');
OPTS.plot_utils_folder = fullfile(OPTS.code_folder, 'Utilities', 'Figures');

% Add necessary paths
addpath(OPTS.plot_folder);
addpath(OPTS.plot_utils_folder);
addpath(genpath(OPTS.plot_utils_folder));

% Also add subfolders for figures
addpath(fullfile(OPTS.plot_folder, 'Figure_1_PM_Comp'));
addpath(fullfile(OPTS.plot_folder, 'Figures_IS2'));
addpath(fullfile(OPTS.plot_folder, 'Figure_SAR_Comp'));

% Add Utilities path if needed (might be needed for some utility functions)
addpath(fullfile(OPTS.code_folder, 'Utilities', 'Processing'));
addpath(fullfile(OPTS.code_folder, 'Utilities', 'Analysis'));

% Output folder for figures
OPTS.plot_save_str = fullfile(OPTS.code_folder, 'Figures', 'Output');

if ~exist(OPTS.plot_save_str, 'dir')
    mkdir(OPTS.plot_save_str);
end

OPTS.plot_save_str = '/Users/chorvat/Library/CloudStorage/Dropbox-Brown/Christopher Horvat/Apps/Overleaf/PM-SIC-JOG/Figures/'

% OPTS.plot_save_str = '/Users/chorvat/Brown Dropbox/Christopher Horvat/Apps/Overleaf/PM-SIC-JOG/Figures/'; 

fprintf('Saving figures to: %s \n', OPTS.plot_save_str);


%% ----------- Passive Microwave Related Figures ---------------


% Display SIC Settings
disp('-----')
disp('Loading Gridded PM-SIC data')
fprintf('Using PM-SIC data located at: %s \n', SIC_data_folder);
disp('-----')

disp('Calculating Statistics of PM-SIC Products');

calc_PM_stats;


%%
% Figure 1 is the base comparison of PM data
plot_PM_global_stats;

%%

plot_PM_bias_map;



%% Loading IS2 Data

% Display the configuration settings.
disp('-----')
disp('Loading MIZ-referenced IS2 data')
fprintf('Using IS2 data located at: %s \n', IS2_data_folder);
disp('-----')


load(IS2_load_str);

%% ----------- IS2 Related Figures ---------------

load_MIZ_waves;

%% Location fo all tracks

create_location_data;

%% Figure on bias between AMSR2-NT2 and CDR

create_composite_figure;

%% Examine differences in MIZ width

create_MIZ_width_panel;

%% Parameteric figure of LIF offset, AMSR2 offset

create_parametric_plots;

%% ----------- SAR Related Figures ---------------

%% The localized SAR comparison with a single IS2 track
% This seems to be standalone or requires separate data?
% Keeping the path add for now.
