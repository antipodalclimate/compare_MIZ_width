% Main driver script for plotting based on analysis of ICESat-2 data.
% This script sets up default paths

clear

% Define the base path for data access. Update this to your local machine.
% OS_string = '/Users/chorvat/Library/CloudStorage/Dropbox-Brown/Christopher Horvat/';


load_str = '~/Code/compare_MIZ_width/Processing/Data/AT_stats_SH_v6_all';


% Display the configuration settings.
disp('-----')
disp('Loading MIZ-referenced IS2 data')
fprintf('Using IS2 data located at: %s \n',load_str);
disp('-----')

load(load_str)


%% Setup analysis options and paths. This will overwrite the pre-existing
% Things that are in OPTS from the load data. Reset these with your local
% paths if you didn't create the data, or moved machines, etc. 

% Location of the analysis code.
OPTS.code_folder = '~/Code/compare_MIZ_width/';
% Location of utility functions.
OPTS.plot_folder = [OPTS.code_folder 'Figures/'];
OPTS.plot_utils_folder = [OPTS.code_folder 'Utilities/Figures/'];

addpath(OPTS.plot_folder);
addpath(OPTS.plot_utils_folder);

OPTS.plot_save_str = '~/Dropbox (Brown)/Apps/Overleaf/PM-SIC-JOG/Figures/';

%% Now start with figures


%% ----------- Passive Microwave Related Figures ---------------
% 
addpath([OPTS.plot_folder 'Figure_1_PM_Comp'])

disp('Calculating Statistics of PM-SIC Products');
calc_PM_stats; 

%% 
% Figure 1 is the base comparison of PM data
plot_PM_global_stats;

%%

plot_PM_bias_map; 

%% ----------- IS2 Related Figures ---------------

addpath([OPTS.plot_folder 'Figures_IS2'])

load_MIZ_waves; 



%% Figure on bias between AMSR2-NT2 and CDR

create_composite_figure; 

%% Examine differences in MIZ width 

create_MIZ_width_panel; 

%% Parameteric figure of LIF offset, AMSR2 offset

create_parametric_plots; 


%% ----------- SAR Related Figures ---------------

%% The localized SAR comparison with a single IS2 track

addpath([OPTS.plot_folder 'Figure_SAR_Comp'])
