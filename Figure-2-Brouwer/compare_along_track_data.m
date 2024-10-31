% Updatable name of output. 


% OPTS.output_str = '/Users/chorvat/Brown Dropbox/Christopher Horvat/Research Projects/Active/Data/ICESat-2/PM-SIC-width/Along_Track_Statistics/AT_stats_SH_v6_all'; 
OPTS.output_str = '/Users/chorvat/Brown Dropbox/Christopher Horvat/Research Projects/Active/Data/ICESat-2/PM-SIC-width/Along_Track_Statistics/AT_stats_KM_v6'; 

%%
load(OPTS.output_str);
%%
close all

addpath([OPTS.code_folder 'Utilities/Plotting'])

% make_figure_MIZ_somewaves(MIZ_DATA,IS2_DATA)
make_figure_MIZ_all;

%%
make_figure_MIZ_waves;

%%

make_figure_MIZ_nowaves;

%% Now evaluate the MIZ width. 

make_figure_MIZ_width; 

