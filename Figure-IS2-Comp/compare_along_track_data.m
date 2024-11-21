% Updatable name of output. 

% OS_string = '/Users/chorvat/Dropbox-Brown/Christopher Horvat/';
OS_string = '/Users/chorvat/Brown Dropbox/Christopher Horvat/';

OPTS.output_str = [OS_string 'Research Projects/Active/Data/ICESat-2/PM-SIC-width/Along_Track_Statistics/AT_stats_SH_v6_all']; 
% OPTS.output_str = [OS_string 'Research Projects/Active/Data/ICESat-2/PM-SIC-width/Along_Track_Statistics/AT_stats_KM_v6_all']; 


%%
load(OPTS.output_str);

%%

cutoff_N = 1000; % need more than 100 segments for analysis to make sense
cutoff_MIZ = 1; % Need more than 1 MIZ segment for analysis. 

close all

addpath([OPTS.code_folder 'Utilities/Plotting'])
 
make_figure_MIZ_all;

%%
make_figure_MIZ_waves;

%%

make_figure_MIZ_somewaves;

make_figure_MIZ_nowaves;




%%
make_location_figure; 



