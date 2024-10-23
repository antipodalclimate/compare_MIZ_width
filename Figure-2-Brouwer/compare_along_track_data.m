% Updatable name of output. 


OPTS.output_str = '/Users/chorvat/Brown Dropbox/Christopher Horvat/Research Projects/Active/Data/ICESat-2/PM-SIC-width/Along_Track_Statistics/AT_stats_Brouwer_v6'; 

%%
load(OPTS.output_str);
%%
close all

addpath([OPTS.code_folder 'Utilities/Plotting'])

make_figure_MIZ_waves(MIZ_DATA,IS2_DATA)

%%
make_figure_MIZ_somewaves(MIZ_DATA,IS2_DATA)

%
make_figure_MIZ_nowaves(MIZ_DATA,IS2_DATA)