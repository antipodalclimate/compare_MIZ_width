% Updatable name of output. 


OPTS.output_str = '/Users/chorvat/Brown Dropbox/Christopher Horvat/Research Projects/Active/Data/ICESat-2/PM-SIC-width/Along_Track_Statistics/AT_stats_SH_v6'; 

%%
load(OPTS.output_str);
%%
make_figure_MIZ_waves(MIZ_DATA,IS2_DATA)

%%
make_figure_MIZ_nowaves(MIZ_DATA,IS2_DATA)