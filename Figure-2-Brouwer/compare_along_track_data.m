% Updatable name of output. 


OPTS.output_str = '/Users/chorvat/Brown Dropbox/Christopher Horvat/Research Projects/Active/Data/ICESat-2/PM-SIC-width/Along_Track_Statistics/AT_stats_NH_v6'; 

%%
load(OPTS.output_str);
%%
make_figure_brouwer_comp(MIZ_DATA,IS2_DATA)