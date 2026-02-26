% OS_string = '/Users/chorvat/Dropbox-Brown/Christopher Horvat/';
OS_string = '/Users/chorvat/Library/CloudStorage/Dropbox-Brown/Christopher Horvat/';

OPTS.output_str = [OS_string 'Research Projects/Active/Data/ICESat-2/PM-SIC-width/Along_Track_Statistics/AT_stats_SH_v6_all']; 
% OPTS.output_str = [OS_string 'Research Projects/Active/Data/ICESat-2/PM-SIC-width/Along_Track_Statistics/AT_stats_KM_v6_all']; 

load(OPTS.output_str);

cutoff_N = 100; % need more than 100 segments for analysis to make sense
% cutoff_MIZ = 1; % Need more than 1 MIZ segment for analysis. 

close all

addpath([OPTS.code_folder 'Utilities/Plotting'])
addpath([OPTS.code_folder 'Figure-IS2-Comp'])

%%

load_MIZ_waves; 
usable = usable_all; %
usable = usable; % & Hvals < 0.1 & Hvals > -.1 & WAFvals < 0.075; 

%%
used_tracks = IS2_DATA.namearray(unique(nameid(usable)));
intersections = unique(nameid(usable) + (beamid(usable)-1)*max(nameid(usable)));
fprintf('Using %2.2f million stencils from %2.0f beams across %2.0f tracks \n',sum(usable)/1e6,length(intersections),length(used_tracks))
make_ancillary_figure_data; 

%%

figure; 


%%

create_composite_figure_pres; 

%%
pos = [4 4]; 
set(gcf,'windowstyle','normal','position',[0 0 pos],'paperposition',[0 0 pos],'papersize',pos,'units','inches','paperunits','inches');
set(gcf,'windowstyle','normal','position',[0 0 pos],'paperposition',[0 0 pos],'papersize',pos,'units','inches','paperunits','inches');

create_bias_figure_pres; 

%% 

pos = [4 4]; 
set(gcf,'windowstyle','normal','position',[0 0 pos],'paperposition',[0 0 pos],'papersize',pos,'units','inches','paperunits','inches');
set(gcf,'windowstyle','normal','position',[0 0 pos],'paperposition',[0 0 pos],'papersize',pos,'units','inches','paperunits','inches');


highbias_fig_pres; 