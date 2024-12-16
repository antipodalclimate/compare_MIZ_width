% OS_string = '/Users/chorvat/Dropbox-Brown/Christopher Horvat/';
OS_string = '/Users/chorvat/Brown Dropbox/Christopher Horvat/';

OPTS.output_str = [OS_string 'Research Projects/Active/Data/ICESat-2/PM-SIC-width/Along_Track_Statistics/AT_stats_SH_v6_all']; 
% OPTS.output_str = [OS_string 'Research Projects/Active/Data/ICESat-2/PM-SIC-width/Along_Track_Statistics/AT_stats_KM_v6_all']; 

load(OPTS.output_str);

cutoff_N = 100; % need more than 100 segments for analysis to make sense
% cutoff_MIZ = 1; % Need more than 1 MIZ segment for analysis. 

close all

addpath([OPTS.code_folder 'Utilities/Plotting'])

%%

load_MIZ_waves; 

%%

usable = usable_all; %
usable = usable & wavytracks >= -1; 

used_tracks = IS2_DATA.namearray(unique(nameid(usable)));
intersections = unique(nameid(usable) + (beamid(usable)-1)*max(nameid(usable)));
fprintf('Using %2.2f million stencils from %2.0f beams across %2.0f tracks \n',sum(usable)/1e6,length(intersections),length(used_tracks))
writematrix(used_tracks,'Track_Lists/out_all.txt')

%%

create_composite_figure; 


pos = [6.5 6]; 
set(gcf,'windowstyle','normal','position',[0 0 pos],'paperposition',[0 0 pos],'papersize',pos,'units','inches','paperunits','inches');
set(gcf,'windowstyle','normal','position',[0 0 pos],'paperposition',[0 0 pos],'papersize',pos,'units','inches','paperunits','inches');
print([OS_string 'Apps/Overleaf/IS2-Waves-PM/Figures/MIZ-SIC-comp-all'],'-dpdf','-r600');

%%