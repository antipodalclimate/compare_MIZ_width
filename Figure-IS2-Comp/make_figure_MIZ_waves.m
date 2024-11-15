% function make_figure_MIZ_waves(MIZ_DATA,IS2_DATA)

load_MIZ_waves; 

%%
% Criteria for selection

usable_all = (Nsegvals > cutoff_N) & usable_all; 
usable_all = usable_all & SICvals > 0.15 & LIFvals > 0.15;
usable_all = usable_all & npoints > 1; 

%%
usable = usable_all; %
usable = usable & wavytracks > 0; 

used_tracks = IS2_DATA.namearray(unique(nameid(usable)));
intersections = unique(nameid(usable) + (beamid(usable)-1)*max(nameid(usable)));
fprintf('Using %2.2f million stencils from %2.0f beams across %2.0f tracks \n',sum(usable)/1e6,length(intersections),length(used_tracks))
writematrix(used_tracks,'Track_Lists/out_waves.txt')

%%
create_MIZ_wave_figure; 



%%
pos = [6.5 6]; 
set(gcf,'windowstyle','normal','position',[0 0 pos],'paperposition',[0 0 pos],'papersize',pos,'units','inches','paperunits','inches');
set(gcf,'windowstyle','normal','position',[0 0 pos],'paperposition',[0 0 pos],'papersize',pos,'units','inches','paperunits','inches');
print([OS_string 'Apps/Overleaf/IS2-Waves-PM/Figures/MIZ-SIC-comp-waves'],'-dpdf','-r600');
%%