% function make_figure_MIZ_waves(MIZ_DATA,IS2_DATA)

load_MIZ_waves; 

%%
% Criteria for selection

usable_all = (Nsegvals > 1000) & usable_all; 
usable_all = usable_all & SICvals > 0.1 & LIFvals > 0.1;
usable_all = usable_all & npoints > 1; 

%%
usable = usable_all; %
usable = usable & wavytracks > 0; 

create_MIZ_wave_figure; 

used_tracks = IS2_DATA.namearray(unique(nameid(usable)));

%%
pos = [6.5 6]; 
set(gcf,'windowstyle','normal','position',[0 0 pos],'paperposition',[0 0 pos],'papersize',pos,'units','inches','paperunits','inches');
set(gcf,'windowstyle','normal','position',[0 0 pos],'paperposition',[0 0 pos],'papersize',pos,'units','inches','paperunits','inches');
print('/Users/chorvat/Brown Dropbox/Christopher Horvat/Apps/Overleaf/IS2-Waves-PM/Figures/MIZ-SIC-comp-waves','-dpdf','-r600');
%%