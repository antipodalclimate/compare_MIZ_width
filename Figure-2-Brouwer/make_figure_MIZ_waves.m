function make_figure_MIZ_waves(MIZ_DATA,IS2_DATA)

load_MIZ_waves; 

% Criteria for selection
usable_all = (Nvals > 100) & ~isnan(Dvals) &~isinf(SICvals) & Dvals < max(Dbins) & Dvals > min(Dbins); 
usable_all = usable_all; % & (timeval > 7 & timeval < 10); 
usable_all = usable_all & npoints > 1; 

usable = usable_all & wavytracks == 1; 

create_MIZ_wave_figure; 

%%
pos = [6.5 6]; 
set(gcf,'windowstyle','normal','position',[0 0 pos],'paperposition',[0 0 pos],'papersize',pos,'units','inches','paperunits','inches');
set(gcf,'windowstyle','normal','position',[0 0 pos],'paperposition',[0 0 pos],'papersize',pos,'units','inches','paperunits','inches');
print('/Users/chorvat/Brown Dropbox/Christopher Horvat/Apps/Overleaf/IS2-Waves-PM/Figures/MIZ-SIC-comp-waves','-dpdf','-r600');
%%