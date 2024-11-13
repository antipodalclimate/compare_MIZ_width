function make_figure_MIZ_nowaves(MIZ_DATA,IS2_DATA)

load_MIZ_waves; 


usable_all = (Nvals > 100) & ~isnan(Dvals) &~isinf(SICvals) & Dvals < max(Dbins) & Dvals > min(Dbins); 
usable_all = usable_all & (timeval > 7 & timeval < 10); 
usable_all = usable_all & npoints > 1; 

usable = usable_all & wavytracks < 1;

create_MIZ_wave_figure; 

used_tracks = IS2_DATA.namearray(unique(nameid(usable)));
writematrix(used_tracks,'Track_Lists/out_somewaves.txt')


pos = [6.5 6]; 
set(gcf,'windowstyle','normal','position',[0 0 pos],'paperposition',[0 0 pos],'papersize',pos,'units','inches','paperunits','inches');
set(gcf,'windowstyle','normal','position',[0 0 pos],'paperposition',[0 0 pos],'papersize',pos,'units','inches','paperunits','inches');
print([OS_string 'Apps/Overleaf/IS2-Waves-PM/Figures/MIZ-SIC-comp-somewaves'],'-dpdf','-r600');
%%