load_MIZ_waves; 

usable_all = (Nsegvals > 100) & usable_all; 
usable_all = usable_all & SICvals > 0.1 & LIFvals > 0.1;
usable_all = usable_all; %  & (timeval > 7 & timeval < 10); 
usable_all = usable_all; %  & (isstrong == 1); 
usable_all = usable_all; % & npoints > 1; 

usable = usable_all;

create_MIZ_wave_figure; 


used_tracks = IS2_DATA.namearray(unique(nameid(usable)));
writematrix(used_tracks,'Track_Lists/out_all.txt')


pos = [6.5 6]; 
set(gcf,'windowstyle','normal','position',[0 0 pos],'paperposition',[0 0 pos],'papersize',pos,'units','inches','paperunits','inches');
set(gcf,'windowstyle','normal','position',[0 0 pos],'paperposition',[0 0 pos],'papersize',pos,'units','inches','paperunits','inches');
print([OS_string 'Apps/Overleaf/IS2-Waves-PM/Figures/MIZ-SIC-comp-all'],'-dpdf','-r600');

%%