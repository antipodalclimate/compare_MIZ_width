% function make_figure_MIZ_waves(MIZ_DATA,IS2_DATA)

load_MIZ_waves; 

%%
usable = usable_all; %
usable = usable & wavytracks >= -1; 

used_tracks = IS2_DATA.namearray(unique(nameid(usable)));
intersections = unique(nameid(usable) + (beamid(usable)-1)*max(nameid));
fprintf('Using %2.2f million stencils from %2.0f beams across %2.0f tracks \n',sum(usable)/1e6,length(intersections),length(used_tracks))
writematrix(used_tracks,'Track_Lists/out_nowaves.txt')

%%

create_MIZ_nowave_figure; 


%%

create_MIZ_width_panel;


%%

allAxesInFigure = findall(gcf,'type','axes');
letter = {'(a)','(b)','(c)','(d)','(e)','(f)','(g)','(e)','(c)'};

for i = 1:length(allAxesInFigure)
    
 posy = get(allAxesInFigure(i),'position');

    set(allAxesInFigure(i),'fontname','times','fontsize',8,'xminortick','on','yminortick','on')
    
    annotation('textbox',[posy(1) - .025 posy(2)+posy(4) + .035 .025 .025], ...
        'String',letter{i},'LineStyle','none','FontName','Helvetica', ...
        'FontSize',8,'Tag','legtag');

end

pos = [6 4.5]; 
set(gcf,'windowstyle','normal','position',[0 0 pos],'paperposition',[0 0 pos],'papersize',pos,'units','inches','paperunits','inches');
set(gcf,'windowstyle','normal','position',[0 0 pos],'paperposition',[0 0 pos],'papersize',pos,'units','inches','paperunits','inches');
print([OS_string 'Apps/Overleaf/IS2-Waves-PM/Figures/MIZ-SIC-comp-nowaves'],'-dpdf','-r600');
%%