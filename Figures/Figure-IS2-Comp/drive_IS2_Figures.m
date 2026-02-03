% OS_string = '/Users/chorvat/Dropbox-Brown/Christopher Horvat/';
OS_string = '/Users/chorvat/Brown Dropbox/Christopher Horvat/';

OPTS.output_str = [OS_string OPTS.output_folder OPTS.savename]; 
% OPTS.output_str = [OS_string 'Research Projects/Active/Data/ICESat-2/PM-SIC-width/Along_Track_Statistics/AT_stats_KM_v6_all']; 

load(OPTS.output_str);

cutoff_N = 100; % need more than 100 segments for analysis to make sense
% cutoff_MIZ = 1; % Need more than 1 MIZ segment for analysis. 

close all

addpath([OPTS.code_folder 'Utilities/Plotting'])

%%

load_MIZ_waves; 

usable = usable_all; %
usable = usable & wavytracks >= -1; 

used_tracks = IS2_DATA.namearray(unique(nameid(usable)));
intersections = unique(nameid(usable) + (beamid(usable)-1)*max(nameid(usable)));
fprintf('Using %2.2f million stencils from %2.0f beams across %2.0f tracks \n',sum(usable)/1e6,length(intersections),length(used_tracks))
writematrix(used_tracks,'Track_Lists/out_all.txt')

create_composite_figure; 
create_MIZ_width_panel; 
%% 
%add_classified_panel; 

%%

allAxesInFigure = findall(gcf,'type','axes');
letter = {'(c)','(b)','(a)','(d)','(e)','(f)','(g)','(e)','(c)'};

for i = 1:length(allAxesInFigure)
    
 posy = get(allAxesInFigure(i),'position');

    set(allAxesInFigure(i),'fontname','times','fontsize',8,'xminortick','on','yminortick','on')
    
    annotation('textbox',[posy(1) - .05 posy(2)+posy(4)-.015 .025 .025], ...
        'String',letter{i},'LineStyle','none','FontName','Helvetica', ...
        'FontSize',8,'Tag','legtag');

end

pos = [6.5 4.5]; 
set(gcf,'windowstyle','normal','position',[0 0 pos],'paperposition',[0 0 pos],'papersize',pos,'units','inches','paperunits','inches');
set(gcf,'windowstyle','normal','position',[0 0 pos],'paperposition',[0 0 pos],'papersize',pos,'units','inches','paperunits','inches');
print([OS_string 'Apps/Overleaf/IS2-PM-SIC/Figures/MIZ-SIC-comp-all'],'-dpdf','-r600');

%%