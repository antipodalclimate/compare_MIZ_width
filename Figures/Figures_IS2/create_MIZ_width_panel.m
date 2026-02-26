% create_MIZ_width_panel

close


%%

Wbins = linspace(0,250,25);
ratbins = logspace(-2,1,20);

wCDR = MIZ_width_CDR(usable_beams);
wAMSR = MIZ_width_amsr(usable_beams);
wWAF = WAF_width(usable_beams); 

wCDR_all = MIZ_width_dist_CDR(beams_all)/1000;
wAMSR_all = MIZ_width_dist_amsr(beams_all)/1000;
wWAF_all = WAF_width_dist(beams_all)/1000; 

wCDR_wav = MIZ_width_dist_CDR(beams_waves)/1000;
wAMSR_wav = MIZ_width_dist_amsr(beams_waves)/1000;
wWAF_wav = WAF_width_dist(beams_waves)/1000; 

wCDR_nowav = MIZ_width_dist_CDR(beams_nowaves)/1000;
wAMSR_nowav = MIZ_width_dist_amsr(beams_nowaves)/1000;
wWAF_nowav = WAF_width_dist(beams_nowaves)/1000; 


% subplot('position',[.1 .1 .8 .15])

histogram(wCDR_wav,Wbins,'FaceColor',[.8 .2 .2],'FaceAlpha',0.5);
hold on

histogram(wAMSR_wav,Wbins,'FaceColor',[.2 .2 .8],'FaceAlpha',.5)
grid on; box on; 
% 
histogram(wWAF_wav,Wbins,'FaceColor',[.2 .8 .2],'FaceAlpha',.5)



legend('CDR','AMSR2','WAF')
hold off
xlim([0 max(Wbins)])
xlabel('km','interpreter','latex')
title('MIZ Width','interpreter','latex')

% %%
% subplot('position',[.55 .05 .35 .2])
% 
% 
% histogram(wCDR_wav,Wbins,'FaceColor',[.8 .2 .2],'FaceAlpha',0.5);
% hold on
% 
% histogram(wAMSR_wav,Wbins,'FaceColor',[.2 .2 .8],'FaceAlpha',.5)
% grid on; box on; 
% 
% histogram(wWAF_wav,Wbins,'FaceColor',[.2 .8 .2],'FaceAlpha',.5)


%%

allAxesInFigure = findall(gcf,'type','axes');
letter = {'(a)','(a)','(a)','(d)','(e)','(f)','(g)','(e)','(c)'};

for i = 1:length(allAxesInFigure)
    
 posy = get(allAxesInFigure(i),'position');

    set(allAxesInFigure(i),'fontname','times','fontsize',8,'xminortick','on','yminortick','on')
    
    annotation('textbox',[posy(1) - .05 posy(2)+posy(4)-.015 .025 .025], ...
        'String',letter{i},'LineStyle','none','FontName','Helvetica', ...
        'FontSize',8,'Tag','legtag');

end

pos = [6.5 2]; 
set(gcf,'windowstyle','normal','position',[0 0 pos],'paperposition',[0 0 pos],'papersize',pos,'units','inches','paperunits','inches');
set(gcf,'windowstyle','normal','position',[0 0 pos],'paperposition',[0 0 pos],'papersize',pos,'units','inches','paperunits','inches');
print([OPTS.plot_save_str 'MIZ-width-comp'],'-dpdf','-r600');
