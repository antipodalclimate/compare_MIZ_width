% This function compares PM data across the MIZ regions. 
% It is set up as a function. 

close all


xax = datetime(datestr(IS2_mutual));

plotcolors = [27,158,119
217,95,2
117,112,179]/256; 

%%
close all
xlimmer = [-3.5 3.5];
bins = -5:.2:5;

Ax{1} = subplot('position',[.075 .65 .65 .225]);
plot(xax,SIE_CDR,'k','linewidth',1);
hold on
plot(xax,SIE_AMSR_NT,'color',plotcolors(1,:),'linewidth',1)
% plot(xax,SIE_SSMI_BS,'--b','linewidth',1);
plot(xax,SIE_AMSR_BS,'color',plotcolors(2,:),'linewidth',1);
plot(xax,SIE_SSMI_BS,'color',plotcolors(3,:),'linewidth',1);
% plot(xax,SIE_ASI,'m','linewidth',1);

grid on; box on;
title('Antarctic Sea Ice Extent','interpreter','latex')
ylabel('Million km$^2$','interpreter','latex')
set(gca,'xticklabel','')

% yline(median(AMIZ_CDR-AMIZ_AMSR_NT),'--r','linewidth',1,'label',sprintf('%2.2f',median(AMIZ_CDR - AMIZ_AMSR_NT)),'fontsize',12)

legend('NSIDC-CDR','AMSR-NT2','AMSR2-BS','SSMI-BS','location',[.3 .95 .4 .025],'orientation','horizontal')

Ax{2} = subplot('position',[.8 .65 .125 .225]);

histogram(SIE_AMSR_NT-SIE_CDR,bins,'FaceColor',plotcolors(1,:),'Normalization','pdf')
hold on
% histogram(SIE_ASI - SIE_CDR,bins,'FaceColor','m','Normalization','pdf')
histogram(SIE_AMSR_BS - SIE_CDR,bins,'Facecolor',plotcolors(2,:),'Normalization','pdf')
histogram(SIE_SSMI_BS - SIE_CDR,bins,'Facecolor',plotcolors(3,:),'Normalization','pdf')
xlim(xlimmer)
grid on; box on;
view(90,-90)
title('$\Delta$ from CDR','interpreter','latex')
set(gca,'yticklabel','')

%

Ax{3} = subplot('position',[.075 .3625 .65 .225]);
plot(xax,SIA_CDR,'k','linewidth',1);
hold on
plot(xax,SIA_AMSR_NT,'color',plotcolors(1,:),'linewidth',1)
plot(xax,SIA_AMSR_BS,'color',plotcolors(2,:),'linewidth',1);
plot(xax,SIA_SSMI_BS,'color',plotcolors(3,:),'linewidth',1);
% plot(xax,SIA_ASI,'-m','linewidth',1);


hold off
grid on; box on;
title('Antarctic Sea Ice Area','interpreter','latex')
ylabel('Million km$^2$','interpreter','latex')
set(gca,'xticklabel','')

% legend('NSIDC-CDR','AMSR-NT2','AMSR2-BS','location',[.875 .8 .075 .1])

Ax{4} = subplot('position',[.8 .3625 .125 .225]);

histogram(SIA_AMSR_NT - SIA_CDR,bins,'FaceColor',plotcolors(1,:),'Normalization','pdf')
hold on
histogram(SIA_AMSR_BS - SIA_CDR,bins,'Facecolor',plotcolors(2,:),'Normalization','pdf')
histogram(SIA_SSMI_BS - SIA_CDR,bins,'Facecolor',plotcolors(3,:),'Normalization','pdf')
% histogram(SIA_ASI - SIA_CDR,bins,'Facecolor','m','Normalization','pdf')
xlim(xlimmer)
%
% xline(median(SIA_AMSR_NT - SIA_CDR),plotcolors(1,:),'linewidth',2)
% xline(median(SIA_AMSR_BS-SIA_CDR),'b','linewidth',2)
%
hold off
grid on; box on;
view(90,-90)
title('$\Delta$ from CDR','interpreter','latex')
set(gca,'yticklabel','')


Ax{5} = subplot('position',[.075 .075 .65 .225]);
plot(xax,AMIZ_CDR,'k','linewidth',1);
hold on
plot(xax,AMIZ_AMSR_NT,'color',plotcolors(1,:),'linewidth',1)
plot(xax,AMIZ_AMSR_BS,'color',plotcolors(2,:),'linewidth',1);
plot(xax,AMIZ_SSMI_BS,'color',plotcolors(3,:),'linewidth',1);
% plot(xax,AMIZ_ASI,'--m','linewidth',1);
hold off
grid on; box on;
title('Antarctic MIZ Extent','interpreter','latex')
ylabel('Million km$^2$','interpreter','latex')

% yline(median(AMIZ_AMSR_BS-AMIZ_CDR),'b','linewidth',2,'label',sprintf('%2.2f',median(AMIZ_CDR - AMIZ_AMSR_BS)),'fontsize',12)
%

% legend('NSIDC-CDR','AMSR-NT2','AMSR2-BS','location',[.875 .8 .075 .1])

Ax{6} = subplot('position',[.8 .075 .125 .225]);

histogram(AMIZ_AMSR_NT-AMIZ_CDR,bins,'FaceColor',plotcolors(1,:),'Normalization','pdf')
hold on
histogram(AMIZ_AMSR_BS-AMIZ_CDR,bins,'Facecolor',plotcolors(2,:),'Normalization','pdf')
histogram(AMIZ_SSMI_BS-AMIZ_CDR,bins,'Facecolor',plotcolors(3,:),'Normalization','pdf')
% histogram(AMIZ_ASI-AMIZ_CDR,bins,'Facecolor','m','Normalization','pdf')
xlim(xlimmer)
view(90,-90)
hold off
grid on; box on;
title('$\Delta$ from CDR','interpreter','latex')
set(gca,'yticklabel','')

letter = {'(a)','(c)','(b)','(d)','(c)','(e)','(g)','(e)','(c)'};

delete(findall(gcf,'Tag','legtag'))

for i = 1:length(Ax)
    set(Ax{i},'fontname','helvetica','fontsize',9,'xminortick','on','yminortick','on')
    posy = get(Ax{i},'position');
    annotation('textbox',[posy(1)-.04 posy(2)+posy(4)+.04 .025 .025], ...
        'String',letter{i},'LineStyle','none','FontName','Helvetica', ...
        'FontSize',8,'Tag','legtag');

end


pos = [6.5 3.5];
set(gcf,'windowstyle','normal','position',[0 0 pos],'paperposition',[0 0 pos],'papersize',pos,'units','inches','paperunits','inches');
set(gcf,'windowstyle','normal','position',[0 0 pos],'paperposition',[0 0 pos],'papersize',pos,'units','inches','paperunits','inches');
print([OPTS.plot_save_str 'SIE-SIA-comp'],'-dpdf','-r600');

%%


