% Must run the Drive_Figures before this. 
plotter = AMSR_NT_SH - AMSR_BS_SH;
plotter(ICE_AMSR_NT == 0) = nan;
plotter(ICE_AMSR_BS == 0) = nan;

plotter_MIZ = AMSR_NT_SH - AMSR_BS_SH;
plotter_MIZ(MIZ_AMSR_BS ~= 1) = nan;
plotter_MIZ(ICE_AMSR_NT == 0) = nan;
plotter_MIZ(ICE_AMSR_BS == 0) = nan;

plotter_MIZ(repmat(sum(MIZ_AMSR_NT,3) < 6,[1 1 size(ICE_AMSR_BS,3)])) = nan;

%% One day bias
climmer = [-.1 .3];

plot_index = 1844; 

worldmap([-90 -55],[-180 180]);
pcolorm(lat_SH,lon_SH,plotter(:,:,plot_index))
make_HR_coastlines([.6 .6 .6]);
set(gca,'clim',climmer)
% set(gca,'clim',[-.2 .2])
% colorbar
title(['NT2-BS: ' datestr(datenum_coincedent(plot_index))],'interpreter','latex','fontsize',24)

nc = 33;
cmapper = brewermap(nc,'PuOr');
cmapper(1:11,:) = [];
colormap(cmapper);

pos = [8 8];
set(gcf,'windowstyle','normal','position',[0 0 pos],'paperposition',[0 0 pos],'papersize',pos,'units','inches','paperunits','inches');
set(gcf,'windowstyle','normal','position',[0 0 pos],'paperposition',[0 0 pos],'papersize',pos,'units','inches','paperunits','inches');
print(['oneday-bias'],'-dpdf','-r600');

ignore = (AMSR_NT_SH(:,:,plot_index) == 0) | ( CDR_SIC_SH(:,:,plot_index) == 0);

darea = plotter(:,:,plot_index).*area_SH;
darea(ignore) = nan; 

area1 = AMSR_NT_SH(:,:,plot_index).*area_SH;
area1(ignore) = nan; 


area2 = CDR_SIC_SH(:,:,plot_index).*area_SH;
area2(ignore) = nan; 

darea = nansum(darea,[1 2])
area1 = nansum(area1,[1 2])
area2 = nansum(area2,[1 2])
100*darea/area2


%%



worldmap([-90 -55],[-180 180]);
pcolorm(lat_SH,lon_SH,median(plotter,3,'omitnan'))
make_HR_coastlines([.6 .6 .6]);
set(gca,'clim',climmer)
% set(gca,'clim',[-.2 .2])
% colorbar
title('Bias: NT2-BS - All','interpreter','latex')

nc = 33;
cmapper = brewermap(nc,'PuOr');
cmapper(1:11,:) = [];
colormap(cmapper);

pos = [8 8];
set(gcf,'windowstyle','normal','position',[0 0 pos],'paperposition',[0 0 pos],'papersize',pos,'units','inches','paperunits','inches');
set(gcf,'windowstyle','normal','position',[0 0 pos],'paperposition',[0 0 pos],'papersize',pos,'units','inches','paperunits','inches');
print(['global-bias'],'-dpdf','-r600');


%%

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

subplot('position',[.1 .1 .65 .7]);

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




subplot('position',[.8 .1 .125 .7]);

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


pos = [6.5 3.5];
set(gcf,'windowstyle','normal','position',[0 0 pos],'paperposition',[0 0 pos],'papersize',pos,'units','inches','paperunits','inches');
set(gcf,'windowstyle','normal','position',[0 0 pos],'paperposition',[0 0 pos],'papersize',pos,'units','inches','paperunits','inches');
print(['SIE'],'-dpdf','-r600');

%%
close all
xlimmer = [-3.5 3.5];
bins = -5:.2:5;

subplot('position',[.1 .1 .65 .7]);

plot(xax,SIA_CDR,'k','linewidth',1);
hold on
plot(xax,SIA_AMSR_NT,'color',plotcolors(1,:),'linewidth',1)
% plot(xax,SIA_SSMI_BS,'--b','linewidth',1);
plot(xax,SIA_AMSR_BS,'color',plotcolors(2,:),'linewidth',1);
plot(xax,SIA_SSMI_BS,'color',plotcolors(3,:),'linewidth',1);
% plot(xax,SIA_ASI,'m','linewidth',1);

grid on; box on;


title('Antarctic Sea Ice Area','interpreter','latex')
ylabel('Million km$^2$','interpreter','latex')
set(gca,'xticklabel','')

% yline(median(AMIZ_CDR-AMIZ_AMSR_NT),'--r','linewidth',1,'label',sprintf('%2.2f',median(AMIZ_CDR - AMIZ_AMSR_NT)),'fontsize',12)

legend('NSIDC-CDR','AMSR-NT2','AMSR2-BS','SSMI-BS','location',[.3 .95 .4 .025],'orientation','horizontal')



subplot('position',[.8 .1 .125 .7]);

histogram(SIA_AMSR_NT-SIA_CDR,bins,'FaceColor',plotcolors(1,:),'Normalization','pdf')
hold on
% histogram(SIA_ASI - SIA_CDR,bins,'FaceColor','m','Normalization','pdf')
histogram(SIA_AMSR_BS - SIA_CDR,bins,'Facecolor',plotcolors(2,:),'Normalization','pdf')
histogram(SIA_SSMI_BS - SIA_CDR,bins,'Facecolor',plotcolors(3,:),'Normalization','pdf')
xlim(xlimmer)
grid on; box on;
view(90,-90)
title('$\Delta$ from CDR','interpreter','latex')
set(gca,'yticklabel','')


pos = [6.5 3.5];
set(gcf,'windowstyle','normal','position',[0 0 pos],'paperposition',[0 0 pos],'papersize',pos,'units','inches','paperunits','inches');
set(gcf,'windowstyle','normal','position',[0 0 pos],'paperposition',[0 0 pos],'papersize',pos,'units','inches','paperunits','inches');
print(['SIA'],'-dpdf','-r600');


%%


close all
xlimmer = [-3.5 3.5];
bins = -5:.2:5;

subplot('position',[.1 .1 .65 .7]);

plot(xax,AMIZ_CDR,'k','linewidth',1);
hold on
plot(xax,AMIZ_AMSR_NT,'color',plotcolors(1,:),'linewidth',1)
% plot(xax,AMIZ_SSMI_BS,'--b','linewidth',1);
plot(xax,AMIZ_AMSR_BS,'color',plotcolors(2,:),'linewidth',1);
plot(xax,AMIZ_SSMI_BS,'color',plotcolors(3,:),'linewidth',1);
% plot(xax,AMIZ_ASI,'m','linewidth',1);

grid on; box on;


title('Antarctic MIZ Extent','interpreter','latex')
ylabel('Million km$^2$','interpreter','latex')
set(gca,'xticklabel','')

% yline(median(AMIZ_CDR-AMIZ_AMSR_NT),'--r','linewidth',1,'label',sprintf('%2.2f',median(AMIZ_CDR - AMIZ_AMSR_NT)),'fontsize',12)

legend('NSIDC-CDR','AMSR-NT2','AMSR2-BS','SSMI-BS','location',[.3 .95 .4 .025],'orientation','horizontal')




subplot('position',[.8 .1 .125 .7]);

histogram(AMIZ_AMSR_NT-AMIZ_CDR,bins,'FaceColor',plotcolors(1,:),'Normalization','pdf')
hold on
% histogram(AMIZ_ASI - AMIZ_CDR,bins,'FaceColor','m','Normalization','pdf')
histogram(AMIZ_AMSR_BS - AMIZ_CDR,bins,'Facecolor',plotcolors(2,:),'Normalization','pdf')
histogram(AMIZ_SSMI_BS - AMIZ_CDR,bins,'Facecolor',plotcolors(3,:),'Normalization','pdf')
xlim(xlimmer)
grid on; box on;
view(90,-90)
title('$\Delta$ from CDR','interpreter','latex')
set(gca,'yticklabel','')


pos = [6.5 3.5];
set(gcf,'windowstyle','normal','position',[0 0 pos],'paperposition',[0 0 pos],'papersize',pos,'units','inches','paperunits','inches');
set(gcf,'windowstyle','normal','position',[0 0 pos],'paperposition',[0 0 pos],'papersize',pos,'units','inches','paperunits','inches');
print(['AMIZ'],'-dpdf','-r600');


%%

close all

jbfill(Bincent,plot_up_del',plot_dn_del',[.2 .2 .8],[1 1 1],1,.5);
hold on
plot(Bincent,meanval_del,'-k','linewidth',2)


%%

xline(.8,'linewidth',1,'label','MIZ');
yline(0,'--','linewidth',1)
xlim([.15 1]);
ylim([-.3 .4]);

% line([0 1],[0 1],'Color',[1 .4 .4],'linestyle','--')

plot(Bincent,1 - Bincent,'--r','linewidth',1)
plot(Bincent,-Bincent,'--r','linewidth',1)

ylabel('Bias','interpreter','latex');

xlabel('CDR SIC','interpreter','latex')
title('Inversion Uncertainty','interpreter','latex');
grid on; box on;


hold on

jbfill(Bincent,abs(meanval_std)',-abs(meanval_std)',[.8 .4 .4],[1 1 1],1,.5);
set(gca,'fontsize',18)


pos = [6.5 4.5];
set(gcf,'windowstyle','normal','position',[0 0 pos],'paperposition',[0 0 pos],'papersize',pos,'units','inches','paperunits','inches');
set(gcf,'windowstyle','normal','position',[0 0 pos],'paperposition',[0 0 pos],'papersize',pos,'units','inches','paperunits','inches');

print(['Alg-bias-0'],'-dpdf','-r600');
