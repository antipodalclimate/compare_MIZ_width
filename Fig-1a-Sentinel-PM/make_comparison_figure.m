
% Which IS-2 beam to choose. 
beam_pickout = 6; 

close all
horvat_colors; 

Ax{1} = subplot('position',[.05 .455 .45 .5]);

% Make the map. 
latlim = sort([0.5 .25] + lat_span); 
lonlim = sort([-1.75 2.25] + lon_span);
worldmap((latlim),(lonlim)); 

setm(gca,'PlabelLocation',-70,'PlineLocation',-70,'PLabelRound',0);
setm(gca,'MLabelLocation',[-84 -82],'MLineLocation',[-84 -82],'MLabelRound',0,'MLabelParallel','south');
setm(gca,'grid','on','GLineWidth',0.5,'GLineStyle','--','GColor','w','fontsize',8);

% Reduce size of amplitude plot. 
skipper = 4; 

inds_plot = lat_S1(:) < latlim(2) & lat_S1(:) > latlim(1) & lon_S1(:) < lonlim(2) & lon_S1(:) > lonlim(1);
Amp_plot = nan*Amp_S1; 
Amp_plot(inds_plot) = Amp_S1(inds_plot); 
Amp_plot(Amp_plot > 1000) = 1000; 
Amp_plot = Amp_plot / 1000; 

[a,b] = ind2sub(size(Amp_S1),find(inds_plot));
l1 = min(a):skipper:max(a); 
l2 = min(b):skipper:max(b); 

pcolorm(lat_S1(l1,l2),lon_S1(l1,l2),Amp_plot(l1,l2));

set(gca,'clim',[0.5 1])
colormap((brewermap(9,'-greys')))
plotm(IS2_obj{beam_pickout}.lat(IS2_obj{beam_pickout}.is_ice),IS2_obj{beam_pickout}.lon(IS2_obj{beam_pickout}.is_ice),'--r','linewidth',1)

% s = scatterm(IS2_obj{beam_pickout}.lat(IS2_obj{beam_pickout}.is_ice),IS2_obj{beam_pickout}.lon(IS2_obj{beam_pickout}.is_ice),20,'k','s','linewidth',3);

scatterm(lat_CIZ_WAF(beam_pickout),lon_CIZ_WAF(beam_pickout),100,clabs(1,:),'s','linewidth',1.5);
scatterm(lat_CIZ_LIF(beam_pickout),lon_CIZ_LIF(beam_pickout),100,clabs(2,:),'s','linewidth',1.5);
scatterm(lat_CIZ_S1(beam_pickout),lon_CIZ_S1(beam_pickout),100,clabs(3,:),'s','linewidth',1.5);
scatterm(lat_CIZ_ISPM(beam_pickout),lon_CIZ_ISPM(beam_pickout),100,clabs(5,:),'s','linewidth',1.5);
scatterm(lat_CIZ_PM(beam_pickout),lon_CIZ_PM(beam_pickout),100,clabs(7,:),'s','linewidth',1.5);

%%

% add_coastlines; 
%
AT_dist = IS2_obj{beam_pickout}.dist/1000; 
AT_isice = IS2_obj{beam_pickout}.is_ice; 

Ax{2} = subplot('position',[.1 .15 .8 .2],'replace');
plot(AT_dist(AT_isice),100*AT_stats{beam_pickout}.height_adj(AT_isice),'--','linewidth',0.05,'color',[.7 .7 .7])
hold on
plot(AT_dist(AT_isice),100*height_smooth{beam_pickout},'k','linewidth',1)
plot(AT_dist(AT_isice),100*height_smooth{beam_pickout} + 100*height_std{beam_pickout},'--k','linewidth',1)
plot(AT_dist(AT_isice),100*height_smooth{beam_pickout} - 100*height_std{beam_pickout},'--k','linewidth',1)

ylim([-.2 .5]*100)
xlim([0 100]); 

xline(X_CIZ_WAF(beam_pickout),'color',clabs(1,:),'linewidth',3)
xline(X_CIZ_LIF(beam_pickout),'color',clabs(2,:),'linewidth',3)
xline(X_CIZ_S1(beam_pickout),'color',clabs(3,:),'linewidth',3)
xline(X_CIZ_ISPM(beam_pickout),'color',clabs(5,:),'linewidth',3)
xline(X_CIZ_PM(beam_pickout),'color',clabs(7,:),'linewidth',3)

grid on; box on; 
ylabel('cm','interpreter','latex');
xlabel('Distance from 1st ice segment','interpreter','latex');

%
Ax{3} = subplot('position',[.525 .5 .35 .45],'replace');


plot(AT_dist,AT_stats{beam_pickout}.WAF,'linewidth',1,'color',clabs(1,:))
hold on
plot(AT_dist,AT_stats{beam_pickout}.LIF,'linewidth',1,'color',clabs(2,:))
plot(AT_dist,AT_S1_LIF{beam_pickout},'linewidth',1,'color',clabs(3,:))
plot(AT_dist,AT_stats{beam_pickout}.SIC,'linewidth',1,'color',clabs(5,:))
plot(AT_dist,AT_stats{beam_pickout}.SIC_amsr,'--','linewidth',1,'color',clabs(5,:))
plot(AT_dist,AT_sic_PM{beam_pickout},'linewidth',1,'color',clabs(7,:))


xlim([0 150]); 
ylim([0 1])
grid on; box on; 
xlabel('Distance from 1st ice segment','interpreter','latex');

% yyaxis right
% set(gca,'ycolor','k')
% plot(xvals,100*height_smooth,'k','linewidth',1)
xlim([0 100]); 
grid on; box on; 
ylabel('cm','interpreter','latex');
% plot(xvals,100*height_smooth + 100*height_std,'--k','linewidth',1)
% plot(xvals,100*height_smooth - 100*height_std,'--k','linewidth',1)
xline(X_CIZ_WAF(beam_pickout),'color',clabs(1,:),'linewidth',3)
xline(X_CIZ_LIF(beam_pickout),'color',clabs(2,:),'linewidth',3)
xline(X_CIZ_S1(beam_pickout),'color',clabs(3,:),'linewidth',3)
xline(X_CIZ_ISPM(beam_pickout),'color',clabs(5,:),'linewidth',3)
xline(X_CIZ_PM(beam_pickout),'color',clabs(7,:),'linewidth',3)

h = legend('AT-WAF','AT-LIF','SAR-SIC','IS2-AMSR','IS2-SSMI','CDR','location','southeast');
set(h,'ItemTokenSize',[20 20]);
%

letter = {'(a)','(c)','(b)','(d)','(e)','(f)','(g)','(e)','(c)'};

delete(findall(gcf,'Tag','legtag'))

for i = 1:length(Ax)
    set(Ax{i},'fontname','helvetica','fontsize',9,'xminortick','on','yminortick','on')
    posy = get(Ax{i},'position');
    annotation('textbox',[posy(1)-.04 posy(2)+posy(4)+.04 .025 .025], ...
        'String',letter{i},'LineStyle','none','FontName','Helvetica', ...
        'FontSize',8,'Tag','legtag');
    
end


pos = [6.5 3.25]; 
set(gcf,'windowstyle','normal','position',[0 0 pos],'paperposition',[0 0 pos],'papersize',pos,'units','inches','paperunits','inches');
set(gcf,'windowstyle','normal','position',[0 0 pos],'paperposition',[0 0 pos],'papersize',pos,'units','inches','paperunits','inches');


print('/Users/chorvat/Dropbox (Brown)/Apps/Overleaf/IS2-Waves-PM/Figures/SAR-Comp','-dpdf','-r600');
% print('/Users/chorvat/Dropbox (Brown)/Christopher Horvat/Apps/Overleaf/2024-NASA-ROSES-PM/Proposal/Figures/SAR-Comp','-dpdf','-r600');

