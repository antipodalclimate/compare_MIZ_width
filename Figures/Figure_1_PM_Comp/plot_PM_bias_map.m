% function fit_and_plot_PM_bias

close all
clear Ax

% This is the fit function
fitted = @(c) -1.639*(c-.6122).^2 + 0.2316;

disp('We are using the following fit')
disp(fitted)

%% Here's the fields we compare


sic_1 = AMSR_NT_SH(ICE_AMSR_NT & ICE_AMSR_BS);
sic_2 = AMSR_BS_SH(ICE_AMSR_NT & ICE_AMSR_BS);
cdr_std = CDR_std_daily_SH(ICE_AMSR_NT & ICE_AMSR_BS); 

% In this case we compare just NT2 and BS from the AMSR2 data. 

plotter = AMSR_NT_SH - AMSR_BS_SH;
plotter(ICE_AMSR_NT == 0) = nan;
plotter(ICE_AMSR_BS == 0) = nan;

plotter_MIZ = AMSR_NT_SH - AMSR_BS_SH;

% filter by some threshold

MIZfilt = AMSR_BS_SH >= 0.8;

plotter_MIZ(MIZfilt) = nan;
plotter_MIZ(ICE_AMSR_NT == 0) = nan;
plotter_MIZ(ICE_AMSR_BS == 0) = nan;

plotter_MIZ(repmat(sum(MIZfilt,3) < 6,[1 1 size(ICE_AMSR_BS,3)])) = nan;

%%

SICbins = linspace(0,1,51);
Bincent_SIC = 0.5*(SICbins(1:end-1) + SICbins(2:end));



[a,b,c] = histcounts(sic_2,SICbins);

upval = @(x) prctile(x,75);
dnval = @(x) prctile(x,25);

meanval = accumarray(c,sic_1,[length(SICbins)-1 1],@nanmedian);
plot_up = accumarray(c,sic_1,[length(SICbins)-1 1],upval);
plot_dn = accumarray(c,sic_1,[length(SICbins)-1 1],dnval);


meanval_del= accumarray(c,sic_1 - sic_2,[length(SICbins)-1 1],@nanmedian);
plot_up_del = accumarray(c,sic_1 - sic_2,[length(SICbins)-1 1],upval);
plot_dn_del = accumarray(c,sic_1 - sic_2,[length(SICbins)-1 1],dnval);

overmean_del = accumarray(c,(sic_1 - sic_2)./(1-sic_2),[length(SICbins)-1 1],@nanmedian);
overmean_del(isinf(overmean_del)) = 1;


%% Look at the uncertainty from the CDR


meanval_std= accumarray(c,cdr_std,[length(SICbins)-1 1],@nanmedian);
plot_up_std = accumarray(c,cdr_std,[length(SICbins)-1 1],upval);
plot_dn_std = accumarray(c,cdr_std,[length(SICbins)-1 1],dnval);

%% Here is the difference plot. 

%% 

climmer = [-.1 .3];


Ax{1} = subplot('position',[0 .5 .25 .45]);
worldmap([-90 -55],[-180 180]);
pcolorm(lat_SH,lon_SH,median(plotter,3,'omitnan'))
make_HR_coastlines([.6 .6 .6]);
set(gca,'clim',climmer)
% set(gca,'clim',[-.2 .2])
colorbar('position',[.55 .55 .025 .35]);
title('Bias: All Days','interpreter','latex')

Ax{2} = subplot('position',[.275 .5 .25 .45]);

worldmap([-90 -55],[-180 180]);
pcolorm(lat_SH,lon_SH,median(plotter_MIZ,3,'omitnan'))
make_HR_coastlines([.6 .6 .6]);
set(gca,'clim',climmer)
% set(gca,'clim',[-.2 .2])
% colorbar
title('Bias: Non-compact Days','interpreter','latex')

nc = 33;
cmapper = brewermap(nc,'PuOr');
cmapper(1:11,:) = [];
colormap(cmapper);

%%

Ax{3} = subplot('position',[.625 .5 .25 .45]);
worldmap([-90 -55],[-180 180]);

pcolorm(lat_SH,lon_SH,iqr(plotter_MIZ,3)./median(plotter_MIZ,3,'omitnan'))
make_HR_coastlines([.6 .6 .6]);

set(gca,'clim',[0 2])
colorbar('position',[.9 .55 .025 .35]);
title('IQR/Bias','interpreter','latex')
colormap(gca,brewermap(11,'paired'));
%
%%
Ax{4} =subplot('position',[.075 .1 .375 .35]);

jbfill(Bincent_SIC,plot_up',plot_dn',[.4 .4 .4],[1 1 1],1,.3);
hold on
plot(Bincent_SIC,meanval,'-k','linewidth',2)
xline(.8,'linewidth',1,'label','MIZ')
yline(.8,'linewidth',1,'label','MIZ')
line([0 1],[0 1],'Color',[.3 .3 .3],'linestyle','--')
title('SIC Comparison','interpreter','latex');

xlim([.15 1]);
ylim([.15 1]);

ylabel('AMSR2-NT2 SIC','interpreter','latex');
xlabel('AMSR2-BS SIC','interpreter','latex')
grid on; box on;

yyaxis right
set(gca,'ycolor','k','yticklabel','')
bar(Bincent_SIC,a/sum(a),1,'FaceColor',[.4 .4 .8],'FaceAlpha',.5,'EdgeColor','none');

xlim([.15 1]);
% ylim([.15 1]);

%%

Ax{5} = subplot('position',[.575 .1 .375 .35]);
cla

jbfill(Bincent_SIC,plot_up_del',plot_dn_del',[.4 .4 .4],[1 1 1],1,.3);
hold on
plot(Bincent_SIC,meanval_del,'-k','linewidth',2)
plot(Bincent_SIC(Bincent_SIC > 0.15),fitted(Bincent_SIC(Bincent_SIC > 0.15)),'--b','linewidth',1);


xline(.8,'linewidth',1,'label','MIZ');
yline(0,'--','linewidth',1)
xlim([.15 1]);
ylim([-.3 .4]);

% line([0 1],[0 1],'Color',[1 .4 .4],'linestyle','--')

plot(Bincent_SIC,1 - Bincent_SIC,'--r','linewidth',1)
plot(Bincent_SIC,-Bincent_SIC,'--r','linewidth',1)

ylabel('$\Delta$ SIC','interpreter','latex');

xlabel('AMSR2-BS SIC','interpreter','latex')
title('Algorithmic Bias','interpreter','latex');
grid on; box on;
%%

letter = {'(a)','(b)','(c)','(d)','(e)','(F)','(g)','(e)','(c)'};

delete(findall(gcf,'Tag','legtag'))

for i = 1:length(Ax)
    set(Ax{i},'fontname','helvetica','fontsize',9,'xminortick','on','yminortick','on')
    posy = get(Ax{i},'position');
    annotation('textbox',[posy(1) posy(2)+posy(4)-.025 .025 .025], ...
        'String',letter{i},'LineStyle','none','FontName','Helvetica', ...
        'FontSize',8,'Tag','legtag');

end

pos = [6.5 3.75];
set(gcf,'windowstyle','normal','position',[0 0 pos],'paperposition',[0 0 pos],'papersize',pos,'units','inches','paperunits','inches');
set(gcf,'windowstyle','normal','position',[0 0 pos],'paperposition',[0 0 pos],'papersize',pos,'units','inches','paperunits','inches');
print([OPTS.plot_save_str 'bias-BS-BS'],'-dpdf','-r600');
