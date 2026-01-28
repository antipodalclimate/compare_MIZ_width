
clear
close all

AMSR_loc = '/Users/chorvat/Dropbox-Brown/Christopher Horvat/Research Projects/Active/Data/SIC-Data/AMSR2-NT/AMSR2_SIC_daily.mat'; 
SSMI_loc = '/Users/chorvat/Dropbox-Brown/Christopher Horvat/Research Projects/Active/Data/SIC-Data/NSIDC-CDR/Daily/NSIDC-CDR_daily.mat'; 
% ASI_loc = '/Users/chorvat/Brown Dropbox/Christopher Horvat/Research Projects/Active/Data/SIC-Data/AMSR2-ASI/AMSR2_ASI_daily.mat'; 

load(AMSR_loc,'AMSR_datenum','AMSR_NT_SH','AMSR_BS_SH');
load(SSMI_loc,'CDR_SIC_SH','lat_SH','lon_SH','CDR_time_SH','BS_SIC_SH','NT_SIC_SH','area_SH');
% load(ASI_loc,'AMSR_ASI_SH','ASI_time_SH','lat_ASI_SH','lon_ASI_SH','area_ASI_SH');

load('dist_to_coast'); 

%% Ignore close-to-coast 

coastmask = 1*(dist_to_coast > 25); 
coastmask(coastmask == 0) = nan; 

% coastmask_ASI = 1*(dist_to_coast_ASI > 25); 
% coastmask_ASI(coastmask_ASI == 0) = nan; 
% 
%%

% Same dates in the IS2 period. 
IS2_datenum = datenum('01-Oct-2018'):datenum('31-Dec-2024'); 

[mutual_datenum,ia,ib] = intersect(AMSR_datenum,CDR_time_SH);
[IS2_mutual,ic,~] = intersect(mutual_datenum,IS2_datenum);
% [ASI_mutual,id,ie] = intersect(IS2_mutual,ASI_time_SH); 

% NASATEAM2 AMSR value
AMSR_NT_SH = bsxfun(@times,AMSR_NT_SH(:,:,ia(ic)),coastmask);
% CDR values
CDR_SIC_SH = bsxfun(@times,CDR_SIC_SH(:,:,ib(ic)),coastmask);
SSMI_BS_SH = bsxfun(@times,BS_SIC_SH(:,:,ib(ic)),coastmask);
SSMI_NT_SH = bsxfun(@times,NT_SIC_SH(:,:,ib(ic)),coastmask);

% AMSR_ASI_SH = bsxfun(@times,AMSR_ASI_SH(:,:,ie),coastmask_ASI');

clear BS_SIC_SH NT_SIC_SH

% BOOTSTRAP AMSR values
AMSR_BS_SH =  bsxfun(@times,AMSR_BS_SH(:,:,ia(ic)),coastmask);

AMSR_NT_SH(AMSR_NT_SH > 1) = nan; 
AMSR_BS_SH(AMSR_BS_SH > 1) = nan; 

%% Now build more fields from the gridded data

MIZ_AMSR_NT = AMSR_NT_SH > 0.15 & AMSR_NT_SH < 0.8; 
MIZ_AMSR_BS = AMSR_BS_SH > 0.15 & AMSR_BS_SH < 0.8; 
MIZ_CDR = CDR_SIC_SH > 0.15 & CDR_SIC_SH < 0.8; 
MIZ_SSMI_NT = SSMI_NT_SH > 0.15 & SSMI_NT_SH < 0.8; 
MIZ_SSMI_BS = SSMI_BS_SH > 0.15 & SSMI_BS_SH < 0.8; 
% MIZ_ASI = AMSR_ASI_SH > 0.15 & AMSR_ASI_SH < 0.8; 

ICE_AMSR_NT = AMSR_NT_SH > 0.15; 
ICE_AMSR_BS = AMSR_BS_SH > 0.15; 
ICE_CDR = CDR_SIC_SH > 0.15; 
ICE_SSMI_NT = SSMI_NT_SH > 0.15; 
ICE_SSMI_BS = SSMI_BS_SH > 0.15; 
% ICE_ASI = AMSR_ASI_SH > 0.15; 

prefac = 1/1e6; 

AMIZ_AMSR_BS = prefac*squeeze(sum(bsxfun(@times,area_SH,MIZ_AMSR_BS),[1 2],'omitnan'));
AMIZ_AMSR_NT = prefac*squeeze(sum(bsxfun(@times,area_SH,MIZ_AMSR_NT),[1 2],'omitnan'));
AMIZ_CDR = prefac*squeeze(sum(bsxfun(@times,area_SH,MIZ_CDR),[1 2],'omitnan')); 
AMIZ_SSMI_NT = prefac*squeeze(sum(bsxfun(@times,area_SH,MIZ_SSMI_NT),[1 2],'omitnan'));
AMIZ_SSMI_BS = prefac*squeeze(sum(bsxfun(@times,area_SH,MIZ_SSMI_BS),[1 2],'omitnan'));
% AMIZ_ASI = prefac*squeeze(sum(bsxfun(@times,area_ASI_SH,MIZ_ASI),[1 2],'omitnan')); 


SIA_AMSR_NT = prefac*squeeze(sum(bsxfun(@times,area_SH,AMSR_NT_SH),[1 2],'omitnan'));
SIA_AMSR_BS = prefac*squeeze(sum(bsxfun(@times,area_SH,AMSR_BS_SH),[1 2],'omitnan'));
SIA_CDR = prefac*squeeze(sum(bsxfun(@times,area_SH,CDR_SIC_SH),[1 2],'omitnan'));
SIA_SSMI_NT = prefac*squeeze(sum(bsxfun(@times,area_SH,SSMI_NT_SH),[1 2],'omitnan'));
SIA_SSMI_BS = prefac*squeeze(sum(bsxfun(@times,area_SH,SSMI_BS_SH),[1 2],'omitnan'));
% SIA_ASI = prefac*squeeze(sum(bsxfun(@times,area_ASI_SH,AMSR_ASI_SH),[1 2],'omitnan'));

SIE_AMSR_NT = prefac*squeeze(sum(bsxfun(@times,area_SH,ICE_AMSR_NT),[1 2],'omitnan'));
SIE_AMSR_BS = prefac*squeeze(sum(bsxfun(@times,area_SH,ICE_AMSR_BS),[1 2],'omitnan'));
SIE_CDR = prefac*squeeze(sum(bsxfun(@times,area_SH,ICE_CDR),[1 2],'omitnan')); 
SIE_SSMI_NT = prefac*squeeze(sum(bsxfun(@times,area_SH,ICE_SSMI_NT),[1 2],'omitnan'));
SIE_SSMI_BS = prefac*squeeze(sum(bsxfun(@times,area_SH,ICE_SSMI_BS),[1 2],'omitnan'));
% SIE_ASI = prefac*squeeze(sum(bsxfun(@times,area_ASI_SH,ICE_ASI),[1 2],'omitnan')); 


%%
% Six products
% AMSR2 - NT
% AMSR2 - BS
% AMSR2 - ASI
% SSMI - CDR
% SSMI - BS
% SSMI - NT

% Interesected in official product differences
% AMSR2 - SSMI-CDR
dAMIZ = AMIZ_AMSR_NT - AMIZ_CDR; 

% Then interested in same algorithms, different sensors
dAMIZ_BS = AMIZ_AMSR_BS - AMIZ_SSMI_BS; 
dAMIZ_NT = AMIZ_AMSR_NT - AMIZ_SSMI_NT; 

% Then interested in same sensor, different algorithms
dAMIZ_AMSR = AMIZ_AMSR_NT - AMIZ_AMSR_BS; 
dAMIZ_SSMI = AMIZ_SSMI_NT - AMIZ_SSMI_BS; 

% Now for shits, different sensor, different algos
dAMIZ_cross = dAMIZ_NT + dAMIZ_SSMI; 

dSIE = SIE_AMSR_NT - SIE_CDR; 

% Now SIA differences
dSIA = SIA_AMSR_NT - SIA_CDR; 
dSIA_NT = SIA_AMSR_NT - SIA_SSMI_NT;

xax = datetime(datestr(IS2_mutual)); 

%%
close all
xlimmer = [-5 5];
bins = -5:.2:5; 

Ax{1} = subplot('position',[.075 .65 .65 .225]);
plot(xax,SIE_CDR,'k','linewidth',1); 
hold on
plot(xax,SIE_AMSR_NT,'r','linewidth',1)
% plot(xax,SIE_SSMI_BS,'--b','linewidth',1); 
plot(xax,SIE_AMSR_BS,'b','linewidth',1); 
% plot(xax,SIE_ASI,'m','linewidth',1); 

grid on; box on; 
title('Sea Ice Extent','interpreter','latex')
ylabel('Million km$^2$','interpreter','latex')
set(gca,'xticklabel','')

% yline(median(AMIZ_CDR-AMIZ_AMSR_NT),'--r','linewidth',1,'label',sprintf('%2.2f',median(AMIZ_CDR - AMIZ_AMSR_NT)),'fontsize',12)

legend('CDR','AMSR2','AMSR2-BS','location',[.3 .95 .4 .025],'orientation','horizontal')

Ax{2} = subplot('position',[.8 .65 .125 .225]);

histogram(SIE_AMSR_NT-SIE_CDR,bins,'FaceColor','r','Normalization','pdf')
hold on
% histogram(SIE_ASI - SIE_CDR,bins,'FaceColor','m','Normalization','pdf')
histogram(SIE_AMSR_BS - SIE_CDR,bins,'Facecolor','b','Normalization','pdf')
xlim(xlimmer)
grid on; box on; 
view(90,-90)
title('$\Delta$ from CDR','interpreter','latex')
set(gca,'yticklabel','')

% 

Ax{3} = subplot('position',[.075 .3625 .65 .225]);
plot(xax,SIA_CDR,'k','linewidth',1); 
hold on
plot(xax,SIA_AMSR_NT,'r','linewidth',1)
plot(xax,SIA_AMSR_BS,'--b','linewidth',1); 
% plot(xax,SIA_ASI,'-m','linewidth',1); 


hold off
grid on; box on; 
title('Sea Ice Area','interpreter','latex')
ylabel('Million km$^2$','interpreter','latex')
set(gca,'xticklabel','')

% legend('CDR','AMSR2','AMSR2-BS','location',[.875 .8 .075 .1])

Ax{4} = subplot('position',[.8 .3625 .125 .225]);

histogram(SIA_AMSR_NT - SIA_CDR,bins,'FaceColor','r','Normalization','pdf')
hold on
histogram(SIA_AMSR_BS - SIA_CDR,bins,'Facecolor','b','Normalization','pdf')
% histogram(SIA_ASI - SIA_CDR,bins,'Facecolor','m','Normalization','pdf')
xlim(xlimmer)
% 
% xline(median(SIA_AMSR_NT - SIA_CDR),'r','linewidth',2)
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
plot(xax,AMIZ_AMSR_NT,'r','linewidth',1)
plot(xax,AMIZ_AMSR_BS,'--b','linewidth',1); 
% plot(xax,AMIZ_ASI,'--m','linewidth',1); 
hold off
grid on; box on; 
title('MIZ Extent','interpreter','latex')
ylabel('Million km$^2$','interpreter','latex')

% yline(median(AMIZ_AMSR_BS-AMIZ_CDR),'b','linewidth',2,'label',sprintf('%2.2f',median(AMIZ_CDR - AMIZ_AMSR_BS)),'fontsize',12)
% 

% legend('CDR','AMSR2','AMSR2-BS','location',[.875 .8 .075 .1])

Ax{6} = subplot('position',[.8 .075 .125 .225]);

histogram(AMIZ_AMSR_NT-AMIZ_CDR,bins,'FaceColor','r','Normalization','pdf')
hold on
histogram(AMIZ_AMSR_BS-AMIZ_CDR,bins,'Facecolor','b','Normalization','pdf')
% histogram(AMIZ_ASI-AMIZ_CDR,bins,'Facecolor','m','Normalization','pdf')
xlim(xlimmer)
view(90,-90)
hold off
grid on; box on; 
title('$\Delta$ from CDR','interpreter','latex')
set(gca,'yticklabel','')

letter = {'(A)','(B)','(C)','(D)','(E)','(F)','(g)','(e)','(c)'};

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
print('~/Dropbox (Brown)/Apps/Overleaf/IS2-PM-SIC/Figures/SIE-SIA-comp','-dpdf','-r600');

%%

clear Ax

climmer = [-.1 .3];

SICbins = linspace(0,1,51);
Bincent = 0.5*(SICbins(1:end-1) + SICbins(2:end));


sic_1 = AMSR_NT_SH(ICE_AMSR_BS);
sic_2 = AMSR_BS_SH(ICE_AMSR_BS);


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


% fitval = fit(Bincent(meanval_del > 0)',meanval_del(meanval_del > 0),myfit); 
% 
% 
% 

% Want to fit evenly over the range. 



%%
% myfit = fittype("-a*(c-d)^2 + b",...
%             dependent = "y", ...
%             coefficients=["a","b","d"],...
%             independent="c");
% % 
% usesic = sic_2 > 0.15 & ~(isnan(sic_2) | isnan(sic_1)); 
% 
% [nper,~,cat] = histcounts(sic_2(usesic),[0:.05:1]);
% weight = (sum(nper)./nper)/100;
% 
% % Doesn't matter since these don't appear. 
% weight(isinf(weight)) = 1; 
% 
% histweights = nan(size(sic_1(usesic)));
% 
% for i = 1:max(cat)
%     histweights(cat==i) = weight(i); 
% end
% 
% [fitval2,gof,fitinfo] = fit(sic_2(usesic),sic_1(usesic) - sic_2(usesic),myfit,'Startpoint',[1.62,.229,.6151],'weights',histweights); 
% 
%%

% This is fitted to data in proportion to representation in SIC_2 space. 
fitted = @(c) -1.639*(c-.6122).^2 + 0.2316; 

% This is fitting to all data without replacing. 
% fitted = @(c) -1.604*(c-.6138).^2 + 0.229; 


plotter = AMSR_NT_SH - AMSR_BS_SH; 
plotter(ICE_AMSR_NT == 0) = nan; 
plotter(ICE_AMSR_BS == 0) = nan; 

plotter_MIZ = AMSR_NT_SH - AMSR_BS_SH; 
plotter_MIZ(MIZ_AMSR_BS ~= 1) = nan; 
plotter_MIZ(ICE_AMSR_NT == 0) = nan; 
plotter_MIZ(ICE_AMSR_BS == 0) = nan; 

plotter_MIZ(repmat(sum(MIZ_AMSR_NT,3) < 6,[1 1 size(ICE_AMSR_BS,3)])) = nan; 

%%


close 

figure(1); 
clf; 

Ax{1} = subplot('position',[0 .5 .25 .45]);

worldmap([-90 -55],[-180 180]); 
pcolorm(lat_SH,lon_SH,median(plotter,3,'omitnan'))
make_HR_coastlines([.6 .6 .6]);
set(gca,'clim',climmer)
% set(gca,'clim',[-.2 .2])
colorbar('position',[.55 .55 .025 .35]); 
title('Bias: All','interpreter','latex')

Ax{2} = subplot('position',[.275 .5 .25 .45]); 

worldmap([-90 -55],[-180 180]); 
pcolorm(lat_SH,lon_SH,median(plotter_MIZ,3,'omitnan'))
make_HR_coastlines([.6 .6 .6]);
set(gca,'clim',climmer)
% set(gca,'clim',[-.2 .2])
% colorbar
title('Bias: Bootstrap MIZ','interpreter','latex')

nc = 33; 
cmapper = brewermap(nc,'PuOr');
cmapper(1:11,:) = [];
colormap(cmapper); 

Ax{3} = subplot('position',[.625 .5 .25 .45]); 
worldmap([-90 -55],[-180 180]); 

pcolorm(lat_SH,lon_SH,iqr(plotter_MIZ,3)./median(plotter_MIZ,3,'omitnan'))
make_HR_coastlines([.6 .6 .6]);

set(gca,'clim',[0 2])
colorbar('position',[.9 .55 .025 .35]); 
title('IQR/Bias','interpreter','latex')
colormap(gca,brewermap(11,'paired'));
%

Ax{4} =subplot('position',[.075 .1 .375 .35]);

jbfill(Bincent,plot_up',plot_dn',[.4 .4 .4],[1 1 1],1,.3);
hold on
plot(Bincent,meanval,'-k','linewidth',2)
xline(.8,'linewidth',1,'label','MIZ')
yline(.8,'linewidth',1,'label','MIZ')
line([0 1],[0 1],'Color',[.3 .3 .3],'linestyle','--')
title('SIC Comparison','interpreter','latex');

xlim([.15 1]);
ylim([.15 1]);

ylabel('AMSR2 SIC','interpreter','latex');
xlabel('AMSR2-BS SIC','interpreter','latex')
grid on; box on; 

yyaxis right
set(gca,'ycolor','k','yticklabel','')
bar(Bincent,a/sum(a),1,'FaceColor',[.4 .4 .8],'FaceAlpha',.5,'EdgeColor','none'); 

xlim([.15 1]);
% ylim([.15 1]);

Ax{5} = subplot('position',[.575 .1 .375 .35]);
cla

jbfill(Bincent,plot_up_del',plot_dn_del',[.4 .4 .4],[1 1 1],1,.3);
hold on
plot(Bincent,meanval_del,'-k','linewidth',2)
plot(Bincent(Bincent > 0.15),fitted(Bincent(Bincent > 0.15)),'--b','linewidth',1);


xline(.8,'linewidth',1,'label','MIZ');
yline(0,'--','linewidth',1)
xlim([.15 1]);
ylim([-.3 .4]);

% line([0 1],[0 1],'Color',[1 .4 .4],'linestyle','--')

plot(Bincent,1 - Bincent,'--r','linewidth',1)
plot(Bincent,-Bincent,'--r','linewidth',1)

ylabel('Bias','interpreter','latex');

xlabel('AMSR2-BS SIC','interpreter','latex')
title('Algorithmic Bias','interpreter','latex');
grid on; box on; 
%%

letter = {'(A)','(B)','(C)','(D)','(E)','(F)','(g)','(e)','(c)'};

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
print('~/Dropbox (Brown)/Apps/Overleaf/IS2-PM-SIC/Figures/Bias-BS','-dpdf','-r600');

% 