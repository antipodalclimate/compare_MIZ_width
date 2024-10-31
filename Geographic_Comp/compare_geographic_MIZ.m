
clear
close all

AMSR_loc = '/Users/chorvat/Brown Dropbox/Christopher Horvat/Research Projects/Active/Data/SIC-Data/AMSR2-NT/AMSR2_SIC_daily.mat'; 
SSMI_loc = '/Users/chorvat/Brown Dropbox/Christopher Horvat/Research Projects/Active/Data/SIC-Data/NSIDC-CDR/Daily/NSIDC-CDR_daily.mat'; 

load(AMSR_loc,'AMSR_datenum','AMSR_NT_SH','AMSR_BS_SH');
load(SSMI_loc,'CDR_SIC_SH','lat_SH','lon_SH','CDR_time_SH','BS_SIC_SH','NT_SIC_SH','area_SH');

load('dist_to_coast'); 



%% Ignore close-to-coast 

coastmask = 1*(dist_to_coast > 25); 
coastmask(coastmask == 0) = nan; 

%%

% Same dates in the IS2 period. 
IS2_datenum = datenum('01-Oct-2018'):datenum('31-Dec-2024'); 

[mutual_datenum,ia,ib] = intersect(AMSR_datenum,CDR_time_SH);
[IS2_mutual,ic,~] = intersect(mutual_datenum,IS2_datenum);

% NASATEAM2 AMSR value
AMSR_NT_SH = bsxfun(@times,AMSR_NT_SH(:,:,ia(ic)),coastmask);
% CDR values
CDR_SIC_SH = bsxfun(@times,CDR_SIC_SH(:,:,ib(ic)),coastmask);
SSMI_BS_SH = bsxfun(@times,BS_SIC_SH(:,:,ib(ic)),coastmask);
SSMI_NT_SH = bsxfun(@times,NT_SIC_SH(:,:,ib(ic)),coastmask);

clear BS_SIC_SH NT_SIC_SH

% BOOTSTRAP AMSR values
AMSR_BS_SH =  bsxfun(@times,AMSR_BS_SH(:,:,ia(ic)),coastmask);

%%

AMSR_NT_SH(AMSR_NT_SH > 1) = nan; 
AMSR_BS_SH(AMSR_BS_SH > 1) = nan; 

MIZ_AMSR_NT = AMSR_NT_SH > 0.15 & AMSR_NT_SH < 0.8; 
MIZ_AMSR_BS = AMSR_BS_SH > 0.15 & AMSR_BS_SH < 0.8; 
MIZ_CDR = CDR_SIC_SH > 0.15 & CDR_SIC_SH < 0.8; 
MIZ_SSMI_NT = SSMI_NT_SH > 0.15 & SSMI_NT_SH < 0.8; 
MIZ_SSMI_BS = SSMI_BS_SH > 0.15 & SSMI_BS_SH < 0.8; 

ICE_AMSR_NT = AMSR_NT_SH > 0.15; 
ICE_AMSR_BS = AMSR_BS_SH > 0.15; 
ICE_CDR = CDR_SIC_SH > 0.15; 
ICE_SSMI_NT = SSMI_NT_SH > 0.15; 
ICE_SSMI_BS = SSMI_BS_SH > 0.15; 

prefac = 1/1e6; 

AMIZ_AMSR_BS = prefac*squeeze(sum(bsxfun(@times,area_SH,MIZ_AMSR_BS),[1 2],'omitnan'));
AMIZ_AMSR_NT = prefac*squeeze(sum(bsxfun(@times,area_SH,MIZ_AMSR_NT),[1 2],'omitnan'));
AMIZ_CDR = prefac*squeeze(sum(bsxfun(@times,area_SH,MIZ_CDR),[1 2],'omitnan')); 
AMIZ_SSMI_NT = prefac*squeeze(sum(bsxfun(@times,area_SH,MIZ_SSMI_NT),[1 2],'omitnan'));
AMIZ_SSMI_BS = prefac*squeeze(sum(bsxfun(@times,area_SH,MIZ_SSMI_BS),[1 2],'omitnan'));


SIA_AMSR_NT = prefac*squeeze(sum(bsxfun(@times,area_SH,AMSR_NT_SH),[1 2],'omitnan'));
SIA_AMSR_BS = prefac*squeeze(sum(bsxfun(@times,area_SH,AMSR_BS_SH),[1 2],'omitnan'));
SIA_CDR = prefac*squeeze(sum(bsxfun(@times,area_SH,CDR_SIC_SH),[1 2],'omitnan'));
SIA_SSMI_NT = prefac*squeeze(sum(bsxfun(@times,area_SH,SSMI_NT_SH),[1 2],'omitnan'));
SIA_SSMI_BS = prefac*squeeze(sum(bsxfun(@times,area_SH,SSMI_BS_SH),[1 2],'omitnan'));

SIE_AMSR_NT = prefac*squeeze(sum(bsxfun(@times,area_SH,ICE_AMSR_NT),[1 2],'omitnan'));
SIE_AMSR_BS = prefac*squeeze(sum(bsxfun(@times,area_SH,ICE_AMSR_BS),[1 2],'omitnan'));
SIE_CDR = prefac*squeeze(sum(bsxfun(@times,area_SH,ICE_CDR),[1 2],'omitnan')); 
SIE_SSMI_NT = prefac*squeeze(sum(bsxfun(@times,area_SH,ICE_SSMI_NT),[1 2],'omitnan'));
SIE_SSMI_BS = prefac*squeeze(sum(bsxfun(@times,area_SH,ICE_SSMI_BS),[1 2],'omitnan'));

% Five products
% AMSR2 - NT
% AMSR2 - BS
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

delvals = AMSR_NT_SH - AMSR_;
delvals_MIZ = delvals; 
delvals_CIZ = delvals; 

delvals_MIZ(CDR_SIC_SH > 0.8) = nan; 
delvals_CIZ(CDR_SIC_SH < 0.8) = nan; 

delvals_MIZ = reshape(delvals_MIZ,[],length(IS2_mutual)); 
delvals_CIZ = reshape(delvals_CIZ,[],length(IS2_mutual)); 

yr0 = 2018; 

[r,upper,lower] = deal([]);

for yrind = 1:3
    for moind = 1:12
        
    data = delvals_MIZ(:,month(xax) == moind & year(xax) == yrind + yr0); 
    data(isnan(data)) = []; 

    r(1,moind,yrind) = median(data); 
    upper(moind,yrind,1) = prctile(data,75);
    lower(moind,yrind,1)= prctile(data,25); 

    data = delvals_CIZ(:,month(xax) == moind & year(xax) == yrind + yr0); 
    data(isnan(data)) = []; 

    r(2,moind,yrind) = median(data); 
    upper(moind,yrind,2) = prctile(data,75);
    lower(moind,yrind,2)= prctile(data,25); 




    end
end

upper = reshape(upper,36,2);
lower = reshape(lower,36,2);

%%
close all
xlimmer = [-5 5];
bins = -5:.2:5; 
%%

subplot('position',[.075 .65 .65 .225])
plot(xax,SIE_CDR,'k','linewidth',1); 
hold on
plot(xax,SIE_AMSR_NT,'r','linewidth',1)
% plot(xax,SIE_SSMI_BS,'--b','linewidth',1); 
plot(xax,SIE_AMSR_BS,'b','linewidth',1); 

grid on; box on; 
title('Sea Ice Extent','interpreter','latex')
ylabel('Million km$^2$','interpreter','latex')
set(gca,'xticklabel','')

% yline(median(AMIZ_CDR-AMIZ_AMSR_NT),'--r','linewidth',1,'label',sprintf('%2.2f',median(AMIZ_CDR - AMIZ_AMSR_NT)),'fontsize',12)

legend('CDR','AMSR2','AMSR2-BS','location',[.3 .95 .4 .025],'orientation','horizontal')

subplot('position',[.8 .65 .125 .225])

histogram(SIE_AMSR_NT-SIE_CDR,bins,'FaceColor','r','Normalization','pdf')
hold on
% histogram(SIE_CDR - SIE_SSMI_BS,bins,'FaceColor','b','Normalization','pdf')
histogram(SIE_AMSR_BS - SIE_CDR,bins,'Facecolor','b','Normalization','pdf')
xlim(xlimmer)
grid on; box on; 
view(90,-90)
title('$\Delta$ from CDR','interpreter','latex')
set(gca,'yticklabel','')

subplot('position',[.075 .3625 .65 .225])
plot(xax,SIA_CDR,'k','linewidth',1); 
hold on
plot(xax,SIA_AMSR_NT,'r','linewidth',1)
plot(xax,SIA_AMSR_BS,'--b','linewidth',1); 


hold off
grid on; box on; 
title('Sea Ice Area','interpreter','latex')
ylabel('Million km$^2$','interpreter','latex')
set(gca,'xticklabel','')

% legend('CDR','AMSR2','AMSR2-BS','location',[.875 .8 .075 .1])

subplot('position',[.8 .3625 .125 .225])

histogram(SIA_AMSR_NT - SIA_CDR,bins,'FaceColor','r','Normalization','pdf')
hold on
histogram(SIA_AMSR_BS - SIA_CDR,bins,'Facecolor','b','Normalization','pdf')
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

subplot('position',[.075 .075 .65 .225])
plot(xax,AMIZ_CDR,'k','linewidth',1); 
hold on
plot(xax,AMIZ_AMSR_NT,'r','linewidth',1)
plot(xax,AMIZ_AMSR_BS,'--b','linewidth',1); 
hold off
grid on; box on; 
title('MIZ Extent','interpreter','latex')
ylabel('Million km$^2$','interpreter','latex')

% yline(median(AMIZ_AMSR_BS-AMIZ_CDR),'b','linewidth',2,'label',sprintf('%2.2f',median(AMIZ_CDR - AMIZ_AMSR_BS)),'fontsize',12)
% 

% legend('CDR','AMSR2','AMSR2-BS','location',[.875 .8 .075 .1])

subplot('position',[.8 .075 .125 .225])

histogram(AMIZ_AMSR_NT-AMIZ_CDR,bins,'FaceColor','r','Normalization','pdf')
hold on
histogram(AMIZ_AMSR_BS-AMIZ_CDR,bins,'Facecolor','b','Normalization','pdf')
xlim(xlimmer)
view(90,-90)
hold off
grid on; box on; 
title('$\Delta$ from CDR','interpreter','latex')
set(gca,'yticklabel','')


pos = [6.5 4.5]; 
set(gcf,'windowstyle','normal','position',[0 0 pos],'paperposition',[0 0 pos],'papersize',pos,'units','inches','paperunits','inches');
set(gcf,'windowstyle','normal','position',[0 0 pos],'paperposition',[0 0 pos],'papersize',pos,'units','inches','paperunits','inches');
print('/Users/chorvat/Brown Dropbox/Christopher Horvat/Apps/Overleaf/IS2-Waves-PM/Figures/SIE-SIA-comp','-dpdf','-r600');

% 
% %%
% jbfill(1:36,upper(:,1)',lower(:,1)',[.8 .4 .4]); 
% hold on
% 
% jbfill(1:36,upper(:,2)',lower(:,2)',[.4 .4 .8]); 
% hold off
% 
% %%
% figure(3);
% jbfill(1:12,upper,lower,[.5 .5 .5])
% 
% hold on
% plot(xax',r,'r')

%%


%%


close 


climmer = [-.1 .3];

% plotter = AMSR_NT_SH - AMSR_BS_SH; 
% % Eliminate Bootstrap MIZ
% plotter(ICE_AMSR_BS & MIZ_AMSR_BS) = nan; 
% plotter(ICE_AMSR_NT == 0) = nan; 
% plotter(ICE_AMSR_BS == 0) = nan; 
% 
% subplot('position',[0 .5 .3 .4]); 
% worldmap([-90 -55],[-180 180]); 
% pcolorm(lat_SH,lon_SH,mean(plotter,3,'omitnan'))
% make_HR_coastlines([.6 .6 .6]);
% set(gca,'clim',climmer)
% % colorbar
% title('Bias: Bootstrap CIZ','interpreter','latex')

% Eliminate Bootstrap CIZ

subplot('position',[0 .5 .275 .5]); 

plotter = AMSR_NT_SH - AMSR_BS_SH; 
plotter(ICE_AMSR_NT == 0) = nan; 
plotter(ICE_AMSR_BS == 0) = nan; 

worldmap([-90 -55],[-180 180]); 
pcolorm(lat_SH,lon_SH,median(plotter,3,'omitnan'))
make_HR_coastlines([.6 .6 .6]);
set(gca,'clim',climmer)
% set(gca,'clim',[-.2 .2])
colorbar('position',[.55 .55 .025 .4]); 
title('Bias: All','interpreter','latex')

subplot('position',[.275 .5 .25 .5]); 

plotter = AMSR_NT_SH - AMSR_BS_SH; 
plotter(MIZ_AMSR_BS ~= 1) = nan; 
plotter(ICE_AMSR_NT == 0) = nan; 
plotter(ICE_AMSR_BS == 0) = nan; 

plotter(repmat(sum(MIZ_AMSR_NT,3) < 6,[1 1 size(ICE_AMSR_BS,3)])) = nan; 

worldmap([-90 -55],[-180 180]); 
pcolorm(lat_SH,lon_SH,median(plotter,3,'omitnan'))
make_HR_coastlines([.6 .6 .6]);
set(gca,'clim',climmer)
% set(gca,'clim',[-.2 .2])
% colorbar
title('Bias: Bootstrap MIZ','interpreter','latex')

cmapper = cmocean('curl','pivot',0);

nc = 33; 
cmapper = brewermap(nc,'PuOr');
cmapper(1:11,:) = [];
colormap(cmapper); 

subplot('position',[.625 .5 .25 .5]); 
worldmap([-90 -55],[-180 180]); 

pcolorm(lat_SH,lon_SH,iqr(plotter,3)./median(plotter,3,'omitnan'))
make_HR_coastlines([.6 .6 .6]);

set(gca,'clim',[0 2])
colorbar('position',[.9 .55 .025 .4]); 
title('IQR/Bias','interpreter','latex')
colormap(gca,brewermap(11,'paired'));
%
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


subplot('position',[.075 .1 .4 .4])

jbfill(Bincent,plot_up',plot_dn',[.4 .4 .4],[1 1 1],1,.3);
hold on
plot(Bincent,meanval,'-k','linewidth',2)
xline(.8,'linewidth',1,'label','MIZ')
yline(.8,'linewidth',1,'label','MIZ')
line([0 1],[0 1],'Color',[.9 .9 .9],'linestyle','--')
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

meanval = accumarray(c,sic_1 - sic_2,[length(SICbins)-1 1],@nanmedian);
plot_up = accumarray(c,sic_1 - sic_2,[length(SICbins)-1 1],upval);
plot_dn = accumarray(c,sic_1 - sic_2,[length(SICbins)-1 1],dnval);

overmean = accumarray(c,(sic_1 - sic_2)./(1-sic_2),[length(SICbins)-1 1],@nanmedian);
overmean(isinf(overmean)) = 1; 

subplot('position',[.55 .1 .4 .4])
cla

jbfill(Bincent,plot_up',plot_dn',[.4 .4 .4],[1 1 1],1,.3);
hold on
plot(Bincent,meanval,'-k','linewidth',2)
xline(.8,'linewidth',1,'label','MIZ');
yline(0,'--k','linewidth',1)
xlim([.15 1]);
ylim([-.2 .3]);

line([0 1],[0 1],'Color',[1 .4 .4],'linestyle','--')

plot(Bincent,1 - Bincent,'--r','linewidth',1,'label','maximum')

ylabel('Bias','interpreter','latex');

xlabel('AMSR2-BS SIC','interpreter','latex')
title('Algorithmic Bias','interpreter','latex');
grid on; box on; 
%%
pos = [6.5 4.5]; 
set(gcf,'windowstyle','normal','position',[0 0 pos],'paperposition',[0 0 pos],'papersize',pos,'units','inches','paperunits','inches');
set(gcf,'windowstyle','normal','position',[0 0 pos],'paperposition',[0 0 pos],'papersize',pos,'units','inches','paperunits','inches');
print('/Users/chorvat/Brown Dropbox/Christopher Horvat/Apps/Overleaf/IS2-Waves-PM/Figures/Bias-BS','-dpdf','-r600');

% 