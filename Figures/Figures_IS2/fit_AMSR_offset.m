clear
close all

AMSR_loc = '/Users/chorvat/Brown Dropbox/Christopher Horvat/Research Projects/Active/Data/SIC-Data/AMSR2-NT/AMSR2_SIC_daily.mat'; 
SSMI_loc = '/Users/chorvat/Brown Dropbox/Christopher Horvat/Research Projects/Active/Data/SIC-Data/NSIDC-CDR/Daily/NSIDC-CDR_daily.mat'; 
% ASI_loc = '/Users/chorvat/Brown Dropbox/Christopher Horvat/Research Projects/Active/Data/SIC-Data/AMSR2-ASI/AMSR2_ASI_daily.mat'; 

load(AMSR_loc,'AMSR_datenum','AMSR_NT_SH','AMSR_BS_SH');
load(SSMI_loc,'lat_SH','lon_SH','CDR_time_SH');
% load(ASI_loc,'AMSR_ASI_SH','ASI_time_SH','lat_ASI_SH','lon_ASI_SH','area_ASI_SH');

load('dist_to_coast'); 

%% Ignore close-to-coast 

coastmask = 1*(dist_to_coast > 25); 
coastmask(coastmask == 0) = nan; 

%%

% Same dates in the IS2 period. 
IS2_datenum = datenum('01-Oct-2018'):datenum('31-Dec-2024'); 

[mutual_datenum,ia,ib] = intersect(AMSR_datenum,CDR_time_SH);
[IS2_mutual,ic,~] = intersect(mutual_datenum,IS2_datenum);
% [ASI_mutual,id,ie] = intersect(IS2_mutual,ASI_time_SH); 

% NASATEAM2 AMSR value
AMSR_NT_SH = bsxfun(@times,AMSR_NT_SH(:,:,ia(ic)),coastmask);

% BOOTSTRAP AMSR values
AMSR_BS_SH =  bsxfun(@times,AMSR_BS_SH(:,:,ia(ic)),coastmask);

AMSR_NT_SH(AMSR_NT_SH > 1) = nan; 
AMSR_BS_SH(AMSR_BS_SH > 1) = nan; 

%%
SICbins = linspace(0,1,51);
Bincent = 0.5*(SICbins(1:end-1) + SICbins(2:end));

ICE_AMSR_BS = AMSR_BS_SH > 0.15; 

sic_1 = AMSR_NT_SH(ICE_AMSR_BS);
sic_2 = AMSR_BS_SH(ICE_AMSR_BS);

[a,b,c] = histcounts(sic_2,SICbins);

upval = @(x) prctile(x,75); 
dnval = @(x) prctile(x,25); 

meanval = accumarray(c,sic_1,[length(SICbins)-1 1],@nanmedian);
plot_up = accumarray(c,sic_1,[length(SICbins)-1 1],upval);
plot_dn = accumarray(c,sic_1,[length(SICbins)-1 1],dnval);

%%
meanval_del= accumarray(c,sic_1 - sic_2,[length(SICbins)-1 1],@nanmedian);
plot_up_del = accumarray(c,sic_1 - sic_2,[length(SICbins)-1 1],upval);
plot_dn_del = accumarray(c,sic_1 - sic_2,[length(SICbins)-1 1],dnval);

overmean_del = accumarray(c,(sic_1 - sic_2)./(1-sic_2),[length(SICbins)-1 1],@nanmedian);
overmean_del(isinf(overmean_del)) = 1; 

%%
myfit = fittype("-a*(c-d)^2 + b",...
            dependent = "y", ...
            coefficients=["a","b","d"],...
            independent="c");
% 
usesic = sic_2 > 0.15 & ~(isnan(sic_2) | isnan(sic_1)); 

[nper,~,cat] = histcounts(sic_2(usesic),[0:.05:1]);
weight = (sum(nper)./nper)/100;


% Doesn't matter since these don't appear. 
weight(isinf(weight)) = 1; 

histweights = nan(size(sic_1(usesic)));

for i = 1:max(cat)
    histweights(cat==i) = weight(i); 
end

[fitval2,gof,fitinfo] = fit(sic_2(usesic),sic_1(usesic) - sic_2(usesic),myfit,'Startpoint',[1.62,.229,.6151],'weights',histweights); 
