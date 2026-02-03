AMSR_loc = '/Users/chorvat/Dropbox (Brown)/Research Projects/Active/Data/SIC-Data/AMSR2-NT/AMSR2_SIC_daily.mat';
SSMI_loc = '/Users/chorvat/Dropbox (Brown)/Research Projects/Active/Data/SIC-Data/NSIDC-CDR/Daily/NSIDC-CDR_daily.mat';

% ASI_loc = '/Users/chorvat/Brown Dropbox/Christopher Horvat/Research Projects/Active/Data/SIC-Data/AMSR2-ASI/AMSR2_ASI_daily.mat';

load(AMSR_loc,'AMSR_datenum','AMSR_NT_SH','AMSR_BS_SH');
load(SSMI_loc,'CDR_SIC_SH','lat_SH','lon_SH','CDR_time_SH','BS_SIC_SH','NT_SIC_SH','area_SH');
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

%% Now examine the bias offset between NT2 and BS. 

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




