AMSR_str = fullfile(SIC_data_folder,'AMSR2-NT','AMSR2_SIC_daily.mat');
SSMI_str = fullfile(SIC_data_folder,'NSIDC-CDR','CDR_daily_SH.mat');

% ASI_loc = '/Users/chorvat/Brown Dropbox/Christopher Horvat/Research Projects/Active/Data/SIC-Data/AMSR2-ASI/AMSR2_ASI_daily.mat';

load(AMSR_str,'AMSR_datenum','AMSR_NT_SH','AMSR_BS_SH');
load(SSMI_str,'CDR_daily_SH','lat_SH','lon_SH','CDR_time_SH','BS_daily_SH','NT_daily_SH','area_SH','CDR_std_daily_SH');
% load(ASI_loc,'AMSR_ASI_SH','ASI_time_SH','lat_ASI_SH','lon_ASI_SH','area_ASI_SH');


%%
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

datenum_coincedent = AMSR_datenum(ia(ic));

%% Structure Data into SIC_Data

% Define products to process
% Format: {Name, Data Variable, Time Indices, Is_AMSR (for >1 check)}
products = {
    'AMSR_NT', AMSR_NT_SH, ia(ic), true;
    'AMSR_BS', AMSR_BS_SH, ia(ic), true;
    'CDR',     CDR_daily_SH, ib(ic), false;
    'SSMI_NT', NT_daily_SH,  ib(ic), false;
    'SSMI_BS', BS_daily_SH,  ib(ic), false;
    'CDR_std', CDR_std_daily_SH, ib(ic), false;
};

SIC_Data = struct();

% Loop through products to process and store in structure
for i = 1:size(products, 1)
    name = products{i, 1};
    raw_data = products{i, 2};
    indices = products{i, 3};
    is_amsr = products{i, 4};

    % Apply time indexing and coast mask
    processed_data = bsxfun(@times, raw_data(:,:,indices), coastmask);

    % Apply specific AMSR threshold logic
    if is_amsr
        processed_data(processed_data > 1) = nan;
    end

    SIC_Data.(name) = processed_data;
end

% Clear original large variables to free memory
clear AMSR_NT_SH AMSR_BS_SH CDR_daily_SH NT_daily_SH BS_daily_SH CDR_std_daily_SH products

%% Calculate Statistics (PM_Stats)

PM_Stats = struct();
prefac = 1/1e6;

fields = fieldnames(SIC_Data);
% We only want to calculate stats for the main SIC products, not CDR_std
stats_fields = setdiff(fields, {'CDR_std'});

for i = 1:numel(stats_fields)
    name = stats_fields{i};
    data = SIC_Data.(name);

    % Masks
    mask_MIZ = data > 0.15 & data < 0.8;
    mask_ICE = data > 0.15;

    % Calculate Areas
    PM_Stats.AMIZ.(name) = prefac * squeeze(sum(bsxfun(@times, area_SH, mask_MIZ), [1 2], 'omitnan'));
    PM_Stats.SIA.(name)  = prefac * squeeze(sum(bsxfun(@times, area_SH, data), [1 2], 'omitnan'));
    PM_Stats.SIE.(name)  = prefac * squeeze(sum(bsxfun(@times, area_SH, mask_ICE), [1 2], 'omitnan'));
end


%% Calculate Differences
% Six products
% AMSR2 - NT
% AMSR2 - BS
% AMSR2 - ASI
% SSMI - CDR
% SSMI - BS
% SSMI - NT

% Interesected in official product differences
% AMSR2 - SSMI-CDR
dAMIZ = PM_Stats.AMIZ.AMSR_NT - PM_Stats.AMIZ.CDR;

% Then interested in same algorithms, different sensors
dAMIZ_BS = PM_Stats.AMIZ.AMSR_BS - PM_Stats.AMIZ.SSMI_BS;
dAMIZ_NT = PM_Stats.AMIZ.AMSR_NT - PM_Stats.AMIZ.SSMI_NT;

% Then interested in same sensor, different algorithms
dAMIZ_AMSR = PM_Stats.AMIZ.AMSR_NT - PM_Stats.AMIZ.AMSR_BS;
dAMIZ_SSMI = PM_Stats.AMIZ.SSMI_NT - PM_Stats.AMIZ.SSMI_BS;

% Now for shits, different sensor, different algos
dAMIZ_cross = dAMIZ_NT + dAMIZ_SSMI;

dSIE = PM_Stats.SIE.AMSR_NT - PM_Stats.SIE.CDR;

% Now SIA differences
dSIA = PM_Stats.SIA.AMSR_NT - PM_Stats.SIA.CDR;
dSIA_NT = PM_Stats.SIA.AMSR_NT - PM_Stats.SIA.SSMI_NT;

%% Now examine the bias offset between NT2 and BS. 

SICbins = linspace(0,1,51);
Bincent_SIC = 0.5*(SICbins(1:end-1) + SICbins(2:end));

% Use masks consistent with original logic: ICE_AMSR_BS > 0.15
% Note: ICE_AMSR_BS corresponds to SIC_Data.AMSR_BS > 0.15
mask_ICE_AMSR_BS = SIC_Data.AMSR_BS > 0.15;

sic_1 = SIC_Data.AMSR_NT(mask_ICE_AMSR_BS);
sic_2 = SIC_Data.AMSR_BS(mask_ICE_AMSR_BS);

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

% Using the processed CDR_std from SIC_Data
cdr_std = SIC_Data.CDR_std(mask_ICE_AMSR_BS);

meanval_std= accumarray(c,cdr_std,[length(SICbins)-1 1],@nanmedian);
plot_up_std = accumarray(c,cdr_std,[length(SICbins)-1 1],upval);
plot_dn_std = accumarray(c,cdr_std,[length(SICbins)-1 1],dnval);

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
