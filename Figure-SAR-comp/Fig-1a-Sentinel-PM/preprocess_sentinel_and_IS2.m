clear

DBstring = '/Users/chorvat/Dropbox (Brown)/Research Projects/'; 

addpath([DBstring 'Plot-Tools/'])
addpath([DBstring 'Plot-Tools/NE_Coastlines/'])
% addpath('Figure-2-Brouwer/')
addpath('../Utilities/')

%%
S1_fold = [DBstring 'Active/Data/Sentinel-1/'];
S1_file = 'S1A_EW_GRDM_1SSH_20190225T014506_20190225T014611_026079_02E89C_56F1.nc';

S1_fold_KT = [DBstring 'Active/Data/Sentinel-1/Tavri_Classified/'];
S1_file_KT = 'ThurSAR__corrected_mask_with_latlon.tiff';

PM_fold = [DBstring 'Active/Data/SIC-Data/NSIDC-CDR/'];
PM_file = 'Daily/SIA_data/seaice_conc_daily_sh_20190224_f17_v04r00.nc'; 

IS2_fold = [DBstring 'Active/Data/ICESat-2/PM-SIC-width/SAR-comp-tracks/'];
IS2_file = 'ATL07-02_20190224012038_08800201_006_02.h5';

%% Pull in S1 image
lat_S1 = ncread([S1_fold S1_file],'lat'); 
lon_S1 = ncread([S1_fold S1_file],'lon'); 

S1_lat_span = [min(lat_S1(:)) max(lat_S1(:))];
S1_lon_span = [min(lon_S1(:)) max(lon_S1(:))];

Amp_S1 = ncread([S1_fold S1_file],'Amplitude_HH'); 
Int_S1 = ncread([S1_fold S1_file],'Intensity_HH'); 

%% Pull in classified S1 Data 
class_KT = imread([S1_fold_KT S1_file_KT]);
lat_KT = rot90(class_KT(:,:,2),2); 
lon_KT = rot90(class_KT(:,:,3),2);
class_KT = rot90(class_KT(:,:,1),2); 

%% Pull in passive microwave data
sic_PM = ncread([PM_fold PM_file],'cdr_seaice_conc');
lat_PM = load([PM_fold 'NSIDC-CDR_monthly'],'lat_SH').lat_SH; 
lon_PM = load([PM_fold 'NSIDC-CDR_monthly'],'lon_SH').lon_SH;

%% Pull in IS2 data. 
% Reference IS2 along-track values to the grids of the original S1 image,
% the classified image of KT, and the passive microwave grid. 
beamnames = {'/gt1r','/gt1l','/gt2r','/gt2l','/gt3r','/gt3l'};

% track_span = [0 150];
lat_span = [-70.9188  -69.5810];
lon_span = [-83.68 -83.1];

inds_S1 = lat_S1(:) < lat_span(2) & lat_S1(:) > lat_span(1) & lon_S1(:) < lon_span(2) & lon_S1(:) > lon_span(1);
inds_KT = lat_KT(:) < lat_span(2) & lat_KT(:) > lat_span(1) ... 
    & lon_KT(:) < lon_span(2) & lon_KT(:) > lon_span(1);
inds_PM = lat_PM(:) < lat_span(2) & lat_PM(:) > lat_span(1) ... 
    & lon_PM(:) < lon_span(2) & lon_PM(:) > lon_span(1);

lat_X_S1 = lat_S1(inds_S1);  
lon_X_S1 = lon_S1(inds_S1);  
S1_Amp_subset = Amp_S1(inds_S1); 

M_S1 = createns([lat_X_S1,lon_X_S1]); 

lat_X_KT = lat_KT(inds_KT);  
lon_X_KT = lon_KT(inds_KT);  
KT_subset = class_KT(inds_KT); 

M_KT = createns([lat_X_KT,lon_X_KT]); 

lat_X_PM = lat_PM(inds_PM);  
lon_X_PM = lon_PM(inds_PM);  
PM_subset = sic_PM(inds_PM); 

M_PM = createns([lat_X_PM,lon_X_PM]); 

%% Calculate the along-track stats from IS2 as well as the along-track values from PM and SAR
 
AT_window = [6250 6250]; 
AT_resolution = 6250; 

for i = 1:length(beamnames)

IS2_obj{i} = preprocess_single_track([IS2_fold IS2_file],beamnames{i});
[~,AT_stats{i}] = generate_AT_statistics(IS2_obj{i},AT_window,AT_resolution,0,1);

% Nearest neighbor of the IS2 lat/lon on the S1 lat/lon. 
ID_S1{i} = knnsearch(M_S1,[IS2_obj{i}.lat IS2_obj{i}.lon],'K',1); 
S1_lats{i} = lat_X_S1(ID_S1{i}); 
S1_lons{i} = lon_X_S1(ID_S1{i}); 
AT_amp_S1{i} = S1_Amp_subset(ID_S1{i}); 

% Now do the same for the classified KT image
ID_KT{i} = knnsearch(M_KT,[IS2_obj{i}.lat IS2_obj{i}.lon],'K',1); 
KT_lats{i} = lat_X_KT(ID_KT{i});
KT_lons{i} = lon_X_KT(ID_KT{i});
AT_class_KT{i} = KT_subset(ID_KT{i}); 


ID_PM{i} = knnsearch(M_PM,[IS2_obj{i}.lat IS2_obj{i}.lon],'K',1); 
PM_lats{i} = lat_X_PM(ID_PM{i});
PM_lons{i} = lon_X_PM(ID_PM{i});
AT_sic_PM{i} = PM_subset(ID_PM{i}); 

end


%% Do along-track smoothing of non-IS2 fields. 
% All others should be using the same AT_window as specified above. 
for i = 1:length(beamnames)

AT_S1_LIF{i} = AT_stats{i}.use_AT.*movmean(AT_class_KT{i},AT_window,'omitnan','samplepoints',IS2_obj{i}.dist);
AT_PM_SIC{i} = AT_stats{i}.use_AT.*movmean(AT_sic_PM{i},AT_window,'omitnan','samplepoints',IS2_obj{i}.dist);

% height_smooth{i} = movsum(AT_stats{i}.height_adj(IS2_obj{i}.is_ice).*IS2_obj{i}.seg_len(IS2_obj{i}.is_ice),AT_window,'omitnan','samplepoints',IS2_obj{i}.dist(IS2_obj{i}.is_ice)) ...
%     ./ movsum(IS2_obj{i}.seg_len(IS2_obj{i}.is_ice),AT_window,'omitnan','samplepoints',IS2_obj{i}.dist(IS2_obj{i}.is_ice));
% 
% height_var{i} = movsum((AT_stats{i}.height_adj(IS2_obj{i}.is_ice) - height_smooth{i}).^2 .* IS2_obj{i}.seg_len(IS2_obj{i}.is_ice),AT_window,'omitnan','samplepoints',IS2_obj{i}.dist(IS2_obj{i}.is_ice)) ...
%     ./ movsum(IS2_obj{i}.seg_len(IS2_obj{i}.is_ice),AT_window,'omitnan','samplepoints',IS2_obj{i}.dist(IS2_obj{i}.is_ice));
% 
% height_std{i} = sqrt(height_var{i}); 
% 
% 
end

%% Now compute the distance between MIZ and not. 

% Compact ice zone start
[X_CIZ_WAF,X_CIZ_LIF,X_CIZ_PM,X_CIZ_S1,X_CIZ_ISPM,X_CIZ_AMSR] = deal(nan(1,length(beamnames)));
[lat_CIZ_WAF,lat_CIZ_LIF,lat_CIZ_PM,lat_CIZ_S1,lat_CIZ_ISPM,lat_CIZ_AMSR] = deal(nan(1,length(beamnames)));
[lon_CIZ_WAF,lon_CIZ_LIF,lon_CIZ_PM,lon_CIZ_S1,lon_CIZ_ISPM,lon_CIZ_AMSR] = deal(nan(1,length(beamnames)));

% Marginal ice zone start
[X_MIZ_WAF,X_MIZ_LIF,X_MIZ_PM,X_MIZ_S1,X_MIZ_ISPM,X_MIZ_AMSR] = deal(nan(1,length(beamnames)));
[lat_MIZ_WAF,lat_MIZ_LIF,lat_MIZ_PM,lat_MIZ_S1,lat_MIZ_ISPM,lat_MIZ_AMSR] = deal(nan(1,length(beamnames)));
[lon_MIZ_WAF,lon_MIZ_LIF,lon_MIZ_PM,lon_MIZ_S1,lon_MIZ_ISPM,lon_MIZ_AMSR] = deal(nan(1,length(beamnames)));

for i = 1:length(beamnames)
%% 

% Take the first point where WAF drops below 0.075
ind_MIZ_WAF = find(AT_stats{i}.WAF > 0.075,1);
ind_CIZ_WAF = find(AT_stats{i}.WAF < 0.075,1);

if ~isempty(ind_MIZ_WAF)
    X_MIZ_WAF(i)= IS2_obj{i}.dist(ind_MIZ_WAF)/1000;
    lat_MIZ_WAF(i)= IS2_obj{i}.lat(ind_MIZ_WAF);
    lon_MIZ_WAF(i)= IS2_obj{i}.lon(ind_MIZ_WAF);
end

if ~isempty(ind_CIZ_WAF)
    X_CIZ_WAF(i)= IS2_obj{i}.dist(ind_CIZ_WAF)/1000;
    lat_CIZ_WAF(i)= IS2_obj{i}.lat(ind_CIZ_WAF);
    lon_CIZ_WAF(i)= IS2_obj{i}.lon(ind_CIZ_WAF);
end

% Take the first point where IS2 LIF reaches above 15% or 80%
ind_MIZ_LIF = find(AT_stats{i}.LIF > 0.15,1);
ind_CIZ_LIF = find(AT_stats{i}.LIF > 0.8,1);

if ~isempty(ind_MIZ_LIF)
    X_MIZ_LIF(i)= IS2_obj{i}.dist(ind_MIZ_LIF)/1000;
    lat_MIZ_LIF(i)= IS2_obj{i}.lat(ind_MIZ_LIF);
    lon_MIZ_LIF(i)= IS2_obj{i}.lon(ind_MIZ_LIF);
end

if ~isempty(ind_CIZ_LIF)
    X_CIZ_LIF(i)= IS2_obj{i}.dist(ind_CIZ_LIF)/1000;
    lat_CIZ_LIF(i)= IS2_obj{i}.lat(ind_CIZ_LIF);
    lon_CIZ_LIF(i)= IS2_obj{i}.lon(ind_CIZ_LIF);
end

% Take the first point where S1 LIF reaches above 15% or 80%
ind_MIZ_S1 = find(AT_S1_LIF{i} > 0.15,1);
ind_CIZ_S1 = find(AT_S1_LIF{i} > 0.80,1);

if ~isempty(ind_MIZ_S1)
    X_MIZ_S1(i)= IS2_obj{i}.dist(ind_MIZ_S1)/1000;
    lat_MIZ_S1(i)= IS2_obj{i}.lat(ind_MIZ_S1);
    lon_MIZ_S1(i)= IS2_obj{i}.lon(ind_MIZ_S1);
end

if ~isempty(ind_CIZ_S1)
    X_CIZ_S1(i)= IS2_obj{i}.dist(ind_CIZ_S1)/1000;
    lat_CIZ_S1(i)= IS2_obj{i}.lat(ind_CIZ_S1);
    lon_CIZ_S1(i)= IS2_obj{i}.lon(ind_CIZ_S1);
end

% Take the first point where PM reaches above 15% or 80%
ind_MIZ_PM = find(AT_sic_PM{i} > 0.15,1);
ind_CIZ_PM = find(AT_sic_PM{i} > 0.80,1);

if ~isempty(ind_MIZ_PM)
    X_MIZ_PM(i)= IS2_obj{i}.dist(ind_MIZ_PM)/1000;
    lat_MIZ_PM(i)= IS2_obj{i}.lat(ind_MIZ_PM);
    lon_MIZ_PM(i)= IS2_obj{i}.lon(ind_MIZ_PM);
end

if ~isempty(ind_CIZ_PM)
    X_CIZ_PM(i)= IS2_obj{i}.dist(ind_CIZ_PM)/1000;
    lat_CIZ_PM(i)= IS2_obj{i}.lat(ind_CIZ_PM);
    lon_CIZ_PM(i)= IS2_obj{i}.lon(ind_CIZ_PM);
end

% Take the first point where PM reaches above 15% or 80%
ind_MIZ_ISPM = find(AT_stats{i}.SIC > 0.15,1);
ind_CIZ_ISPM = find(AT_stats{i}.SIC > 0.80,1);

if ~isempty(ind_MIZ_ISPM)
    X_MIZ_ISPM(i)= IS2_obj{i}.dist(ind_MIZ_ISPM)/1000;
    lat_MIZ_ISPM(i)= IS2_obj{i}.lat(ind_MIZ_ISPM);
    lon_MIZ_ISPM(i)= IS2_obj{i}.lon(ind_MIZ_ISPM);
end

if ~isempty(ind_CIZ_ISPM)
    X_CIZ_ISPM(i)= IS2_obj{i}.dist(ind_CIZ_ISPM)/1000;
    lat_CIZ_ISPM(i)= IS2_obj{i}.lat(ind_CIZ_ISPM);
    lon_CIZ_ISPM(i)= IS2_obj{i}.lon(ind_CIZ_ISPM);
end

% Take the first point where AMSR2 reaches above 15% or 80%
ind_MIZ_AMSR = find(AT_stats{i}.SIC_amsr > 0.15,1);
ind_CIZ_AMSR = find(AT_stats{i}.SIC_amsr > 0.80,1);

if ~isempty(ind_MIZ_AMSR)
    X_MIZ_AMSR(i)= IS2_obj{i}.dist(ind_MIZ_AMSR)/1000;
    lat_MIZ_AMSR(i)= IS2_obj{i}.lat(ind_MIZ_AMSR);
    lon_MIZ_AMSR(i)= IS2_obj{i}.lon(ind_MIZ_AMSR);
end

if ~isempty(ind_CIZ_AMSR)
    X_CIZ_AMSR(i)= IS2_obj{i}.dist(ind_CIZ_AMSR)/1000;
    lat_CIZ_AMSR(i)= IS2_obj{i}.lat(ind_CIZ_AMSR);
    lon_CIZ_AMSR(i)= IS2_obj{i}.lon(ind_CIZ_AMSR);
end


    
end

W_PM = X_CIZ_PM - X_MIZ_PM; 
W_S1 = X_CIZ_S1 - X_MIZ_S1; 
W_LIF = X_CIZ_LIF - X_MIZ_LIF; 
W_WAF = X_CIZ_WAF - X_MIZ_WAF; 
W_ISPM = X_CIZ_ISPM - X_MIZ_ISPM; 
W_AMSR = X_CIZ_AMSR - X_MIZ_AMSR; 

%%

