function create_coast_dist_data

OPTS.code_folder = '/Users/chorvat/Code/compare_MIZ_width';

SIC_data_folder = fullfile(OPTS.code_folder,'Data','SIC_Data');

SSMI_str = fullfile(SIC_data_folder,'NSIDC-CDR','CDR_daily_SH_v4.mat');
ASI_str = fullfile(SIC_data_folder,'AMSR2-ASI','AMSR2_ASI_daily.mat');
OSI_str = fullfile(SIC_data_folder,'OSI-450','OSISAF_daily.mat');

%% Do with KNN

lander = ncread('NSIDC-0780_SeaIceRegions_PS-S25km_v1.0.nc','sea_ice_region_NASA_surface_mask'); 
is_land = lander > 5; 

%%
load(SSMI_str,'lat_SH','lon_SH');

% Make the KNN searcher
landpoints = [lat_SH(is_land),lon_SH(is_land)]; 

earthellipsoid = referenceSphere('earth','km');
lldist = @(x,y) distance(x,y,earthellipsoid);

M = createns(landpoints,'Distance',lldist);

%% Now the dist to coast for PS NSIDC Grid
[ID,dist_to_coast] = knnsearch(M,[lat_SH(:),lon_SH(:)],'K',1,'Distance',lldist);
dist_to_coast = reshape(dist_to_coast,size(lat_SH)); 
save('dist_to_coast','dist_to_coast'); 

%% Now with ASI 

% load(ASI_str,"lat_ASI_SH",'lon_ASI_SH');
% [ID,dist_to_coast_ASI] = knnsearch(M,[lat_ASI_SH(:),lon_ASI_SH(:)],'K',1,'Distance',lldist);
% dist_to_coast_ASI = reshape(dist_to_coast_ASI,size(lat_ASI_SH)); 
% save('dist_to_coast','dist_to_coast_ASI','-append'); 

%% Now with OSI 

load(OSI_str,"lat_OSI_SH",'lon_OSI_SH');
[ID,dist_to_coast_OSI] = knnsearch(M,[lat_OSI_SH(:),lon_OSI_SH(:)],'K',1,'Distance',lldist);
dist_to_coast_OSI = reshape(dist_to_coast_OSI,size(lat_OSI_SH)); 
save('dist_to_coast','dist_to_coast_OSI','-append'); 

end