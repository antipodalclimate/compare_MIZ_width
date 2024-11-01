
clear

SSMI_loc = '/Users/chorvat/Brown Dropbox/Christopher Horvat/Research Projects/Active/Data/SIC-Data/NSIDC-CDR/Daily/NSIDC-CDR_daily.mat'; 
ASI_loc = '/Users/chorvat/Brown Dropbox/Christopher Horvat/Research Projects/Active/Data/SIC-Data/AMSR2-ASI/AMSR2_ASI_daily.mat'; 

load(SSMI_loc,'lat_SH','lon_SH');

%% Do with KNN

lander = ncread('NSIDC-0780_SeaIceRegions_PS-S25km_v1.0.nc','sea_ice_region_NASA_surface_mask'); 
is_land = lander > 5; 

%%

earthellipsoid = referenceSphere('earth','km');
lldist = @(x,y) distance(x,y,earthellipsoid);

landpoints = [lat_SH(is_land),lon_SH(is_land)]; 

M = createns(landpoints,'Distance',lldist);

%%
[ID,dist_to_coast] = knnsearch(M,[lat_SH(:),lon_SH(:)],'K',1,'Distance',lldist);

dist_to_coast = reshape(dist_to_coast,size(lat_SH)); 
save('dist_to_coast','dist_to_coast'); 

%% Now with ASI 

load(ASI_loc,"lat_ASI_SH",'lon_ASI_SH');
[ID,dist_to_coast_ASI] = knnsearch(M,[lat_ASI_SH(:),lon_ASI_SH(:)],'K',1,'Distance',lldist);
dist_to_coast_ASI = reshape(dist_to_coast_ASI,size(lat_ASI_SH)); 
save('dist_to_coast','dist_to_coast_ASI','-append'); 