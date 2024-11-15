clear

OS_string = '/Users/chorvat/Dropbox-Brown/Christopher Horvat/';
OS_string = '/Users/chorvat/Brown Dropbox/Christopher Horvat/';

OPTS.output_str = '/Users/chorvat/Brown Dropbox/Christopher Horvat/Research Projects/Active/Data/ICESat-2/PM-SIC-width/Along_Track_Statistics/AT_stats_SH_v6_all'; 
% OPTS.output_str = [OS_string 'Research Projects/Active/Data/ICESat-2/PM-SIC-width/Along_Track_Statistics/AT_stats_KM_v6_all']; 

load(OPTS.output_str);

cutoff_N = 100; % need more than 100 segments for analysis to make sense
cutoff_MIZ = 1; % Need more than 1 MIZ segment for analysis. 


% make_location_figure
load_MIZ_waves; 

%%
% Criteria for selection

usable_all = (Nsegvals > cutoff_N) & usable_all; 
usable_all = usable_all & SICvals > 0.15 & LIFvals > 0.15;
usable_all = usable_all & npoints > 1; 
usable_all = usable_all & abs(lonvals) < 180 & latvals > -90; 

time_indices = 12*(yearval-2018) + monthval; 


%%
usable = usable_all; %
usable_waves = usable & wavytracks > 0; 
usable_nowaves = usable & wavytracks == 0; 


lat_waves = latvals(usable_waves); 
lon_waves = lonvals(usable_waves);
timer_waves = time_indices(usable_waves); 


lat_nowaves = latvals(usable_nowaves); 
lon_nowaves = lonvals(usable_nowaves);
timer_nowaves = time_indices(usable_nowaves); 

%%
close all

% Put into a 1/4 degree grid. 

dlat = 2; 

lat_grid = -90:dlat:0; 
lon_grid = -180:dlat:180; 

[lat_even,lon_even] = meshgrid(lat_grid,lon_grid);
lat_even = lat_even'; 
lon_even = lon_even'; 

lat_ind = floor((1/dlat)*(lat_waves - min(lat_grid))) + 1;
lon_ind = floor((1/dlat)*(lon_waves - (min(lon_grid))))+1;
ID_waves = sub2ind(size(lat_even),lat_ind,lon_ind); 
nper_waves = accumarray(ID_waves,0+1*lat_waves,[numel(lat_even) 1],@sum); 

subplot('position',[.05 .05 .4 .5]); 

worldmap([-90 -50],[-180 180])
pcolorm(lat_even,lon_even,reshape(nper_waves,size(lat_even)))
make_HR_coastlines([.6 .6 .6]);

lat_ind = floor((1/dlat)*(lat_nowaves - min(lat_grid))) + 1;
lon_ind = floor((1/dlat)*(lon_nowaves - (min(lon_grid))))+1;


ID_nowaves = sub2ind(size(lat_even),lat_ind,lon_ind); 
nper_nowaves = accumarray(ID_nowaves,0+1*lat_nowaves,[numel(lat_even) 1],@sum); 

subplot('position',[.55 .05 .4 .5]); 

worldmap([-90 -50],[-180 180])
pcolorm(lat_even,lon_even,reshape(nper_nowaves,size(lat_even)))
make_HR_coastlines([.6 .6 .6]);
% Now do time plots

subplot('position',[.1 .7 .8 .25]);

nmonthly_waves = accumarray(timer_waves,1+0*timer_waves,[72 1],@sum);
nmonthly_nowaves = accumarray(timer_nowaves,1+0*timer_nowaves,[72 1],@sum);

plot(nmonthly_waves,'--k'); 
hold on
plot(nmonthly_nowaves,'--r')
legend('Waves','No Waves')
grid on; box on; 