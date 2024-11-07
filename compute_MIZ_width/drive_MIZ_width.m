clear

% Location of the overall data
OPTS.track_folder = '/gpfs/data/epscor/chorvat/IS2/Data/All_Track_Data/v6/SH/'

% Updatable name of input
OPTS.output_folder = '/gpfs/data/epscor/chorvat/IS2/Along_Track_Statistics/'

savename = 'AT_stats_SH_v6_all'; 

% Updatable name of output. 
OPTS.output_str = [OPTS.output_folder savename];

disp('-----')
disp('Creating Along-Track Output Files')
fprintf('Using IS2 tracks located at: %s \n',OPTS.track_folder); 
fprintf('Saving to %s \n',OPTS.output_folder); 
fprintf('Filename: %s \n',savename); 
disp('-----')


%%
% Location of this code. 
OPTS.code_folder = '~/Code/compare_MIZ_width/';
% Location of utility files, output directory, and tracks. 
OPTS.utils_folder = [OPTS.code_folder 'Utilities']; % Location of util files


OPTS.do_weak = 1; % Use the weak beams
OPTS.AT_window = [6250 6250];
OPTS.AT_resolution = 6250; % The resolution of the along-track data
OPTS.beamnames = {'/gt1r','/gt1l','/gt2r','/gt2l','/gt3r','/gt3l'};

%% Necessary components for the location of IS2 data

OPTS.filenames = cat(1,dir([OPTS.track_folder '*.nc']),dir([OPTS.track_folder '*.h5']));
OPTS.nfiles = length(OPTS.filenames);

% Where do we put the output

addpath(OPTS.utils_folder); 

%% This creates a file which has the along-track downscaled data to the specified resolution
create_AT_stat_file(OPTS);

%% Next we realize that each track has both an ascending and descending component
% This needs to be adjusted since each track crosses the pole and therefore
% nominally intersects the MIZ once, on the descending and ascending side. 

OPTS.load_string = OPTS.output_str; 

reference_stats_to_MIZ(OPTS);

%%
