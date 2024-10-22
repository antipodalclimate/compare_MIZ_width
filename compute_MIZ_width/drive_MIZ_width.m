clear

OPTS.data_folder = '/Users/chorvat/Brown Dropbox/Christopher Horvat/Research Projects/Active/Data/ICESat-2/';
OPTS.code_folder = '/Users/chorvat/Code/compare_MIZ_width/'; 
OPTS.utils_folder = [OPTS.code_folder 'Utilities']; % Location of util files
OPTS.output_folder = [OPTS.data_folder 'PM-SIC-width/Brouwer-Tracks/Compiled_Stats/']; 

OPTS.do_weak = 0; % Don't use the weak beams
OPTS.AT_resolution = 6250; % The resolution of the along-track data
OPTS.beamnames = {'/gt1r','/gt1l','/gt2r','/gt2l','/gt3r','/gt3l'};

%% Necessary components for the location of IS2 data

OPTS.folder = [OPTS.data_folder 'PM-SIC-width/Brouwer-Tracks/Raw_Tracks_v6/'];
OPTS.filenames = cat(1,dir([OPTS.folder '*.nc']),dir([OPTS.folder '*.h5']));
OPTS.nfiles = length(OPTS.filenames);

% Where do we put the output
OPTS.output_str = [OPTS.output_folder 'AT_stats_SH_remake_v5'];  

addpath(OPTS.utils_folder); 

%% This creates a file which has the along-track downscaled data to the specified resolution
create_AT_stat_file(OPTS);

%% Next we realize that each track has both an ascending and descending component
% This needs to be adjusted since each track crosses the pole and therefore
% nominally intersects the MIZ once, on the descending and ascending side. 

OPTS.load_string = OPTS.output_str; 

reference_stats_to_MIZ(OPTS);

%% Second Figure

addpath('Figure-2-Brouwer/')
addpath('Utilities/Plotting/')

make_figure_brouwer_comp(OPTS)

%
% make_figure_brouwer_split;
