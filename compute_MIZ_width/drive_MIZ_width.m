% Main driver script for computing MIZ width from ICESat-2 data.
% This script sets up the analysis parameters, processes the data to generate
% along-track statistics, and then references these statistics to the MIZ.

clear

% Define the base path for data access. Update this to your local machine.
% OS_string = '/Users/chorvat/Library/CloudStorage/Dropbox-Brown/Christopher Horvat/';
OS_string = '/Users/chorvat/Brown Dropbox/Christopher Horvat/';

% Input folder containing the ICESat-2 track data.
OPTS.track_folder = '/gpfs/data/epscor/chorvat/IS2/Data/All_Track_Data/v6/SH/';

% Updatable name of input
OPTS.output_folder = '/Users/chorvat/Brown Dropbox/Christopher Horvat/Research Projects/Active/Data/ICESat-2/PM-SIC-width/Along_Track_Statistics/'

% Filename for the output data.
savename = 'AT_stats_SH_v6_all';

% Construct the full path for the output file.
OPTS.output_str = [OPTS.output_folder savename];

% Display the configuration settings.
disp('-----')
disp('Creating Along-Track Output Files')
fprintf('Using IS2 tracks located at: %s \n',OPTS.track_folder);
fprintf('Saving to %s \n',OPTS.output_folder);
fprintf('Filename: %s \n',savename);
disp('-----')

%% Setup analysis options and paths.

% Location of the analysis code.
OPTS.code_folder = '~/Code/compare_MIZ_width/';
% Location of utility functions.
OPTS.utils_folder = [OPTS.code_folder 'Utilities'];

% Analysis options.
OPTS.do_weak = 1; % Flag to include weak beams in the analysis (1 = yes, 0 = no).
OPTS.AT_window = [6250 6250]; % Along-track window size in meters.
OPTS.AT_resolution = 6250; % Resolution of the along-track data in meters.
OPTS.beamnames = {'/gt1r','/gt1l','/gt2r','/gt2l','/gt3r','/gt3l'}; % Names of the beams to process.

%% Locate ICESat-2 data files.

% Get a list of all .nc and .h5 files in the track folder.
OPTS.filenames = cat(1,dir([OPTS.track_folder '*.nc']),dir([OPTS.track_folder '*.h5']));
OPTS.nfiles = length(OPTS.filenames);

% Add the utilities folder to the MATLAB path.
addpath(OPTS.utils_folder);

%% Step 1: Create the along-track statistics file.
% This function processes each ICESat-2 track to compute downscaled statistics.
create_AT_stat_file(OPTS);

%% Step 2: Reference the statistics to the MIZ.
% This function adjusts the along-track statistics by referencing them to the
% MIZ, accounting for ascending and descending track components.

% The output from the previous step is used as input for this step.
OPTS.load_string = OPTS.output_str;

% Reference the computed statistics to the MIZ.
reference_stats_to_MIZ(OPTS);
