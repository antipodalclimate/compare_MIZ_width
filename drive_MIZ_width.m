clear

OPTS.main_folder = pwd; % Location of main code
OPTS.code_folder = [OPTS.main_folder '/Utilities']; % Location of util files

OPTS.do_weak = 0; % Don't use the weak beams
OPTS.AT_resolution = 12500; % The resolution of the along-track data
OPTS.beamnames = {'/gt1r','/gt1l','/gt2r','/gt2l','/gt3r','/gt3l'};

%% Necessary components for the location of IS2 data

OPTS.folder = '/Users/chorvat/Library/CloudStorage/Dropbox-Brown/Christopher Horvat/Research Projects/Active/Data/ICESat-2/Brouwer-Tracks-All/';
OPTS.filenames = dir([OPTS.folder '*.h5']);
OPTS.nfiles = length(OPTS.filenames);

% Where do we put the output
OPTS.output_str = [OPTS.main_folder '/Output/AT_stats_Brouwer_test'];  

addpath(OPTS.code_folder); 

%% This creates a file which has the along-track downscaled data to the specified resolution
create_AT_stat_file(OPTS);

%% Next we realize that each track has both an ascending and descending component
% This needs to be adjusted since each track crosses the pole and therefore
% nominally intersects the MIZ once, on the descending and ascending side. 

OPTS.load_string = [OPTS.main_folder '/Output/AT_stats_Brouwer_remake']; 

reference_stats_to_MIZ(OPTS);

%% Second Figure

addpath('Figure-2-Brouwer/')
addpath('Utilities/Plotting/')

make_figure_brouwer_comp(OPTS)

%
% make_figure_brouwer_split;
