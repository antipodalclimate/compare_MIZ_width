% MAKE_DUMMY_DATA Generates a mock ICESat-2 HDF5 file for testing.
%
% This script creates a dummy .h5 file with the structure expected by the
% MIZ width analysis pipeline. It generates synthetic data for one beam
% to allow the pipeline to run through without errors.

clear
close all

%% Setup Paths
[current_path, ~, ~] = fileparts(mfilename('fullpath'));
root_path = fileparts(current_path);
data_folder = fullfile(root_path, 'Data', 'Track_Data');

if ~exist(data_folder, 'dir')
    mkdir(data_folder);
end

filename = fullfile(data_folder, 'dummy_track_data_006_01.h5');
fprintf('Creating dummy data file: %s\n', filename);

% Delete if exists
if exist(filename, 'file')
    delete(filename);
end

%% parameters
n_segments = 1000;
beam_names = {'gt1r', 'gt1l', 'gt2r', 'gt2l', 'gt3r', 'gt3l'};

%% Create File and Attributes
% Create the file (implicitly done by h5create, but we need group structure)
% We'll use h5write to write data, h5create to define datasets.

% Global Attribute: time_coverage_start
% We can't easily write global attributes with high-level functions in older Matlab,
% but h5writeatt works.
% Create a dummy dataset first to ensure file exists?
% Actually h5create creates the file.

% Define datasets
for b = 1:length(beam_names)
    beam = beam_names{b};
    base_path = ['/' beam];

    % Beam Type Attribute
    % We will write this later using h5writeatt

    % Define paths
    p_height = [base_path '/sea_ice_segments/heights/height_segment_height'];
    p_seg_len = [base_path '/sea_ice_segments/heights/height_segment_length_seg'];
    p_lat = [base_path '/sea_ice_segments/latitude'];
    p_lon = [base_path '/sea_ice_segments/longitude'];
    p_photon_rate = [base_path '/sea_ice_segments/stats/photon_rate'];
    p_is_ice = [base_path '/sea_ice_segments/heights/height_segment_type'];
    p_ssh_flag = [base_path '/sea_ice_segments/heights/height_segment_ssh_flag'];
    p_conc_ssmi = [base_path '/sea_ice_segments/stats/ice_conc_ssmi'];
    p_conc_amsr = [base_path '/sea_ice_segments/stats/ice_conc_amsr2'];

    % Create datasets
    h5create(filename, p_height, [n_segments 1]);
    h5create(filename, p_seg_len, [n_segments 1]);
    h5create(filename, p_lat, [n_segments 1]);
    h5create(filename, p_lon, [n_segments 1]);
    h5create(filename, p_photon_rate, [n_segments 1]);
    h5create(filename, p_is_ice, [n_segments 1], 'Datatype', 'int32');
    h5create(filename, p_ssh_flag, [n_segments 1], 'Datatype', 'int32');
    h5create(filename, p_conc_ssmi, [n_segments 1]);
    h5create(filename, p_conc_amsr, [n_segments 1]);

    % Generate Synthetic Data
    % Lat/Lon: Simulate a track moving from south to north or crossing pole
    % Let's do a simple track from -70 to -60 lat
    lat_data = linspace(-70, -60, n_segments)';
    lon_data = ones(n_segments, 1) * 180;

    % Height: Random noise around 0.5m
    height_data = 0.5 + 0.1 * randn(n_segments, 1);

    % Seg len: mostly 30m
    seg_len_data = 30 * ones(n_segments, 1);

    % Type: 1 = ice. Make some ocean (type > 1) at the beginning
    type_data = ones(n_segments, 1, 'int32');
    type_data(1:100) = 0; % invalid
    type_data(101:200) = 2; % Ocean? check preprocess logic: is_ocean = is_ice > 1; is_ice = is_ice == 1;

    % Conc: High concentration
    conc_data = 90 * ones(n_segments, 1);

    % Write Data
    h5write(filename, p_height, height_data);
    h5write(filename, p_seg_len, seg_len_data);
    h5write(filename, p_lat, lat_data);
    h5write(filename, p_lon, lon_data);
    h5write(filename, p_photon_rate, 10 * ones(n_segments, 1));
    h5write(filename, p_is_ice, type_data);
    h5write(filename, p_ssh_flag, zeros(n_segments, 1, 'int32'));
    h5write(filename, p_conc_ssmi, conc_data);
    h5write(filename, p_conc_amsr, conc_data);

    % Write Attributes
    % beam_type needs to be 'strong' or 'weak'
    % is_strong(i, j) = strcmp(beamtype(1:3), 'str');
    % gt1r is usually strong? Let's make all strong for simplicity
    h5writeatt(filename, ['/' beam], 'atlas_beam_type', 'strong');
end

% Global Attribute
h5writeatt(filename, '/', 'time_coverage_start', '2020-01-01T00:00:00');

disp('Dummy data generation complete.');
