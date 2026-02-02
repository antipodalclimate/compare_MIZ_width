function create_AT_stat_file(OPTS)
% CREATE_AT_STAT_FILE Processes ICESat-2 tracks to generate along-track statistics.
%
% This function loops through a list of ICESat-2 data files, preprocesses
% each track to extract and clean the data, and then computes along-track
% statistics. The results are saved to a .mat file.
%
% Inputs:
%   OPTS - A struct containing options and paths for the analysis:
%     .track_folder   - Path to the folder with ICESat-2 track files.
%     .output_str     - Full path for the output .mat file.
%     .do_weak        - Flag to process weak beams (1 = yes, 0 = no).
%     .AT_window      - Along-track window size.
%     .AT_resolution  - Resolution for along-track data.
%     .beamnames      - Cell array of beam names to process.
%     .filenames      - List of track file names.
%     .nfiles         - Number of track files.
%     .code_folder    - Path to the main code directory.

% Add the code folder to the MATLAB path.
addpath(OPTS.code_folder);

% Initialize cell arrays to store statistics for each file and beam.
DS_stats = cell(OPTS.nfiles, 6);
AT_stats = DS_stats;
is_strong = nan(OPTS.nfiles, 6);

% Set flags for data sampling.
do_downsample = true;  % Downsample the data to the specified resolution.
do_fullsample = false; % Do not return all per-segment data.

% Loop through each ICESat-2 track file.
for i = 1:OPTS.nfiles
    if mod(i, 100) == 1
        fprintf('Processing file %d of %d \n', i, OPTS.nfiles);
    end

    % Loop through each beam in the track.
    for j = 1:6
        trackname = [OPTS.filenames(i).folder '/' OPTS.filenames(i).name];

        % Determine if the beam is strong or weak.
        beamtype = char(h5readatt(trackname, OPTS.beamnames{j}, 'atlas_beam_type'));
        is_strong(i, j) = strcmp(beamtype(1:3), 'str');

        % Process the beam if it's strong or if weak beams are included.
        if OPTS.do_weak == 1 || is_strong(i, j)
            % Preprocess the track to extract and clean data.
            IS2_obj = preprocess_single_track(trackname, OPTS.beamnames{j});

            % Generate along-track statistics.
            [DS_stats{i, j}, AT_stats{i, j}] = generate_AT_statistics(IS2_obj, OPTS.AT_window, OPTS.AT_resolution, do_downsample, do_fullsample);
        end
    end
end

% Check the version of the ICESat-2 data.
fieldname = [OPTS.filenames(1).folder '/' OPTS.filenames(1).name];
if ~contains(fieldname, '_006_')
    IS2_DATA.v6 = 0;
else
    IS2_DATA.v6 = 1;
end

% Consolidate results into a single struct.
IS2_DATA.AT_stats = AT_stats;
IS2_DATA.DS_stats = DS_stats;
IS2_DATA.is_strong = is_strong;
IS2_DATA.namearray = string(vertcat(OPTS.filenames(:).name));

% Save the processed data to a .mat file.
save(OPTS.output_str, 'OPTS', 'IS2_DATA', '-v7.3');
