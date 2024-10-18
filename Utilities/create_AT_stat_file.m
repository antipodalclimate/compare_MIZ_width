function create_AT_stat_file(OPTS)
% This code creates a folder which contains a list of IS2 tracks and
% relevant statistics calculated in successive along-track windows.

% The width of the windowing is set by OPTS.AT_resolution
% The location of the output is set at OPTS.output_str
% We also select weak beams with the optional flag OPTS.do_weak
% We take the along-track variables and downscale them to the appropriate
% resolution

% Location of relevant code
addpath(OPTS.code_folder)

% Each file has 6 beams
AT_stats = cell(OPTS.nfiles,6);

for i = 1:OPTS.nfiles % for each individual track

    for j = 1:6 % for either all beams or just strong

        trackname = [OPTS.filenames(i).folder '/' OPTS.filenames(i).name];

        % Is this strong or weak beam
        beamtype = h5readatt(trackname,OPTS.beamnames{j},'atlas_beam_type');
        IS2_DATA.is_strong(i,j) = strcmp(beamtype(1:3),'str');


        % Process this if we don't care about the beam
        % Or if we do care about the 
        if OPTS.do_weak == 1 || IS2_DATA.is_strong(i,j)

            % This pulls out the data from each track. It orders the data
            % and removes poor-quality segments from the ATL07 data. 
            IS2_obj = preprocess_single_track(trackname,OPTS.beamnames{j});

            % This now calculates along-track statistics by binning to the
            % appropriate along-track resolution
            IS2_DATA.AT_stats{i,j} = downscale_AT_statistics(IS2_obj,OPTS.AT_resolution);

        end



    end

end

IS2_DATA.namearray = string(vertcat(OPTS.filenames(:).name));

save(OPTS.output_str,'OPTS','IS2_DATA');
