% compare sentinel and IS2
clear

OS_string = '/Users/chorvat/Brown Dropbox/Christopher Horvat/';


addpath([OS_string 'Research Projects/Plot-Tools/'])
addpath([OS_string 'Research Projects/Plot-Tools/NE_Coastlines/'])

S1_fold = [OS_string 'Research Projects/Active/Data/Sentinel-1/Tavri_Classified/'];
S1_file_list = dir([S1_fold '*.mat']);

IS2_fold = [OS_string 'Research Projects/Active/Data/ICESat-2/PM-SIC-width/Kimia_Tracks/'];

% Location of this code.
OPTS.code_folder = '~/Code/compare_MIZ_width/';
% Location of utility files, output directory, and tracks.
OPTS.utils_folder = [OPTS.code_folder 'Utilities']; % Location of util files


OPTS.do_weak = 1; % Use the weak beams
OPTS.AT_window = [6250 6250];
OPTS.AT_resolution = 6250; % The resolution of the along-track data
OPTS.beamnames = {'/gt1r','/gt1l','/gt2r','/gt2l','/gt3r','/gt3l'};


% Pick an S1 files

for sarind = 1:length(S1_file_list)



    S1_file = S1_file_list(sarind)

    %%

    opts = delimitedTextImportOptions("NumVariables", 1);
    opts.Delimiter = ",";
    overlapping = readtable([OS_string 'Research Projects/Active/Data/Sentinel-1/2018-2023-overlaps-B1.txt'],opts);

    Sarnames = overlapping(:,1).Var1(1:2:end);
    Sarnames = cell2mat(Sarnames);
    IS2names = overlapping(:,1).Var1(2:2:end);

    overlapind = find(startsWith(string(Sarnames),S1_file.name(1:30)));

    IS2_track = IS2names(overlapind);

    %%

    trackname = [IS2_fold IS2_track{:}];

    do_downsample = 1;
    do_fullsample = 1;

    for j = 1:6 % for either all beams or just strong


        % Is this strong or weak beam
        beamtype = char(h5readatt(trackname,OPTS.beamnames{j},'atlas_beam_type'));
        is_strong(j) = strcmp(beamtype(1:3),'str');


        % Process this if we don't care about the beam
        % Or if we do care about the
        if OPTS.do_weak == 1 || is_strong(j)

            % This pulls out the data from each track. It orders the data
            % and removes poor-quality segments from the ATL07 data.
            IS2_obj = preprocess_single_track(trackname,OPTS.beamnames{j});

            % This now calculates along-track statistics by binning to the
            % appropriate along-track resolution
            [DS_stats{j},AT_stats{j}] = generate_AT_statistics(IS2_obj,OPTS.AT_window,OPTS.AT_resolution,do_downsample,do_fullsample);

        end

    end

    %%


    S1_data = load([S1_fold S1_file.name]);

    class_KT = S1_data.class(:,:,1);
    lat_KT = double(S1_data.lat);
    lon_KT = double(S1_data.lon);

    % class_KT = flipud(fliplr(S1_data.class(:,:,1)));
    % lat_KT = double(flipud(fliplr(S1_data.lat)));
    % lon_KT = double(flipud(fliplr(S1_data.lon)));

    %%

    unusable = class_KT == 255 | lat_KT == 0;

    class_KT(unusable) = nan;
    lon_KT(unusable) = nan;
    lat_KT(unusable) = nan;



    %%

    dx = 0.05; 
    xw = (1 - dx*sarind)/length(S1_file_list); 

    subplot('position',[dx*sarind + xw*(sarind-1) 0.4 xw 0.55])

    latlim = [min(lat_KT(:)) max(lat_KT(:))];
    lonlim = [min(lon_KT(:)) max(lon_KT(:))];

    worldmap(latlim,lonlim);
    pcolorm(lat_KT,lon_KT,class_KT);
    make_HR_coastlines([.8 .8 .8]);
    colormap(brewermap(12,'-Greys'))

    for i = 1:6
        plotm(AT_stats{i}.lat,AT_stats{i}.lon,'--');
    end

    %% Add IS2 things
    subplot('position',[dx*sarind + xw*(sarind-1) 0.05 xw 0.3])
    hold on


    for i = 1:6

    plot(DS_stats{i}.D_to_edge,DS_stats{i}.LIF,'b')
    plot(DS_stats{i}.D_to_edge,DS_stats{i}.SIC,'--k')
    plot(DS_stats{i}.D_to_edge,DS_stats{i}.SIC_amsr,'--r')

    end

    hold off

end