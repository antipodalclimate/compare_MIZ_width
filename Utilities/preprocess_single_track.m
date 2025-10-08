function IS2_obj = preprocess_single_track(fieldname, beamname)
% PREPROCESS_SINGLE_TRACK Loads and preprocesses a single ICESat-2 track.
%
% This function reads data from a single ICESat-2 file for a specified beam,
% cleans it by removing poor-quality or duplicate segments, and organizes
% the data into a structured object for further analysis.
%
% Inputs:
%   fieldname - The full path to the ICESat-2 data file (.h5 or .nc).
%   beamname  - The name of the beam to process (e.g., '/gt1r').
%
% Outputs:
%   IS2_obj   - A struct containing the preprocessed track data, including
%               latitude, longitude, distance, height, and other parameters.

%% Load data from HDF5 file
height = h5read(fieldname, [beamname '/sea_ice_segments/heights/height_segment_height']);
seg_len = h5read(fieldname, [beamname '/sea_ice_segments/heights/height_segment_length_seg']);
lat = h5read(fieldname, [beamname '/sea_ice_segments/latitude']);
lon = h5read(fieldname, [beamname '/sea_ice_segments/longitude']);
photon_rate = h5read(fieldname, [beamname '/sea_ice_segments/stats/photon_rate']);
is_ice = h5read(fieldname, [beamname '/sea_ice_segments/heights/height_segment_type']);
ssh_flag = h5read(fieldname, [beamname '/sea_ice_segments/heights/height_segment_ssh_flag']);

% Check data version and load appropriate sea ice concentration data.
if ~contains(fieldname, '_006_')
    v6 = 0;
    conc = h5read(fieldname, [beamname '/sea_ice_segments/stats/ice_conc']);
else
    v6 = 1;
    conc = h5read(fieldname, [beamname '/sea_ice_segments/stats/ice_conc_ssmi']);
    conc_amsr = h5read(fieldname, [beamname '/sea_ice_segments/stats/ice_conc_amsr2']);
end
conc(conc > 1e37) = nan; % Handle fill values.

% Get the start time of the data coverage.
time_str = char(h5readatt(fieldname, '/', 'time_coverage_start'));
time_str = datenum(time_str(1:10));

%% Preprocess the track data
earthellipsoid = referenceSphere('earth', 'm');

% Calculate along-track distance from lat/lon coordinates.
if length(lat) > 1
    dist = distance([lat(1:end-1) lon(1:end-1)], [lat(2:end) lon(2:end)], earthellipsoid);
    dist = [0; cumsum(dist)];
else
    dist = [];
end

if isempty(dist)
    IS2_obj = []; % Return empty if no distance data
    return;
end

% Identify and mark unusable data segments.
unusable = find(abs(height) > 1000 | seg_len > 2000 | is_ice <= 0);
dist(unusable) = [];

% Identify and mark duplicate data points based on distance.
dupes = find(diff(dist) < 0.5) + 1;
dist(dupes) = [];

% Sort the data by along-track distance.
[dist, b] = sort(dist);
dist = dist - dist(1); % Normalize distance to start at zero.

% --- Clean and sort all data vectors based on the valid indices ---
% Surface type
is_ice(unusable) = [];
is_ice(dupes) = [];
is_ice = is_ice(b);

% Classify surface types based on the 'is_ice' flag.
is_ocean = is_ice > 1;
is_dark = is_ice > 5;
is_spec = is_ice > 1 & is_ice < 6;
is_ice = is_ice == 1; % Final logical flag for sea ice.

% Height
height(unusable) = [];
height(dupes) = [];
height = height(b);

% Segment length
seg_len(unusable) = [];
seg_len(dupes) = [];
seg_len = seg_len(b);

% Sea ice concentration
conc(unusable) = [];
conc(dupes) = [];
conc = conc(b);
if v6
    conc_amsr(unusable) = [];
    conc_amsr(dupes) = [];
    conc_amsr = conc_amsr(b);
end

% Geolocation
lat(unusable) = [];
lat(dupes) = [];
lat = lat(b);
lon(unusable) = [];
lon(dupes) = [];
lon = lon(b);

% Other parameters
photon_rate(unusable) = [];
photon_rate(dupes) = [];
photon_rate = photon_rate(b);
ssh_flag(unusable) = [];
ssh_flag(dupes) = [];
ssh_flag = ssh_flag(b);

%% Assemble the output struct
IS2_obj.lat = lat;
IS2_obj.lon = lon;
IS2_obj.dist = dist;
IS2_obj.height = height;
IS2_obj.seg_len = seg_len;
IS2_obj.photon_rate = photon_rate;
IS2_obj.is_ice = is_ice;
IS2_obj.is_ocean = is_ocean;
IS2_obj.is_dark = is_dark;
IS2_obj.is_spec = is_spec;
IS2_obj.ssh_flag = ssh_flag;
IS2_obj.conc = conc;
IS2_obj.timer = time_str;
IS2_obj.v6 = v6;

if v6
    IS2_obj.conc_amsr = conc_amsr;
end

end
