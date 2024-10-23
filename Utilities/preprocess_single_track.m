function IS2_obj = preprocess_single_track(fieldname,beamname)
%% Load files

height = h5read(fieldname,[beamname '/sea_ice_segments/heights/height_segment_height']);
seg_len = h5read(fieldname,[beamname '/sea_ice_segments/heights/height_segment_length_seg']);
lat = h5read(fieldname,[beamname '/sea_ice_segments/latitude']);
lon = h5read(fieldname,[beamname '/sea_ice_segments/longitude']);
photon_rate = h5read(fieldname,[beamname '/sea_ice_segments/stats/photon_rate']); %    'freeboard_beam_segment/beam_freeboard/beam_fb_height']);
is_ice = h5read(fieldname,[beamname '/sea_ice_segments/heights/height_segment_type']);
ssh_flag = h5read(fieldname,[beamname '/sea_ice_segments/heights/height_segment_ssh_flag']);

if ~contains(fieldname,'_006_')
    v6 = 0; 
    conc = h5read(fieldname,[beamname '/sea_ice_segments/stats/ice_conc']);
else
    v6 = 1; 
    conc = h5read(fieldname,[beamname '/sea_ice_segments/stats/ice_conc_ssmi']);
    conc_amsr = h5read(fieldname,[beamname '/sea_ice_segments/stats/ice_conc_amsr2']);
end


% quality_flag = h5read(fieldname,[beamname '/sea_ice_segments/heights/height_segment_fit_quality_flag']);
% quality = h5read(fieldname,[beamname '/sea_ice_segments/heights/height_segment_quality']);

[~,filestr,~] = fileparts(fieldname); 
time_str = h5readatt(fieldname,'/','time_coverage_start');
time_str = datenum(time_str(1:10)); 

% if mean(diff(lat)) > 0
%     height = flipud(height); 
%     seg_len = flipud(seg_len); 
%     lat = flipud(lat);
%     lon = flipud(lon); 
%     photon_rate = flipud(photon_rate);
%     is_ice = flipud(is_ice); 
%     ssh_flag = flipud(ssh_flag); 
%     conc = flipud(conc); 
%     exmax_1 = flipud(exmax_1); 
%     exmax_2 = flipud(exmax_2); 
%     quality_flag = flipud(quality_flag); 
%     quality = flipud(quality); 
% end

%% Remove places with local low surface heights

earthellipsoid = referenceSphere('earth','m');
% First compute distance along track using lat and lon coordinates

%% Preprocess The track

if length(lat) > 1 % along-track distance
    dist = distance([lat(1:end-1) lon(1:end-1)],[lat(2:end) lon(2:end)],earthellipsoid);
else
    dist = [];
end

if ~isempty(dist)
    dist = [0; cumsum(dist)];
end


%% Preprocess The track

% Now remove unusable values along-track
unusable = find(abs(height > 1000 | seg_len > 2000));

dist(unusable) = [];

% Find duplicate values
dupes = find(diff(dist)<0.5)+1; % Find the duplicate points
dist(dupes) = []; % Those with such a distance get cut

[dist,b] = sort(dist); % Sort the distance to be increasing.

if ~isempty(dist)

dist = dist - dist(1); % Start distance at zero. 

end

% Dedupe and sort ice vector
is_ice(unusable) = [];
is_ice(dupes) = [];
is_ice = is_ice(b);

% Ocean is the stuff that isn't ice.
is_ocean = is_ice > 1;
is_ice = is_ice == 1; 

% Dedupe and sort height vector
height(unusable) = [];
height(dupes) = [];
height = height(b);

% Dedupe and sort segment length vector
seg_len(unusable) = [];
seg_len(dupes) = [];
seg_len = seg_len(b);

% Dedupe and sort concentration vector
conc(unusable) = [];
conc(dupes) = [];
conc = conc(b);

if v6
% Dedupe and sort concentration vector
conc_amsr(unusable) = [];
conc_amsr(dupes) = [];
conc_amsr = conc_amsr(b);
end

lat(unusable) = [];
lat(dupes) = [];
lat = lat(b);

lon(unusable) = [];
lon(dupes) = [];
lon = lon(b);

photon_rate(unusable) = [];
photon_rate(dupes) = [];
photon_rate = photon_rate(b);

ssh_flag(unusable) = [];
ssh_flag(dupes) = [];
ssh_flag = ssh_flag(b);

% exmax_1(unusable) = [];
% exmax_1(dupes) = [];
% exmax_1 = exmax_1(b);
% 
% exmax_2(unusable) = [];
% exmax_2(dupes) = [];
% exmax_2 = exmax_2(b);

% quality(unusable) = [];
% quality(dupes) = [];
% quality = quality(b);
% 
% quality_flag(unusable) = [];
% quality_flag(dupes) = [];
% quality_flag = quality_flag(b);

%% Now put into structure
IS2_obj.lat = lat; 
IS2_obj.lon = lon; 
IS2_obj.dist = dist; 
IS2_obj.height = height; 
IS2_obj.seg_len = seg_len; 
IS2_obj.photon_rate = photon_rate; 
IS2_obj.is_ice = is_ice; 
IS2_obj.is_ocean = is_ocean; 
IS2_obj.ssh_flag = ssh_flag; 
IS2_obj.conc = conc; 

if v6
    IS2_obj.conc_amsr = conc_amsr; 
end

% IS2_obj.exmax_1 = exmax_1; 
% IS2_obj.exmax_2 = exmax_2; 
% IS2_obj.quality = quality; 
% IS2_obj.quality_flag = quality_flag; 
IS2_obj.timer = time_str; 

IS2_obj.v6 = v6; 

end