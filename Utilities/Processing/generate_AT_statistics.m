function [DS_stats, AT_stats] = generate_AT_statistics(IS2_obj, AT_window, AT_resolution, do_downsample, do_fullsample)
% GENERATE_AT_STATISTICS Computes along-track statistics from preprocessed ICESat-2 data.
%
% This function takes a preprocessed ICESat-2 track object and calculates
% various along-track statistics, such as Linear Ice Fraction (LIF),
% Wave-Affected Fraction (WAF), and Floe Size Distribution (FSD) metrics.
% The statistics are computed over a moving window and can be downsampled
% to a coarser resolution.
%
% Inputs:
%   IS2_obj         - A struct with preprocessed ICESat-2 data.
%   AT_window       - The size of the moving window for statistics (meters).
%   AT_resolution   - The resolution for downsampled data (meters).
%   do_downsample   - Flag to perform downsampling (true/false).
%   do_fullsample   - Flag to return high-resolution stats (true/false).
%
% Outputs:
%   DS_stats        - A struct with downsampled along-track statistics.
%   AT_stats        - A struct with high-resolution along-track statistics.

%% Initialize and Extract Data
AT_stats = struct(); % For high-resolution along-track stats
DS_stats = struct(); % For downsampled stats

height = IS2_obj.height;
is_ice = IS2_obj.is_ice;
is_ocean = IS2_obj.is_ocean;
is_dark = IS2_obj.is_dark;
is_spec = IS2_obj.is_spec;
seg_len = IS2_obj.seg_len;
lat = IS2_obj.lat;
lon = IS2_obj.lon;
dist = IS2_obj.dist;

%% Calculate Distance to Ice Edge
% This handles tracks that may cross the pole.
revdist = max(dist) - dist;
D_to_edge = dist; % Default assumption

if length(lat) > 1
    if abs(max(lon) - min(lon)) > 90 % Pole-crossing track
        D_to_edge = min(dist, revdist);
        swap_ind = find(revdist < dist, 1);
        D_to_edge(swap_ind:end) = -D_to_edge(swap_ind:end);
    else % Non-pole-crossing track
        if abs(lat(end)) < abs(lat(1)) % Equatorward track
            D_to_edge = -revdist;
        end
    end
end

%% Extract Additional Data
timer = IS2_obj.timer;
conc = IS2_obj.conc;
if IS2_obj.v6
    conc_amsr = IS2_obj.conc_amsr;
end

%% Define Analysis Windows and Parameters
window_1k = 1000;  % 1 km window
window_10k = 10000; % 10 km window
max_seg_size = 200; % Max segment size to be considered valid

%% Estimate Sea Surface Height (SSH)
% Interpolate SSH from available open water points.
if sum(is_ocean) > 1
    ssh_interp = interp1(dist(is_ocean), height(is_ocean), dist);
else
    ssh_interp = nan * dist;
end

height_adjusted = height - ssh_interp; % Height relative to local SSH

%% Calculate Moving Statistics for Wave Analysis
height_moving_std = movstd(height, window_10k, 'samplepoints', dist);
moving_ssh_no = movsum(is_ocean, window_10k, 'samplepoints', dist);
ssh_moving_std = movstd(ssh_interp, window_10k, 'samplepoints', dist);

isneg = height_adjusted < 0;
moving_neg = movsum(isneg, window_1k, 'samplepoints', dist);
moving_pos = movsum(~isneg, window_1k, 'samplepoints', dist);

%% Identify Floes
% Find transitions from ocean (1) to ice (0) and vice-versa.
up = strfind([0, is_ocean'], [0 1]);
down = strfind([is_ocean', 0], [1 0]);

if ~isempty(up)
    toosmall = intersect(up, down);
    up = setxor(up, toosmall)';
    down = setxor(down, toosmall)';
    floelen = dist(down) - dist(up);
    floeind = round(0.5 * (down + up));
    floe_seglength = floelen ./ (down - up);
    floe_nsegs = down - up;
else
    floelen = [];
    floeind = [];
end

% Filter for usable floes based on length and number of segments.
if ~isempty(floeind)
    usable_floe = logical((floelen > 30) .* (floe_seglength < 100) .* (floe_nsegs >= 3));
    usable_floe([1, end]) = 0; % Exclude endpoints
    floeind = floeind(usable_floe);
    floe_length = floelen(usable_floe);
else
    floe_length = [];
end

%% Define Criteria for Data Inclusion and Wave-Affected Segments
not_too_long = seg_len < max_seg_size;
close_to_ssh = (moving_ssh_no >= 2);
wave_cutoff_height = max(max(ssh_moving_std, .1), max(height_moving_std, .1));

% Define different surface types
is_ice = is_ice == 1;
is_not_spec = is_ice | is_dark;
is_not_dark = is_ice | is_spec;

% Define criteria for including segments in statistics
is_real = ~isnan(ssh_interp);
is_real_conc = conc < 1e37;
close_to_positive = moving_pos >= 2;
is_included = logical(not_too_long .* is_ice .* close_to_positive .* is_real);
is_included_all = logical(not_too_long .* close_to_positive .* is_real);
is_included_conc = is_included_all .* is_real_conc;

% Define wave-affected segments
close_to_negative = moving_neg >= 2;
is_wave_candidate = logical(is_included .* close_to_ssh .* close_to_negative);
is_under_both_var = logical((height_adjusted < -wave_cutoff_height) .* is_wave_candidate);

%% Calculate High-Resolution Along-Track (AT) Statistics
% Determine which points have enough data density for reliable stats.
AT_N_all = movsum(is_included_all, AT_window, 'SamplePoints', dist);
use_AT = 1 * (AT_N_all > sum(AT_window) / 30);
use_AT(use_AT == 0) = nan;
use_AT(dist < AT_window(1)) = nan;
use_AT(abs(max(dist) - dist) < AT_window(2)) = nan;
AT_N = use_AT .* AT_N_all;

% Wave-Affected Fraction (WAF)
wave_area_frac_both = 2 * movsum(seg_len .* is_under_both_var, AT_window, 'samplepoints', dist) ./ movsum(seg_len .* is_ice, AT_window, 'samplepoints', dist);
AT_WAF = use_AT .* wave_area_frac_both;

% Linear Ice Fraction (LIF) and its variants
AT_LIF = use_AT .* movsum(seg_len .* is_included, AT_window, 'samplepoints', dist) ./ movsum(seg_len .* is_included_all, AT_window, 'samplepoints', dist);
AT_LIF_spec = use_AT .* movsum(seg_len .* is_not_spec .* is_included_all, AT_window, 'samplepoints', dist) ./ movsum(seg_len .* is_included_all, AT_window, 'samplepoints', dist);
AT_LIF_dark = use_AT .* movsum(seg_len .* is_not_dark .* is_included_all, AT_window, 'samplepoints', dist) ./ movsum(seg_len .* is_included_all, AT_window, 'samplepoints', dist);

% Sea Ice Concentration (SIC) from passive microwave data
AT_SIC = use_AT .* (1/100) .* movsum(seg_len .* conc .* is_included_conc, AT_window, 'samplepoints', dist) ./ movsum(seg_len .* is_included_conc, AT_window, 'samplepoints', dist);
if IS2_obj.v6
    AT_SIC_amsr = use_AT .* (1/100) .* movsum(seg_len .* conc_amsr .* is_included_all, AT_window, 'samplepoints', dist) ./ movsum(seg_len .* is_included_all, AT_window, 'samplepoints', dist);
end

% Floe Size Distribution (FSD) metrics
AT_RFSD = use_AT(floeind) .* movsum(floe_length.^3, AT_window, 'samplepoints', dist(floeind)) ./ movsum(floe_length.^2, AT_window, 'samplepoints', dist(floeind));
AT_MFSD = use_AT(floeind) .* movsum(floe_length, AT_window, 'samplepoints', dist(floeind)) ./ movsum(floe_length.^0, AT_window, 'samplepoints', dist(floeind));

% Along-track height variance (Energy)
AT_E = use_AT .* movsum(seg_len .* height_adjusted.^2 .* is_included_all, AT_window, 'samplepoints', dist) ./ movsum(seg_len .* is_included_all, AT_window, 'samplepoints', dist);

% Along-track averaged height
AT_height = use_AT .* movsum(height_adjusted .* seg_len .* is_included, AT_window, 'omitnan', 'samplepoints', dist) ./ movsum(seg_len .* is_included, AT_window, 'omitnan', 'samplepoints', dist);

%% Store High-Resolution Stats if Requested
if do_fullsample
    AT_stats.timer = timer;
    AT_stats.lat = lat;
    AT_stats.lon = lon;
    AT_stats.use_AT = use_AT;
    AT_stats.LIF = AT_LIF;
    AT_stats.WAF = AT_WAF;
    AT_stats.SIC = AT_SIC;
    AT_stats.FSD = AT_RFSD;
    AT_stats.N = AT_N;
    AT_stats.H = AT_height;
    if IS2_obj.v6
        AT_stats.SIC_amsr = AT_SIC_amsr;
    end
    AT_stats.height_adj = height_adjusted;
    AT_stats.height_adj(AT_stats.height_adj > 10) = nan;
end

%% Downsample Statistics to Coarser Grid if Requested
if do_downsample
    % Map high-resolution points to the coarser grid
    dist_ind = floor(dist / AT_resolution / 2) + 1;
    [~, ~, ind_mapper] = unique(dist_ind);
    floe_mapper = ind_mapper(floeind);

    DS_stats.timer = timer;

    % Use accumarray to bin the data
    DS_stats.N_strict = accumarray(ind_mapper, use_AT, [], @nansum);
    DS_stats.lat = accumarray(ind_mapper, use_AT .* lat, [], @nanmean);
    DS_stats.lon = accumarray(ind_mapper, use_AT .* lon, [], @nanmean);
    DS_stats.D_to_edge = accumarray(ind_mapper, use_AT .* D_to_edge, [], @nanmean);
    DS_stats.H = accumarray(ind_mapper, use_AT .* height_adjusted, [], @nanmean);
    DS_stats.H_var = accumarray(ind_mapper, use_AT .* height_adjusted, [], @nanstd);
    DS_stats.WAF = accumarray(ind_mapper, AT_WAF, [], @nanmean);
    DS_stats.LIF = accumarray(ind_mapper, AT_LIF, [], @nanmean);
    DS_stats.LIF_spec = accumarray(ind_mapper, AT_LIF_spec, [], @nanmean);
    DS_stats.LIF_dark = accumarray(ind_mapper, AT_LIF_dark, [], @nanmean);
    DS_stats.SIC = accumarray(ind_mapper, AT_SIC, [], @nanmean);
    if IS2_obj.v6
        DS_stats.SIC_amsr = accumarray(ind_mapper, AT_SIC_amsr, [], @nanmean);
    end
    DS_stats.E = accumarray(ind_mapper, AT_E, [], @nanmean);

    if ~isempty(floe_mapper)
        DS_stats.Nfloe = accumarray(floe_mapper, use_AT(floeind), [], @nansum);
        DS_stats.MFSD = accumarray(floe_mapper, AT_MFSD, [], @nanmean);
        DS_stats.RFSD = accumarray(floe_mapper, AT_RFSD, [], @nanmean);
    else
        DS_stats.Nfloe = nan(size(DS_stats.H));
        DS_stats.MFSD = nan(size(DS_stats.H));
        DS_stats.RFSD = nan(size(DS_stats.H));
    end
end

end

