function [DS_stats,AT_stats] = generate_AT_statistics(IS2_obj,AT_window,AT_resolution,do_downsample,do_fullsample)
%% Pull out details from the object

% One is simple along-track statistics
AT_stats = struct();
% Second is along-track statistics downsampled to a coarser grid.
DS_stats = struct();




height = IS2_obj.height;
is_ice = IS2_obj.is_ice;
is_ocean = IS2_obj.is_ocean;
is_dark = IS2_obj.is_dark;
is_spec = IS2_obj.is_spec;
seg_len = IS2_obj.seg_len;
lat = IS2_obj.lat;
lon = IS2_obj.lon;

%%
% This is the distance from the first point, and from the last point in the
% along-track distance.
dist = IS2_obj.dist;
revdist = max(dist) - dist;

% Initially, we assume first track is edge. Distance is then how far from
% that first ice point.
D_to_edge = dist;

% D_to_edge should be positive if moving inside the edge in the direction of the track.
% It will then be negative if moving towards the edge in the direction of
% the track.

% But this doesn't have to be the outer edge. The track can be
% ascending in lat or descending. It can also cross the continent.

% If the track crosses the entire pole, we need to consider
% that the distance from the edge switches sign. The first point will be
% one edge, the last point will be the last edge.

if length(lat) > 1

    if abs(max(lon) - min(lon)) > 90

        % Minimum of difference from either side.
        D_to_edge = min(dist,revdist);
        % At a certain point, closer to end than beginning.
        swap_ind = find(revdist < dist,1);
        D_to_edge(swap_ind:end) = -D_to_edge(swap_ind:end);

        % If the track does not cross the entire continent, we need to check to be
        % sure that the edge is not actually the last point.

    else

        if abs(lat(end)) < abs(lat(1))
            % equatorially-going track. First latitude is actually closer to the
            % pole
            % This means that the final point is actually the "edge" closer to the
            % equator
            D_to_edge = -revdist;

        else
            % This is a pole-going track that doesn't cross substantially over
            % the pole. First point is the most equatorward, therefore the
            % distance in the edge should be positive increasing
        end

    end

end

%%
timer = IS2_obj.timer;
conc = IS2_obj.conc;

if IS2_obj.v6

    conc_amsr = IS2_obj.conc_amsr;

end

%%
% Identify local moving average
window_1k = 1000;
window_10k = 10000; % meters - size of moving window
window_25k = 25000; % meters - size of moving window

slide_25k = [12500 12500];

% slide_10k = [0 10000];
% slide_50k = [0 50000];
%

max_seg_size = 200; % Maximum size for individual segments
% earthellipsoid = referenceSphere('earth','m');
% exmax_diff_threshold = 1; % meters - difference in two surfaces tracked by the ex-max algorithm.

%% Now worry about SSH points

% This is the index of the nearest neighbor to each point that has open
% water
if isempty(dist(is_ocean))

    dist_to_ssh = nan*dist;

else

    [~,dist_to_ssh] = knnsearch(dist(is_ocean),dist);
end

if sum(is_ocean) > 1

    % This is the interpolated SSH field, on all of the dist places
    ssh_interp = interp1(dist(is_ocean),height(is_ocean),dist);


else

    ssh_interp = nan*dist;

end

%% Objects interpolated on the 50km moving window

% Moving average wave energy
% Number of SSH points within that window

% moving_en = movmean((height-ssh_interp).^2 .* seg_len,window_50k,'samplepoints',dist);


%% Objects interpolated on the 10km moving window
% Standard deviation of Height

% height_moving_avg = cat(1,height_moving_avg,movmean(height,window_variance,'samplepoints',dist));
height_moving_std = movstd(height,window_10k,'samplepoints',dist);

moving_ssh_no = movsum(is_ocean,window_10k,'samplepoints',dist);

ssh_moving_std = movstd(ssh_interp,window_10k,'samplepoints',dist);


%% Objects interpolated on the 1km moving window

% Negative values of adjusted SSH
isneg = height - ssh_interp < 0;

moving_neg = movsum(isneg,window_1k,'samplepoints',dist);
moving_pos = movsum(~isneg,window_1k,'samplepoints',dist);

%% FSD-related things
up = strfind([0,is_ocean'],[0 1]);
down = strfind([is_ocean',0],[1 0]);

if ~isempty(up)

    toosmall = intersect(up,down);
    up = setxor(up,toosmall)';
    down = setxor(down,toosmall)';

    Uloc = up - 1;
    Uloc(Uloc==0) = 1;

    floelen = dist(down) - dist(up);
    floeind = round(.5*(down + up));
    floe_seglength = floelen./(down - up);
    floe_nsegs = down - up;


else

    floelen = [];
    floeind = [];
    floe_seglength = [];
    floe_nsegs = [];

end

%%

naive_floe = (floe_nsegs >= 2);
naive_floe(1) = 0;
naive_floe(end) = 0;

usable_floe = logical((floelen > 30).* (floe_seglength < 100).*(floe_nsegs >= 3));
usable_floe(1) = 0; % exclude endpoints
usable_floe(end) = 0; % exclude endpoints

% Now look at whether - at each naive location, there is a usable
% floe.

goodfloes = usable_floe(naive_floe);

% Naive - need at least 3 segments
floelen_0 = floelen(naive_floe);
floeind_0 = floeind(naive_floe);
floeind = floeind(usable_floe);
floe_length = floelen(usable_floe);

floe_seglength_0 = floe_seglength(naive_floe);
floe_seglength = floe_seglength(usable_floe);

%% Wave data
% Multiple for WAL
% Mval = (.5 - (1/pi)*asin(ssh_interp./moving_en + 1/sqrt(2))).^(-1);
% Mval(abs(imag(Mval)) > 0) = 0;
% Mval(isnan(Mval)) = 0;
% Mval = min(Mval,10);

% Compute the fraction of all measurements (by length or number) that are
% "wave-affected"

% First eliminate all segments larger than 1 km - may be artificially big
not_too_long =seg_len<max_seg_size;

% Need to have two ssh points within 10 km
close_to_ssh = (moving_ssh_no >= 2);%  .* (ssh_neighbors <= window_10k);

% Need to know that the multi-surfaces are near each other.
% has_single_surface = abs(exmax_1 - exmax_2) ...
%     < exmax_diff_threshold;

% Minimum threshold on moving std
wave_cutoff_ssh = max(ssh_moving_std,.1);
wave_cutoff_height = max(height_moving_std,.1);
both_cutoff_height = max(wave_cutoff_ssh,wave_cutoff_height);

% ice if tagged as sea ice by ATL07
% surf_type = is_ice;


is_ice = is_ice == 1;
is_not_spec = is_ice | is_dark;
is_not_dark = is_ice | is_spec;

% adjust for deviation from local ssh
height_adjusted = height - ssh_interp;

% Not too long, is ice, and has positive moving average
% naive_reasonable = logical(not_too_long .* is_ice);

% Has positive points nearby
close_to_positive = moving_pos >= 2;
% Has negative points nearby - only for those values that are negative
close_to_negative = moving_neg >= 2;

% Also require that there is an interpolated ssh field. May not always be
% true especially at endpoints.
is_real = ~isnan(ssh_interp);

is_real_conc = conc < 1e37; 

% Included points have positive points nearby, aren't too long, and are
% identified as ice. Criteria I1-I2.
is_included = logical(not_too_long .* is_ice .* close_to_positive .* is_real);
is_included_all = logical(not_too_long .* close_to_positive .* is_real);
is_included_conc = is_included_all .* is_real_conc; 

% Included points have positive points nearby, and aren't too long. Criteria I1.
% is_included_lead = logical(not_too_long .* close_to_positive);

%%
% Wave candidates are those i|ncluded segs close to ssh points, and
% close to other negative values.
is_wave_candidate = logical(is_included .*close_to_ssh.* ...
    close_to_negative);

% Only include if there is a single surface
% is_single_candidate = is_wave_candidate .* has_single_surface;

% Different metrics for being wave-affected

% Just negative, ice, and not too long.
% naive_under = logical((height < 0).*is_wave_candidate);

% Adjusted height is negative
% is_under = logical((height_adjusted < 0).*is_wave_candidate);

% Height is negative beyond ssh or height variance
% is_under_ssh_var = logical((height_adjusted < -wave_cutoff_ssh).*is_wave_candidate);
% is_under_height_var = logical((height_adjusted < -wave_cutoff_height).*is_wave_candidate);
% is_under_both_ex = logical((height_adjusted < -both_cutoff_height).*is_single_candidate);
is_under_both_var = logical((height_adjusted < -both_cutoff_height).*is_wave_candidate);

%% Along-track number of points in each window.

% The number of segments that go into each moving average window.
% This excludes adjusted height points.
AT_N_all = movsum(is_included_all,AT_window,'SamplePoints',dist);
AT_N_conc = movsum(is_included_conc,AT_window,'SamplePoints',dist);
% Need at least a density of 1/50 meters for us to consider including one of the variables computed on the moving average window.

use_AT = 1*(AT_N_all > sum(AT_window) / 30);
use_AT(use_AT == 0) = nan;
use_AT(dist < AT_window(1)) = nan;
use_AT(abs(max(dist) - dist) < AT_window(2)) = nan;

% This is now the number along-track in appropriate windows.
AT_N = use_AT.*AT_N_all;


%% Along-track WAF

% wave_area_frac_naive = 2 * movsum(seg_len .* naive_under,window_25k,'samplepoints',dist) ./ movsum(seg_len .* is_ice,window_25k,'samplepoints',dist);
% wave_area_frac = 2 * movsum(seg_len .* is_under,window_25k,'samplepoints',dist) ./ movsum(seg_len .* is_ice,window_25k,'samplepoints',dist);
% wave_area_frac_ssh = 2 * movsum(seg_len .* is_under_ssh_var,window_25k,'samplepoints',dist) ./ movsum(seg_len .* is_ice,window_25k,'samplepoints',dist);
% wave_area_frac_height = 2 * movsum(seg_len .* is_under_height_var,window_25k,'samplepoints',dist) ./ movsum(seg_len .* is_ice,window_25k,'samplepoints',dist);
wave_area_frac_both = 2 * movsum(seg_len .* is_under_both_var,AT_window,'samplepoints',dist) ./ movsum(seg_len .* is_ice,AT_window,'samplepoints',dist);
% wave_area_frac_both_ex = 2 * movsum(seg_len .* is_under_both_ex,window_25k,'samplepoints',dist) ./ movsum(seg_len .* is_ice,window_25k,'samplepoints',dist);

AT_WAF = use_AT.*wave_area_frac_both;

%% Collecting data at each segment of the IS2 track.
% This can be saved into AT_stats
% Or it will be binned for DS_stats.

% Linear ice fraction. Ignore any along-track values where we don't have
% enough point density.
AT_LIF = use_AT.*movsum(seg_len.*is_included,AT_window,'samplepoints',dist) ./ ...
    movsum(seg_len.*is_included_all,AT_window,'samplepoints',dist);

% Also an "adjusted" LIF where negative heights go to zero.
AT_positive = height_adjusted > 0;
AT_LIF_adj = use_AT.*movsum(seg_len.*is_included.*AT_positive,AT_window,'samplepoints',dist) ...
    ./ movsum(seg_len.*is_included_all,AT_window,'samplepoints',dist);

% LIF computed by setting specular returns to ocean or dark leads to ocean,
% only.
AT_LIF_spec = use_AT.*movsum(seg_len.*is_not_spec.*is_included_all,AT_window,'samplepoints',dist) ./ movsum(seg_len.*is_included_all,AT_window,'samplepoints',dist);
AT_LIF_dark = use_AT.* movsum(seg_len.*is_not_dark.*is_included_all,AT_window,'samplepoints',dist) ./ movsum(seg_len.*is_included_all,AT_window,'samplepoints',dist);

% SIC is segment length weighted mean
% because there are huge nan-like values of conc, we omit them. 
AT_SIC =  use_AT.* (1/100).*movsum(seg_len.*conc.*is_included_conc,AT_window,'samplepoints',dist) ./ movsum(seg_len.*is_included_conc,AT_window,'samplepoints',dist);

if IS2_obj.v6

    % SIC is segment length weighted mean
    AT_SIC_amsr =  use_AT.*(1/100).*movsum(seg_len.*conc_amsr.*is_included_all,AT_window,'samplepoints',dist) ./ movsum(seg_len.*is_included_all,AT_window,'samplepoints',dist);

end

% FSD is segment length ratio. Ignore floes on the boundaries that wouldn't
% have enough points to make them.
AT_RFSD = use_AT(floeind).*movsum(floe_length.^3,AT_window,'samplepoints',dist(floeind)) ./ movsum(floe_length.^2,AT_window,'samplepoints',dist(floeind));

% FSD is mean floe length
AT_MFSD =  use_AT(floeind).*movsum(floe_length,AT_window,'samplepoints',dist(floeind)) ./ movsum(floe_length.^0,AT_window,'samplepoints',dist(floeind));

% Along-track variance. Again ignore those points without enough segments
% in them.
AT_E =  use_AT.*movsum(seg_len.*height_adjusted.^2.*is_included_all,AT_window,'samplepoints',dist) ./ movsum(seg_len.*is_included_all,AT_window,'samplepoints',dist);

%% Along-track averaged height
% Exclude points where there isn't enough segments to make the moving
% average work.

AT_height = use_AT.*movsum(height_adjusted.*seg_len.*is_included,AT_window,'omitnan','samplepoints',dist) ...
    ./ movsum(seg_len.*is_included,AT_window,'omitnan','samplepoints',dist);

% Exclude locations where there's a big deviation for the purposes of
% computing along-track variance.
include_variance = (height_adjusted - AT_height < 1);

% Along-track variance. Ditto for points.
AT_var = use_AT.*movsum((height_adjusted - AT_height).^2 .* seg_len .* is_included.*include_variance,AT_window,'omitnan','samplepoints',dist) ...
    ./ movsum(seg_len.*is_included.*include_variance,AT_window,'omitnan','samplepoints',dist);


%%

if do_fullsample

    AT_stats.timer = timer;

    AT_stats.lat = lat;
    AT_stats.lon = lon;
    AT_stats.use_AT = use_AT;
    AT_stats.LIF = AT_LIF;
    AT_stats.LIF_adj = AT_LIF_adj;
    AT_stats.WAF = AT_WAF;
    AT_stats.SIC = AT_SIC;
    AT_stats.FSD = AT_RFSD;

    AT_stats.N = AT_N;
    AT_stats.N_all = AT_N_all;
    AT_stats.H = AT_height;
    AT_stats.H_var = AT_var;

    if IS2_obj.v6
        AT_stats.SIC_amsr = AT_SIC_amsr;
    end

    % This field is not on the same moving window as the others. Therefore
    % we don't need to remove any adjustments.
    AT_stats.height_adj = height_adjusted;
    % Remove outlier values
    AT_stats.height_adj(AT_stats.height_adj > 10) = nan;

end


%% Downsample to other grid
% Creates downsampled grid files.

if do_downsample

    % Downsampling has to take raw data and bring it to a downsampled grid.
    % It

    dist_ind = floor(dist/AT_resolution/2)+1;

    [downscale_inds,~,ind_mapper] = unique(dist_ind);

    floe_mapper = ind_mapper(floeind);

    DS_stats.timer = timer;

    % Take mean number of segments used to produce an along-track grid.

    % Have to be careful here. AT_N is the number of points used in each
    % moving window. It can be too small and so we exclude it from other
    % computations.

    % The strict number that fall into each bin. Whether for non-ice or for ice.
    % This is for fields that are being collected into each bin.
    DS_stats.N_strict = accumarray(ind_mapper,use_AT,[length(downscale_inds) 1],@nansum);
    DS_stats.N_strict_ice = accumarray(ind_mapper,use_AT.*is_included_all,[length(downscale_inds) 1],@nansum);

    % Lats, lons, distance-based metrics.
    DS_stats.lat = accumarray(ind_mapper,use_AT.*lat,[length(downscale_inds) 1],@nanmean);
    DS_stats.lon = accumarray(ind_mapper,use_AT.*lon,[length(downscale_inds) 1],@nanmean);
    DS_stats.D_to_edge = accumarray(ind_mapper,use_AT.*D_to_edge,[length(downscale_inds) 1],@nanmean);

    % The following are of ice fields and therefore need to include the
    % number of included segments.
    DS_stats.H = accumarray(ind_mapper,use_AT.*height_adjusted,[length(downscale_inds) 1],@nanmean);
    DS_stats.H_var = accumarray(ind_mapper,use_AT.*height_adjusted,[length(downscale_inds) 1],@nanstd);

    % CH-derived statistics
    % Since these include fields derived on moving averages, we need to
    % appropriately weight the fact that each value uses an estimate of
    % a certain number.
    DS_stats.WAF = accumarray(ind_mapper,AT_WAF,[length(downscale_inds) 1],@nanmean);
    DS_stats.LIF = accumarray(ind_mapper,AT_LIF,[length(downscale_inds) 1],@nanmean);
    DS_stats.LIF_spec = accumarray(ind_mapper,AT_LIF_spec,[length(downscale_inds) 1],@nanmean);
    DS_stats.LIF_dark = accumarray(ind_mapper,AT_LIF_dark,[length(downscale_inds) 1],@nanmean);

    DS_stats.SIC = accumarray(ind_mapper,AT_SIC,[length(downscale_inds) 1],@nanmean);

    if IS2_obj.v6

        DS_stats.SIC_amsr = accumarray(ind_mapper,AT_SIC_amsr,[length(downscale_inds) 1],@nanmean);

    end

    % Along-track statistics

    DS_stats.E = accumarray(ind_mapper,AT_E,[length(downscale_inds) 1],@nanmean);

    if ~isempty(floe_mapper)

        DS_stats.Nfloe = accumarray(floe_mapper,use_AT(floeind),[length(downscale_inds) 1],@nansum);
        DS_stats.MFSD = accumarray(floe_mapper,AT_MFSD,[length(downscale_inds) 1],@nanmean);
        DS_stats.RFSD = accumarray(floe_mapper,AT_RFSD,[length(downscale_inds) 1],@nanmean);

    else
        DS_stats.Nfloe = nan*DS_stats.H;
        DS_stats.MFSD = nan*DS_stats.H;
        DS_stats.RFSD = nan*DS_stats.H;
    end

end

