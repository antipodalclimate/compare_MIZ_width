function AT_stats = downscale_AT_statistics(IS2_obj,window_output)

%% Pull out details from the object
dist = IS2_obj.dist;
revdist = max(dist) - dist;

swap_ind = find(revdist < dist,1);

D_to_edge = min(dist,revdist);
D_to_edge(swap_ind:end) = -D_to_edge(swap_ind:end);

height = IS2_obj.height;
is_ice = IS2_obj.is_ice;
is_ocean = IS2_obj.is_ocean;
seg_len = IS2_obj.seg_len;
lat = IS2_obj.lat;
lon = IS2_obj.lon;
timer = IS2_obj.timer;
conc = IS2_obj.conc;

%%
% Identify local moving average
window_1k = 1000;
window_10k = 10000; % meters - size of moving window
window_25k = 25000; % meters - size of moving window

slide_25k = [0 25000];
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

% adjust for deviation from local ssh
height_adjusted = height - ssh_interp;

% Not too long, is ice, and has positive moving average
% naive_reasonable = logical(not_too_long .* is_ice);

% Has positive points nearby
close_to_positive = moving_pos >= 2;
% Has negative points nearby - only for those values that are negative
close_to_negative = moving_neg >= 2;

% Included points have positive points nearby, aren't too long, and are
% identified as ice. Criteria I1-I2.
is_included = logical(not_too_long .* is_ice .* close_to_positive);

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

%% Along-track WAF

% wave_area_frac_naive = 2 * movsum(seg_len .* naive_under,window_25k,'samplepoints',dist) ./ movsum(seg_len .* is_ice,window_25k,'samplepoints',dist);
% wave_area_frac = 2 * movsum(seg_len .* is_under,window_25k,'samplepoints',dist) ./ movsum(seg_len .* is_ice,window_25k,'samplepoints',dist);
% wave_area_frac_ssh = 2 * movsum(seg_len .* is_under_ssh_var,window_25k,'samplepoints',dist) ./ movsum(seg_len .* is_ice,window_25k,'samplepoints',dist);
% wave_area_frac_height = 2 * movsum(seg_len .* is_under_height_var,window_25k,'samplepoints',dist) ./ movsum(seg_len .* is_ice,window_25k,'samplepoints',dist);
wave_area_frac_both = 2 * movsum(seg_len .* is_under_both_var,window_25k,'samplepoints',dist) ./ movsum(seg_len .* is_ice,window_25k,'samplepoints',dist);
% wave_area_frac_both_ex = 2 * movsum(seg_len .* is_under_both_ex,window_25k,'samplepoints',dist) ./ movsum(seg_len .* is_ice,window_25k,'samplepoints',dist);

AT_WAF = wave_area_frac_both;



%% Along-track LIF, SIC, mean floe size
AT_LIF = movsum(seg_len.*is_ice,slide_25k,'samplepoints',dist) ./ movsum(seg_len.*(is_ice + is_ocean),slide_25k,'samplepoints',dist);

% SIC is segment length weighted mean
AT_SIC = (1/100)*movsum(seg_len.*conc,slide_25k,'samplepoints',dist) ./ movsum(seg_len,slide_25k,'samplepoints',dist);

% FSD is segment length ratio
AT_RFSD = movsum(floe_length.^3,slide_25k,'samplepoints',dist(floeind)) ./ movsum(floe_length.^2,slide_25k,'samplepoints',dist(floeind));

% FSD is mean floe length
AT_MFSD = movsum(floe_length,slide_25k,'samplepoints',dist(floeind)) ./ movsum(floe_length.^0,slide_25k,'samplepoints',dist(floeind));

%%
AT_E = movsum(seg_len.*height_adjusted.^2,slide_25k,'samplepoints',dist) ./ movsum(seg_len,slide_25k,'samplepoints',dist);

%% Downsample to smaller grid

dist_ind = floor(dist/window_output/2)+1;

[downscale_inds,~,ind_mapper] = unique(dist_ind);
ninds = length(downscale_inds);

floe_mapper = ind_mapper(floeind);


AT_stats.timer = timer;
AT_stats.N = accumarray(ind_mapper,1,[length(downscale_inds) 1],@sum);
AT_stats.lat = accumarray(ind_mapper,lat,[length(downscale_inds) 1],@sum)./AT_stats.N;
AT_stats.lon = accumarray(ind_mapper,lon,[length(downscale_inds) 1],@sum)./AT_stats.N;
AT_stats.D_to_edge = accumarray(ind_mapper,D_to_edge,[length(downscale_inds) 1],@sum)./AT_stats.N;

% CH-derived statistics
AT_stats.WAF = accumarray(ind_mapper,AT_WAF,[length(downscale_inds) 1],@sum)./AT_stats.N;
AT_stats.LIF = accumarray(ind_mapper,AT_LIF,[length(downscale_inds) 1],@sum)./AT_stats.N;
AT_stats.SIC = accumarray(ind_mapper,AT_SIC,[length(downscale_inds) 1],@sum)./AT_stats.N;

% Along-track statistics
AT_stats.E = accumarray(ind_mapper,AT_E,[length(downscale_inds) 1],@sum)./AT_stats.N;
AT_stats.H = accumarray(ind_mapper,height_adjusted,[length(downscale_inds) 1],@sum)./AT_stats.N;

if ~isempty(floe_mapper)

    AT_stats.Nfloe = accumarray(floe_mapper,1,[length(downscale_inds) 1],@sum);
    AT_stats.MFSD = accumarray(floe_mapper,AT_MFSD,[length(downscale_inds) 1],@sum)./AT_stats.Nfloe;
    AT_stats.RFSD = accumarray(floe_mapper,AT_RFSD,[length(downscale_inds) 1],@sum)./AT_stats.Nfloe;

else
    AT_stats.Nfloe = nan*AT_stats.H;
    AT_stats.MFSD = nan*AT_stats.H;
    AT_stats.RFSD = nan*AT_stats.H;
end


