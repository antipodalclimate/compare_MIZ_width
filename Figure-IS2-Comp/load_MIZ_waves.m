% Load in the segmented statistics. Each is an array of stats which is
% indexed by
% nT - number of tracks
% nB - number of beams (this is frequently double the actual number of
% beams because we count Northward and Southward halves of a single track
% the same way

nT = size(MIZ_DATA.timer,1);
nB = size(MIZ_DATA.timer,2);

wave_thresh = 0.075;

%%

% Find out if a wavy track at some point

[iswavy,monthval,yearval,npoints,isstrong,nameid,beamid,hasMIZ,hasMIZ_amsr,int_id] ...
    = deal(MIZ_DATA.WAF);


%%
num_beams = size(IS2_DATA.is_strong,2);

search_dir = {'first','last'};
[MIZ_width,MIZ_width_amsr,WAF_width,track_is_strong,track_timeval,track_nMIZ] = deal(zeros(nT,nB));

for i = 1:nT
    for j = 1:nB % Both forward and reverse tracks, all beams
        %%

        if ~isempty(MIZ_DATA.lat{i,j})

            beaminds = (nB/2)*floor(j/2)+1:(nB/2)*floor(j/2) + nB/2;

            og_beam = mod(j,num_beams);

            if og_beam == 0
                og_beam = num_beams;
            end

            lat = MIZ_DATA.lat{i,j}; 
            WAF = MIZ_DATA.WAF{i,j};
            D = MIZ_DATA.D_to_MIZ{i,j};
            SIC = MIZ_DATA.SIC{i,j}; 
            N = MIZ_DATA.N_strict{i,j};
            LIF = MIZ_DATA.LIF{i,j};
            H = MIZ_DATA.H{i,j}; 

            dummy = 0*WAF; 

            keepvals = N > cutoff_N & SIC > 0.15 & LIF > 0.15 ... 
                & ~isnan(D) & ~isinf(SIC) & ~isnan(H);

            %%
            AT_usable{i,j} = keepvals; 
            WAF(~keepvals) = nan; 
            SIC(~keepvals) = nan; 
            N(~keepvals) = nan; 
  
            %% Now calculations abiyt the MIZ width,etc

            MIZ_width(i,j) = sum(SIC < 0.8 & D <= 0);
            MIZ_width_dist(i,j) = -min(D);

            track_nMIZ(i,j) = sum(SIC < 0.8);
            
            haswaves = sum(WAF > wave_thresh & N > cutoff_N);
            haswaves_MIZ = sum(WAF > wave_thresh & D <= 0 & N > cutoff_N);

            WAF_width(i,j) = haswaves;

            track_is_strong(i,j) = IS2_DATA.is_strong(i,og_beam);

            track_timeval(i,j) = 0 + month(MIZ_DATA.timer{i,j});
            
            % Number of points with sufficient values of N
            track_npoints(i,j) = sum(N > cutoff_N);

            if IS2_DATA.v6
                              
                first80 = find(MIZ_DATA.SIC_amsr{i,j} > 0.8,1,search_dir{MIZ_DATA.is_reversed(i,j)+1});
                if ~isempty(first80)
                    MIZ_width_amsr(i,j) = sum(MIZ_DATA.SIC_amsr{i,j} < 0.8 & D <= 0);
                    MIZ_width_dist_amsr(i,j) = -(min(D) - D(first80));
                end
            end
            % 
            % if i == 528
            %     disp('oj')
            % end
            % 
            if MIZ_width(i,j) > 0
                hasMIZ{i,j} = MIZ_width(i,j) + 0*dummy; 
            else
                hasMIZ{i,j} = 0 + 0*dummy; 
            end

            if MIZ_width_amsr(i,j) > 0
                hasMIZ_amsr{i,j} = MIZ_width_amsr(i,j) + 0*dummy; 
            else
                hasMIZ_amsr{i,j} = 0 + 0*dummy; 
            end


            if haswaves_MIZ > 0
                iswavy{i,j} = 1 + 0*dummy;
            else
                if haswaves > 0
                    iswavy{i,j} = -1 + 0*dummy;
                else
                    iswavy{i,j} = 0 + 0*dummy;
                end
            end

            isstrong{i,j} = IS2_DATA.is_strong(i,og_beam) + 0*dummy;

            monthval{i,j} = month(MIZ_DATA.timer{i,j}) + 0*dummy;
            yearval{i,j} = year(MIZ_DATA.timer{i,j}) + 0*dummy;
        

            npoints{i,j} = sum(D<=0) + 0*dummy;
            nMIZ{i,j} = track_nMIZ(i,j) + 0*dummy; 
            nameid{i,j} = i + 0*dummy; 
            beamid{i,j} = j + 0*dummy;
            int_id{i,j} = i + (j-1)*nT + 0*dummy; 

        end

    end
end

%%
Nsegvals = vertcat(MIZ_DATA.N_strict{:});
SICvals = vertcat(MIZ_DATA.SIC{:});

if IS2_DATA.v6

    SICvals_amsr = vertcat(MIZ_DATA.SIC_amsr{:});
    biasvals = (SICvals_amsr - SICvals);

end

LIFvals = vertcat(MIZ_DATA.LIF{:});
LIF_spec_vals = vertcat(MIZ_DATA.LIF_spec{:});
LIF_dark_vals = vertcat(MIZ_DATA.LIF_dark{:});

biasvals_LIF = LIFvals - SICvals; 

Hvals = vertcat(MIZ_DATA.H{:});
Evals = vertcat(MIZ_DATA.E{:});
WAFvals = vertcat(MIZ_DATA.WAF{:});
wavytracks = vertcat(iswavy{:});
monthval = vertcat(monthval{:});
yearval = vertcat(yearval{:});
npoints = vertcat(npoints{:});
isstrong = vertcat(isstrong{:});
nameid = vertcat(nameid{:});
beamid = vertcat(beamid{:});
int_id = vertcat(int_id{:}); 
latvals = vertcat(MIZ_DATA.lat{:}); 
lonvals = vertcat(MIZ_DATA.lon{:}); 
keepvals = vertcat(AT_usable{:});

num_MIZ = vertcat(nMIZ{:});
hasMIZ = vertcat(hasMIZ{:});
hasMIZ_amsr = vertcat(hasMIZ_amsr{:});

Dvals = (vertcat(MIZ_DATA.D_to_MIZ{:})/1000);

spacer = 12.5;

Dbins = -1000+spacer/2:spacer:1000;

Bincent = 0.5*(Dbins(1:end-1) + Dbins(2:end));
% Bincent(end+1) = Bincent(end) + Bincent(2) - Bincent(1);

%% Create usable vector

% Problem is that as we eliminate tracks, we also kill of zero-centering.
% Some stencils don't pass, which is fine, but if these stencils are D=0
% then they mess up our counting.  

all_intersections = unique(int_id);

% Don't want things with infinite SIC or nan values.
usable_all = keepvals & ~isnan(nameid); 
used_tracks = IS2_DATA.namearray(unique(nameid(usable_all)));

fprintf('----- \n')
fprintf('Total of %2.0f intersections possible from %2.0f tracks \n',length(all_intersections),length(unique(nameid)));
fprintf('Total of %2.2f million post-processed stencils from %2.0f intersections over %2.0f tracks \n',sum(usable_all)/1e6,length(intersections),length(used_tracks))



%%

% Take the unique values that are NOT usable.
% We want to remove entire intersections that have screwy values. 
% lost_ints = unique(int_id(~usable_all));
% lost_locs = ismember(int_id,lost_ints,'rows');
% usable_all(lost_locs) = 0; 
% used_tracks = IS2_DATA.namearray(unique(nameid(usable_all)));
% intersections = unique(int_id(usable_all));

%% Now include other requirements. 

usable_all = usable_all & npoints > 0;
used_tracks = IS2_DATA.namearray(unique(nameid(usable_all)));
intersections = unique(nameid(usable_all) + (beamid(usable_all)-1)*max(nameid(usable_all)));
fprintf('Has a CIZ: %2.2f million from %2.0f beams across %2.0f tracks \n',sum(usable_all)/1e6,length(intersections),length(used_tracks))

% Need to have values of 
usable_all = usable_all & num_MIZ > 0;
used_tracks = IS2_DATA.namearray(unique(nameid(usable_all)));
intersections = unique(nameid(usable_all) + (beamid(usable_all)-1)*max(nameid(usable_all)));
fprintf('Has a MIZ point: %2.2f million from %2.0f beams across %2.0f tracks \n',sum(usable_all)/1e6,length(intersections),length(used_tracks))

usable_all = usable_all & hasMIZ > 1; 
used_tracks = IS2_DATA.namearray(unique(nameid(usable_all)));
intersections = unique(nameid(usable_all) + (beamid(usable_all)-1)*max(nameid(usable_all)));
fprintf('Has a MIZ before CIZ: %2.2f million from %2.0f beams across %2.0f tracks \n',sum(usable_all)/1e6,length(intersections),length(used_tracks))

% Remove unplottable things. 
usable_all = usable_all & Dvals < max(Dbins) & Dvals > min(Dbins);
% S
