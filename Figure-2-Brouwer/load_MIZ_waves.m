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
iswavy = MIZ_DATA.WAF;
timeval = MIZ_DATA.WAF;
npoints = MIZ_DATA.WAF;
isstrong = MIZ_DATA.WAF;
nameid = MIZ_DATA.WAF; 

num_beams = size(IS2_DATA.is_strong,2);


[MIZ_width,WAF_width,track_is_strong,track_timeval,track_nMIZ] = deal(zeros(nT,nB));

for i = 1:nT
    for j = 1:nB % Both forward and reverse tracks, all beams
        %%

        if ~isempty(MIZ_DATA.lat{i,j})

            beaminds = (nB/2)*floor(j/2)+1:(nB/2)*floor(j/2) + nB/2;

            og_beam = mod(j,num_beams);

            if og_beam == 0
                og_beam = num_beams;
            end

            WAF = MIZ_DATA.WAF{i,j};
            D = MIZ_DATA.D_to_MIZ{i,j};
            N = MIZ_DATA.Nseg{i,j};

            MIZ_width(i,j) = sum(MIZ_DATA.SIC{i,j} < 0.8 & MIZ_DATA.D_to_MIZ{i,j} <= 0);

            haswaves = sum(WAF > wave_thresh & N > 10);
            haswaves_MIZ = sum(WAF > wave_thresh & D <= 0 & N > 10);

            WAF_width(i,j) = haswaves;

            track_is_strong(i,j) = IS2_DATA.is_strong(i,og_beam);

            track_timeval(i,j) = 0 + month(MIZ_DATA.timer{i,j});

            track_nMIZ(i,j) = sum(D<=0);

            track_npoints(i,j) = sum(N > 10);

            if IS2_DATA.v6
                              
                first80 = find(MIZ_DATA.SIC_amsr{i,j} > 0.8,1);
                if ~isempty(first80)
                    MIZ_width_amsr(i,j) = sum(MIZ_DATA.SIC_amsr{i,j}(1:first80) < 0.8);
                end
            end


            if haswaves_MIZ > 1
                iswavy{i,j} = 1 + 0*WAF;
            else
                if haswaves > 0
                    iswavy{i,j} = -1 + 0*WAF;
                else
                    iswavy{i,j} = 0 + 0*WAF;
                end
            end

            isstrong{i,j} = IS2_DATA.is_strong(i,og_beam) + 0*WAF;

            timeval{i,j} = month(MIZ_DATA.timer{i,j}) + 0*WAF;
            npoints{i,j} = sum(D<=0) + 0*WAF;

            nameid{i,j} = i + 0*WAF; 

        end
        %     if ~isempty(MIZ_DATA.SIC{i,j})
        %
        %         % Distance from the edge
        %         var1 = vertcat(MIZ_DATA.D_to_edge{i,beaminds});
        %         var2 = vertcat(MIZ_DATA.SIC{i,beaminds});
        %
        %         Nvals = vertcat(MIZ_DATA.N{i,beaminds});
        %
        %
        %         % Usable values to do the correlation
        %         usable = var1 < 2e5 & Nvals > 100;
        %
        %         if sum(usable) > 1
        %
        %         cmat = corrcoef(var1(usable),var2(usable));
        %
        %         CC(i,j) = cmat(1,2);
        %
        %         end
        %
        %     else
        %         CC(i,j) = nan;
        %     end
        %
    end
end

%%
Nsegvals = vertcat(MIZ_DATA.Nseg{:});



SICvals = vertcat(MIZ_DATA.SIC{:});


if IS2_DATA.v6

    SICvals_amsr = vertcat(MIZ_DATA.SIC_amsr{:});
    biasvals = abs(SICvals_amsr - SICvals);

end


% namevals = vertcat(MIZ_DATA.names(:));

LIFvals = vertcat(MIZ_DATA.LIF{:});
LIF_spec_vals = vertcat(MIZ_DATA.LIF_spec{:});
LIF_dark_vals = vertcat(MIZ_DATA.LIF_dark{:});


Hvals = vertcat(MIZ_DATA.H{:});
Evals = vertcat(MIZ_DATA.E{:});
WAFvals = vertcat(MIZ_DATA.WAF{:});
wavytracks = vertcat(iswavy{:});
timeval = vertcat(timeval{:});
npoints = vertcat(npoints{:});
isstrong = vertcat(isstrong{:});
nameid = vertcat(nameid{:});

Dvals = (vertcat(MIZ_DATA.D_to_MIZ{:})/1000);

Dbins = -1000+12.5:25:1000;
Bincent = 0.5*(Dbins(1:end-1) + Dbins(2:end));
% Bincent(end+1) = Bincent(end) + Bincent(2) - Bincent(1);

%% Create usable vector

% Don't want things outside of our bins, or infinite SIC or nan values.
usable_all = ~isnan(Dvals) &~isinf(SICvals) & ~isnan(Hvals) & Dvals < max(Dbins) & Dvals > min(Dbins);