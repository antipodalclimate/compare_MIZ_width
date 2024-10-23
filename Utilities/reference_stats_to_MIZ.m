function reference_stats_to_MIZ(OPTS)

load(OPTS.load_string,'IS2_DATA')

%%

if OPTS.do_weak
    maxB = 6;
else
    maxB = 3;
end

for i = 1:size(IS2_DATA.AT_stats,1)

    beamct = 0;

    for j = 1:6

        if ~isempty((IS2_DATA.AT_stats{i,j}))

            beamct = beamct + 1;

            % Since the satellite passed over the edge twice - two
            % additional measurements

            % Front side
            side_indices{1} = IS2_DATA.AT_stats{i,j}.D_to_edge > 0;
            % Back side
            side_indices{2} = IS2_DATA.AT_stats{i,j}.D_to_edge < 0;

            % can be an error if we are going up/down in lat but lon is
            % swapping across 180 I think. Just need to make sure
            side_swap = [1 sum(diff(sign(IS2_DATA.AT_stats{i,j}.lon))) > 0];
            side_mult = [1 -1];

            search_dir = {'first','last'};


            for side_index = 1:2

                % If we are sure we actually have a side of the pole.
                if side_swap(side_index)

                    indices = side_indices{side_index};

                    offset = maxB * (side_index - 1);

                    % Geographic info
                    MIZ_DATA.timer{i,offset+beamct} = IS2_DATA.AT_stats{i,j}.timer;
                    MIZ_DATA.lat{i,offset+beamct} = IS2_DATA.AT_stats{i,j}.lat(indices);
                    MIZ_DATA.lon{i,offset+beamct} = IS2_DATA.AT_stats{i,j}.lon(indices);
                    
                    MIZ_DATA.D_to_edge{i,offset+beamct} = IS2_DATA.AT_stats{i,j}.D_to_edge(indices);

                    % Number of segments
                    MIZ_DATA.N{i,offset+beamct} = IS2_DATA.AT_stats{i,j}.N(indices);

                    % Mean height of record
                    MIZ_DATA.H{i,offset+beamct} =  IS2_DATA.AT_stats{i,j}.H(indices);

                    % Along-track variance
                    MIZ_DATA.E{i,offset+beamct} =  IS2_DATA.AT_stats{i,j}.E(indices);

                    % Horvat-derived quantities
                    MIZ_DATA.LIF{i,offset+beamct} = IS2_DATA.AT_stats{i,j}.LIF(indices);
                    MIZ_DATA.LIF_spec{i,offset+beamct} = IS2_DATA.AT_stats{i,j}.LIF_spec(indices);
                    MIZ_DATA.LIF_dark{i,offset+beamct} = IS2_DATA.AT_stats{i,j}.LIF_dark(indices);
                    MIZ_DATA.WAF{i,offset+beamct} = IS2_DATA.AT_stats{i,j}.WAF(indices);

                    % Estimate of floe size distribution
                    
                    try
                        MIZ_DATA.RFSD{i,offset+beamct} = IS2_DATA.AT_stats{i,j}.RFSD(indices);
                        MIZ_DATA.MFSD{i,offset+beamct} = IS2_DATA.AT_stats{i,j}.MFSD(indices);
                    catch 
                        
                    end

                    % PM Sea ice concentration
                    MIZ_DATA.SIC{i,offset+beamct} =  IS2_DATA.AT_stats{i,j}.SIC(indices);

                    % PM Sea ice concentration
                    MIZ_DATA.SIC{i,offset+beamct} =  IS2_DATA.AT_stats{i,j}.SIC(indices);

                    if IS2_DATA.v6 == 1

                    % PM Sea ice concentration
                    MIZ_DATA.SIC_amsr{i,offset+beamct} =  IS2_DATA.AT_stats{i,j}.SIC_amsr(indices);

                    end

                    % Look for the location of the MIZ edge. Search from
                    % the front if lat is declining, search from the back
                    % if lat is increasing
                    MIZ_edge = find(MIZ_DATA.SIC{i,offset+beamct} > 0.8,1,search_dir{side_index});

                    % If we can find the MIZ edge, then compute distance to
                    % it
                    if ~isempty(MIZ_edge)
                        MIZ_DATA.D_to_MIZ{i,offset+beamct} = side_mult(side_index)*(MIZ_DATA.D_to_edge{i,offset+beamct} - MIZ_DATA.D_to_edge{i,offset+beamct}(MIZ_edge));
                    else
                        MIZ_DATA.D_to_MIZ{i,offset+beamct} = nan*MIZ_DATA.D_to_edge{i,offset+beamct};
                    end

                    % Have we reversed
                    MIZ_DATA.is_reversed(i,offset+beamct) = side_index - 1;

                end

            end

        end

    end

end

save(OPTS.load_string,'MIZ_DATA','-append');


