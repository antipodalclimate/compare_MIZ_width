
if ~exist('AT_stat_filename')

    load('AT_stats_Brouwer.mat');

else

    load(AT_stat_filename);

end

%%

if do_strong
    maxB = 3;
else
    maxB = 6;
end

for i = 1:size(AT_stats,1)

    beamcount = 0;

    for j = 1:6

        if ~isempty((AT_stats{i,j}))


            % Since the satellite passed over the edge twice - two
            % additional measurements
            front_side = AT_stats{i,j}.D_to_edge > 0;
            reverse_side = AT_stats{i,j}.D_to_edge < 0;

            side_swap = sum(diff(sign(AT_stats{i,j}.lon))) > 0;


            if sum(front_side) > 0

                beamcount = beamcount + 1;

                % Geographic info
                timer{i,beamcount} = AT_stats{i,j}.timer;
                lat{i,beamcount} = AT_stats{i,j}.lat(front_side);
                lon{i,beamcount,3} = AT_stats{i,j}.lon(front_side);

                % Distance from the beginning of the track
                D_to_edge{i,beamcount} = AT_stats{i,j}.D_to_edge(front_side);


                % Rirect fields from height

                % Number of segments
                N{i,beamcount} = AT_stats{i,j}.N(front_side);

                % Mean height of record
                H{i,beamcount} =  AT_stats{i,j}.H(front_side);

                % Along-track variance
                E{i,beamcount} =  AT_stats{i,j}.E(front_side);

                % Horvat-derived quantities
                LIF{i,beamcount} = AT_stats{i,j}.LIF(front_side);
                WAF{i,beamcount} = AT_stats{i,j}.WAF(front_side);

                % Estimate of floe size distribution
                RFSD{i,beamcount} = AT_stats{i,j}.RFSD(front_side);
                MFSD{i,beamcount} = AT_stats{i,j}.MFSD(front_side);

                % PM Sea ice concentration
                SIC{i,beamcount} =  AT_stats{i,j}.SIC(front_side);

                MIZ_edge = find(SIC{i,beamcount} > 0.8,1);

                if ~isempty(MIZ_edge)

                    D_to_MIZ{i,beamcount} = D_to_edge{i,beamcount} - D_to_edge{i,beamcount}(MIZ_edge);

                else
                    D_to_MIZ{i,beamcount} = nan*D_to_edge{i,beamcount};
                end
                

                reversed(i,beamcount) = 0;


            end

            if sum(reverse_side > 0) && side_swap

                %%
                % Geographic info
                timer{i,maxB+beamcount} = AT_stats{i,j}.timer;
                lat{i,maxB+beamcount} = AT_stats{i,j}.lat(reverse_side);
                lon{i,maxB+beamcount,3} = AT_stats{i,j}.lon(reverse_side);

                % Distance from the reverse end of the track
                % Keep negative for future analysis
                D_to_edge{i,maxB+beamcount} = (AT_stats{i,j}.D_to_edge(reverse_side));

                % Rirect fields from height

                % Number of segments
                N{i,maxB+beamcount} = AT_stats{i,j}.N(reverse_side);

                % Mean height of record
                H{i,maxB+beamcount} =  AT_stats{i,j}.H(reverse_side);

                % Along-track variance
                E{i,maxB+beamcount} =  AT_stats{i,j}.E(reverse_side);

                % Horvat-derived quantities
                LIF{i,maxB+beamcount} = AT_stats{i,j}.LIF(reverse_side);
                WAF{i,maxB+beamcount} = AT_stats{i,j}.WAF(reverse_side);

                % Estimate of floe size distribution
                RFSD{i,maxB+beamcount} = AT_stats{i,j}.RFSD(reverse_side);
                MFSD{i,maxB+beamcount} = AT_stats{i,j}.MFSD(reverse_side);

                % PM Sea ice concentration
                SIC{i,maxB+beamcount} =  AT_stats{i,j}.SIC(reverse_side);

                MIZ_edge = find(SIC{i,maxB+beamcount} > 0.8,1,'last');

                if ~isempty(MIZ_edge)

                    D_to_MIZ{i,maxB+beamcount} = -(D_to_edge{i,maxB+beamcount} - D_to_edge{i,maxB+beamcount}(MIZ_edge));

                else
                    D_to_MIZ{i,maxB+beamcount} = nan*D_to_edge{i,maxB+beamcount};
                end

                reversed(i,maxB+beamcount) = 1;

            end



        end

    end
end
