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

for i = 1:nT
    for j = 1:nB % Both forward and reverse tracks, all strong beams
        
        beaminds = (nB/2)*floor(j/2)+1:(nB/2)*floor(j/2) + nB/2; 
        
        WAF = MIZ_DATA.WAF{i,j}; 
        D = MIZ_DATA.D_to_MIZ{i,j};
        N = MIZ_DATA.N{i,j}; 

        haswaves = sum(WAF > wave_thresh & N > 10);
        haswaves_MIZ = sum(WAF > wave_thresh & D < 0 & N > 10);

  
        
        if haswaves_MIZ > 0
            iswavy{i,j} = 1 + 0*WAF; 
        else 
           if haswaves > 0
               iswavy{i,j} = -1 + 0*WAF; 
           else
               iswavy{i,j} = 0 + 0*WAF; 
           end
        end
         
        timeval{i,j} = month(MIZ_DATA.timer{i,j}) + 0*WAF;
        npoints{i,j} = sum(D<=0) + 0*WAF; 


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
Nvals = vertcat(MIZ_DATA.N{:}); 

if IS2_DATA.v6

    SICvals_amsr = vertcat(MIZ_DATA.SIC_amsr{:}); 

end

SICvals = vertcat(MIZ_DATA.SIC{:}); 
LIFvals = vertcat(MIZ_DATA.LIF{:}); 
Hvals = vertcat(MIZ_DATA.H{:}); 
Evals = vertcat(MIZ_DATA.E{:}); 
WAFvals = vertcat(MIZ_DATA.WAF{:}); 
wavytracks = vertcat(iswavy{:}); 
timeval = vertcat(timeval{:}); 
npoints = vertcat(npoints{:}); 

Dvals = (vertcat(MIZ_DATA.D_to_MIZ{:})/1000); 

Dbins = -1000+12.5:25:1000; 
Bincent = 0.5*(Dbins(1:end-1) + Dbins(2:end));
Bincent(end+1) = Bincent(end) + Bincent(2) - Bincent(1); 
