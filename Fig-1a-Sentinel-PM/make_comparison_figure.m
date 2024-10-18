%% Threshold and plot the S1 image
Amp2 = Amp;
Amp_subset = nan(size(Amp2));
Amp_subset(inds) = Amp2(inds); 

% Amp(Amp > 1000) = 1000; 
Amp_norm = Amp2 / nanmax(Amp2(:)); 
Int_norm = Int / nanmax(Int(:));
ampvals = Amp_norm(ID); 
T = graythresh(ampvals);

amp_bin = ampvals; 
amp_bin(ampvals > T) = 1; 
amp_bin(ampvals <= T) = 0; 
amp_smooth = movmean(amp_bin,mov_window,'omitnan','samplepoints',IS2_obj.dist);
