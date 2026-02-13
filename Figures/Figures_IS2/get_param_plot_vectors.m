
nbins = 100*floor(round(sqrt(sum(param_use)))/250) + 1;
nbins = max(nbins,50);

nbins_H = nbins+1; 
nbins_C = min(nbins,100)+1;
nbins_W = nbins; 

Hbins = linspace(-4,8,nbins_H);
Bincent_H = 0.5*(Hbins(1:end-1) + Hbins(2:end));

Wbins = linspace(0,2,nbins_W);
Wbins = Wbins - 0.5*Wbins(2);
Bincent_W = 0.5*(Wbins(1:end-1) + Wbins(2:end));

Cbins = linspace(0,1,nbins_C);
dC = Cbins(2) - Cbins(1); 
Cbins = Cbins + dC/2; 
Bincent_C = 0.5*(Cbins(1:end-1) + Cbins(2:end));


[nH,~,mapper_H] = histcounts(Hvals(param_use),Hbins);
[nW,~,mapper_W] = histcounts(WAFvals(param_use),Wbins);
[nC,~,mapper_C] = histcounts(SICvals_CDR(param_use),Cbins);

mapper_C(SICvals_CDR(param_use) > 1) = length(Cbins)-1;

ncutoff = 10; 

% xlimmer_C = [Cbins(find(nC>ncutoff,1)) Cbins(find(nC>ncutoff,1,'last'))];
% xlimmer_H = [Hbins(find(nH>ncutoff,1)) Hbins(find(nH>ncutoff,1,'last'))];
% xlimmer_W = [Wbins(find(nW>ncutoff,1)) Wbins(find(nW>ncutoff,1,'last'))];
% 

%%

LIF_by_C = accumarray(mapper_C,LIFvals(param_use),[length(Cbins)-1 1],@nanmedian);
LIF_by_C_up = accumarray(mapper_C,LIFvals(param_use),[length(Cbins)-1 1],upval);
LIF_by_C_dn = accumarray(mapper_C,LIFvals(param_use),[length(Cbins)-1 1],dnval);


CDR_by_C = accumarray(mapper_C,SICvals_CDR(param_use),[length(Cbins)-1 1],@nanmedian);
CDR_by_C_up = accumarray(mapper_C,SICvals_CDR(param_use),[length(Cbins)-1 1],upval);
CDR_by_C_dn = accumarray(mapper_C,SICvals_CDR(param_use),[length(Cbins)-1 1],dnval);

AMSR_by_C = accumarray(mapper_C,SICvals_AMSR(param_use),[length(Cbins)-1 1],@nanmedian);
AMSR_by_C_up = accumarray(mapper_C,SICvals_AMSR(param_use),[length(Cbins)-1 1],upval);
AMSR_by_C_dn = accumarray(mapper_C,SICvals_AMSR(param_use),[length(Cbins)-1 1],dnval);


bias_AMSR_by_C = accumarray(mapper_C,biasvals_AMSR(param_use),[length(Cbins)-1 1],@nanmedian);
bias_AMSR_by_C_up = accumarray(mapper_C,biasvals_AMSR(param_use),[length(Cbins)-1 1],upval);
bias_AMSR_by_C_dn = accumarray(mapper_C,biasvals_AMSR(param_use),[length(Cbins)-1 1],dnval);

bias_LIF_by_C = accumarray(mapper_C,biasvals_LIF(param_use),[length(Cbins)-1 1],@nanmedian);
bias_LIF_by_C_up = accumarray(mapper_C,biasvals_LIF(param_use),[length(Cbins)-1 1],upval);
bias_LIF_by_C_dn = accumarray(mapper_C,biasvals_LIF(param_use),[length(Cbins)-1 1],dnval);



%%

LIF_by_H = accumarray(mapper_H,LIFvals(param_use),[length(Hbins)-1 1],@nanmedian);
LIF_by_H_up = accumarray(mapper_H,LIFvals(param_use),[length(Hbins)-1 1],upval);
LIF_by_H_dn = accumarray(mapper_H,LIFvals(param_use),[length(Hbins)-1 1],dnval);


CDR_by_H = accumarray(mapper_H,SICvals_CDR(param_use),[length(Hbins)-1 1],@nanmedian);
CDR_by_H_up = accumarray(mapper_H,SICvals_CDR(param_use),[length(Hbins)-1 1],upval);
CDR_by_H_dn = accumarray(mapper_H,SICvals_CDR(param_use),[length(Hbins)-1 1],dnval);
CDR_std_by_H = accumarray(mapper_H,SICvals_CDR(param_use),[length(Hbins)-1 1],@nanstd);

AMSR_by_H = accumarray(mapper_H,SICvals_AMSR(param_use),[length(Hbins)-1 1],@nanmedian);
AMSR_by_H_up = accumarray(mapper_H,SICvals_AMSR(param_use),[length(Hbins)-1 1],upval);
AMSR_by_H_dn = accumarray(mapper_H,SICvals_AMSR(param_use),[length(Hbins)-1 1],dnval);

bias_AMSR_by_H = accumarray(mapper_H,biasvals_AMSR(param_use),[length(Hbins)-1 1],@nanmedian);
bias_AMSR_by_H_up = accumarray(mapper_H,biasvals_AMSR(param_use),[length(Hbins)-1 1],upval);
bias_AMSR_by_H_dn = accumarray(mapper_H,biasvals_AMSR(param_use),[length(Hbins)-1 1],dnval);

bias_LIF_by_H = accumarray(mapper_H,biasvals_LIF(param_use),[length(Hbins)-1 1],@nanmedian);
bias_LIF_by_H_up = accumarray(mapper_H,biasvals_LIF(param_use),[length(Hbins)-1 1],upval);
bias_LIF_by_H_dn = accumarray(mapper_H,biasvals_LIF(param_use),[length(Hbins)-1 1],dnval);

%%
LIF_by_W = accumarray(mapper_W,LIFvals(param_use),[length(Wbins)-1 1],@nanmedian);
LIF_by_W_up = accumarray(mapper_W,LIFvals(param_use),[length(Wbins)-1 1],upval);
LIF_by_W_dn = accumarray(mapper_W,LIFvals(param_use),[length(Wbins)-1 1],dnval);


CDR_by_W = accumarray(mapper_W,SICvals_CDR(param_use),[length(Wbins)-1 1],@nanmedian);
CDR_by_W_up = accumarray(mapper_W,SICvals_CDR(param_use),[length(Wbins)-1 1],upval);
CDR_by_W_dn = accumarray(mapper_W,SICvals_CDR(param_use),[length(Wbins)-1 1],dnval);
CDR_std_by_W = accumarray(mapper_W,SICvals_CDR(param_use),[length(Wbins)-1 1],@nanstd);


AMSR_by_W = accumarray(mapper_W,SICvals_AMSR(param_use),[length(Wbins)-1 1],@nanmedian);
AMSR_by_W_up = accumarray(mapper_W,SICvals_AMSR(param_use),[length(Wbins)-1 1],upval);
AMSR_by_W_dn = accumarray(mapper_W,SICvals_AMSR(param_use),[length(Wbins)-1 1],dnval);

bias_AMSR_by_W = accumarray(mapper_W,biasvals_AMSR(param_use),[length(Wbins)-1 1],@nanmedian);
bias_AMSR_by_W_up = accumarray(mapper_W,biasvals_AMSR(param_use),[length(Wbins)-1 1],upval);
bias_AMSR_by_W_dn = accumarray(mapper_W,biasvals_AMSR(param_use),[length(Wbins)-1 1],dnval);

bias_LIF_by_W = accumarray(mapper_W,biasvals_LIF(param_use),[length(Wbins)-1 1],@nanmedian);
bias_LIF_by_W_up = accumarray(mapper_W,biasvals_LIF(param_use),[length(Wbins)-1 1],upval);
bias_LIF_by_W_dn = accumarray(mapper_W,biasvals_LIF(param_use),[length(Wbins)-1 1],dnval);

%% Now do 2D bias maps

% Here by W and H

nseg_sub = Nsegvals(param_use); 
bias_sub_LIF = biasvals_LIF(param_use); 
bias_sub_AMSR = biasvals_AMSR(param_use); 

[nWH,~,~,iW,iH] = histcounts2(WAFvals(param_use), Hvals(param_use), Wbins, Hbins);

good_indices = (iW>=1 & iW<=nbins_W) & (iH>=1 & iH<=nbins_H);

% Make W the *column* index so W plots on x-axis
lindex = sub2ind([nbins_H-1, nbins_W-1], iH(good_indices), iW(good_indices));   % rows=H, cols=W

% counts array in same orientation: rows=H, cols=W
nWH = nWH.';   % histcounts2 returns nWH as nW-by-nH (W rows, H cols)

bias_LIF_by_WH = accumarray(lindex, bias_sub_LIF(good_indices),[(nbins_H-1)*(nbins_W-1) 1], @nanmedian);
bias_LIF_by_WH = reshape(bias_LIF_by_WH, [nbins_H-1, nbins_W-1]);  % rows=H, cols=W
% apply count cutoff
bias_LIF_by_WH(nWH < 2*sqrt(ncutoff)) = NaN;

bias_AMSR_by_WH = accumarray(lindex, bias_sub_AMSR(good_indices),[(nbins_H-1)*(nbins_W-1) 1], @nanmedian);
bias_AMSR_by_WH = reshape(bias_AMSR_by_WH, [nbins_H-1, nbins_W-1]);  % rows=H, cols=W
% apply count cutoff
bias_AMSR_by_WH(nWH < 2*sqrt(ncutoff)) = NaN;

%% Here by W and C

[nWC,~,~,iW,iC] = histcounts2(WAFvals(param_use), SICvals_CDR(param_use), Wbins, Cbins);

% Make W the *column* index so W plots on x-axis

good_indices = (iW>=1 & iW<=nbins_W) & (iC>=1 & iC<=nbins_C);

lindex = sub2ind([nbins_C-1, nbins_W-1], iC(good_indices), iW(good_indices));   % rows=H, cols=W

% counts array in same orientation: rows=H, cols=W
nWC = nWC.';   % histcounts2 returns nWC as nW-by-nH (W rows, H cols)

bias_LIF_by_WC = accumarray(lindex, bias_sub_LIF(good_indices),[(nbins_W-1)*(nbins_C-1) 1], @nanmedian);
bias_LIF_by_WC = reshape(bias_LIF_by_WC, [nbins_C-1, nbins_W-1]);  % rows=H, cols=W
% apply count cutoff
bias_LIF_by_WC(nWC < sqrt(ncutoff)) = NaN;

bias_AMSR_by_WC = accumarray(lindex, bias_sub_AMSR(good_indices),[(nbins_W-1)*(nbins_C-1) 1], @nanmedian);
bias_AMSR_by_WC = reshape(bias_AMSR_by_WC, [nbins_C-1, nbins_W-1]);  % rows=H, cols=W
% apply count cutoff
bias_AMSR_by_WC(nWC < sqrt(ncutoff)) = NaN;

%% Now H and C

[nHC,~,~,iH,iC] = histcounts2(Hvals(param_use), SICvals_CDR(param_use), Hbins, Cbins);
good_indices = (iH>=1 & iH<=nbins_H) & (iC>=1 & iC<=nbins_C);

% Make W the *column* index so W plots on x-axis
lindex = sub2ind([nbins_C-1, nbins_H-1], iC(good_indices), iH(good_indices));   % rows=H, cols=W

% counts array in same orientation: rows=H, cols=W
nHC = nHC.';   % histcounts2 returns nWC as nW-by-nH (W rows, H cols)

bias_LIF_by_HC = accumarray(lindex, bias_sub_LIF(good_indices),[(nbins_H-1)*(nbins_C-1) 1], @nanmedian);
bias_LIF_by_HC = reshape(bias_LIF_by_HC, [nbins_C-1, nbins_H-1]);  % rows=H, cols=W
% apply count cutoff
bias_LIF_by_HC(nHC < sqrt(ncutoff)) = NaN;

bias_AMSR_by_HC = accumarray(lindex, bias_sub_AMSR(good_indices),[(nbins_H-1)*(nbins_C-1) 1], @nanmedian);
bias_AMSR_by_HC = reshape(bias_AMSR_by_HC, [nbins_C-1, nbins_H-1]);  % rows=H, cols=W
% apply count cutoff
bias_AMSR_by_HC(nHC < sqrt(ncutoff)) = NaN;

nseg_by_HC = accumarray(lindex, nseg_sub(good_indices),[(nbins_H-1)*(nbins_C-1) 1], @nanmedian);
nseg_by_HC = reshape(nseg_by_HC, [nbins_C-1, nbins_H-1]);  % rows=H, cols=W
nseg_by_HC(nHC < sqrt(ncutoff)) = NaN;


%% Let's pull out the positive biased values

overbias = bias_LIF_by_HC < -0.1;

is_overbiased_sub = false(size(iH));       
is_overbiased_sub(good_indices) = overbias( sub2ind([nbins_C-1, nbins_H-1], ...
                                         iC(good_indices), iH(good_indices)) );

% If you want a NEW param_use in the *original full* vector length:
idx = find(param_use);                         % original indices used in the subset
param_use_overbiased = false(size(param_use)); % same size as original logical
param_use_overbiased(idx(is_overbiased_sub)) = true;

