
nusable = sum(usable_all);
nusable_segmented = sum(usable); 


if IS2_DATA.v6

usable = usable & ~isnan(SICvals_amsr); 

SICvals_amsr = SICvals_amsr(usable); 
biasvals = biasvals(usable);

end

Nsegvals = Nsegvals(usable); 
SICvals = SICvals(usable); 
LIFvals = LIFvals(usable); 
LIF_spec_vals = LIF_spec_vals(usable); 
LIF_dark_vals = LIF_dark_vals(usable); 
Dvals = Dvals(usable); 
Hvals = Hvals(usable); 
Evals = Evals(usable); 
WAFvals = WAFvals(usable);
wavytracks = wavytracks(usable); 
biasvals_LIF = biasvals_LIF(usable); 

[binct,~,binval] = histcounts(Dvals,Dbins);

disp('------')
fprintf('Fraction of tracks that are segmented here is %2.2f \n',100*nusable_segmented./nusable);
fprintf('Fraction of the tracks that is MIZ is %2.2f \n',100*sum(Dvals < 0)./numel(Dvals));
fprintf('Fraction of the tracks that is MIZ or wavy is %2.2f \n',100*sum(Dvals < 0 | WAFvals > wave_thresh)./numel(Dvals));
fprintf('Fraction of those that have bias is %2.2f \n',100*sum(WAFvals > wave_thresh)./numel(WAFvals));


% binval(binval == 0) = length(Dbins); 

upval = @(x) prctile(x,75); 
dnval = @(x) prctile(x,25); 

Nvec = accumarray(binval,1,[length(Dbins)-1 1],@sum); 
%%

Nsegvec = accumarray(binval,Nsegvals,[length(Dbins)-1 1],@mean); 
Nsegvar = accumarray(binval,Nsegvals,[length(Dbins)-1 1],@std); 

SICvec = accumarray(binval,SICvals,[length(Dbins)-1 1],@median); 

SICup = accumarray(binval,SICvals,[length(Dbins) - 1 1],upval); 
SICup(isinf(SICup)) = 1; 
SICdn = accumarray(binval,SICvals,[length(Dbins) - 1 1],dnval); 

if IS2_DATA.v6

SICvec_amsr = accumarray(binval,SICvals_amsr,[length(Dbins) - 1 1],@median); 
SICup_amsr= accumarray(binval,SICvals_amsr,[length(Dbins) - 1 1],upval); 
SICdn_amsr = accumarray(binval,SICvals_amsr,[length(Dbins) - 1 1],dnval); 

biasvec = accumarray(binval,biasvals,[length(Dbins) - 1 1],@median); 
biasup= accumarray(binval,biasvals,[length(Dbins) - 1 1],upval); 
biasdn = accumarray(binval,biasvals,[length(Dbins) - 1 1],dnval); 


end

bias_LIFvec = accumarray(binval,biasvals_LIF,[length(Dbins) - 1 1],@median); 
bias_LIFup= accumarray(binval,biasvals_LIF,[length(Dbins) - 1 1],upval); 
bias_LIFdn = accumarray(binval,biasvals_LIF,[length(Dbins) - 1 1],dnval); 


WAFvec = accumarray(binval,WAFvals,[length(Dbins) - 1 1],@median); 
WAFup = accumarray(binval,WAFvals,[length(Dbins) - 1 1],upval); 
WAFdn = accumarray(binval,WAFvals,[length(Dbins) - 1 1],dnval); 

Hvec = accumarray(binval,Hvals,[length(Dbins) - 1 1],@median); 
Hup = accumarray(binval,Hvals,[length(Dbins) - 1 1],upval); 
Hdn = accumarray(binval,Hvals,[length(Dbins) - 1 1],dnval); 


LIFvec = accumarray(binval,LIFvals,[length(Dbins) - 1 1],@median); 
LIFup = accumarray(binval,LIFvals,[length(Dbins) - 1 1],upval); 
LIFdn = accumarray(binval,LIFvals,[length(Dbins) - 1 1],dnval); 

LIF_spec_vec = accumarray(binval,LIF_spec_vals,[length(Dbins) - 1 1],@median); 
LIF_spec_up = accumarray(binval,LIF_spec_vals,[length(Dbins) - 1 1],upval); 
LIF_spec_dn = accumarray(binval,LIF_spec_vals,[length(Dbins) - 1 1],dnval); 

LIF_dark_vec = accumarray(binval,LIF_dark_vals,[length(Dbins) - 1 1],@median); 
LIF_dark_up = accumarray(binval,LIF_dark_vals,[length(Dbins) - 1 1],upval); 
LIF_dark_dn = accumarray(binval,LIF_dark_vals,[length(Dbins) - 1 1],dnval); 



Hvar = accumarray(binval,Hvals,[length(Dbins) - 1 1],@std); 

Evec = accumarray(binval,Evals,[length(Dbins) - 1 1],@mean); 
Evar = accumarray(binval,Evals,[length(Dbins) - 1 1],@std); 

%
xfirst = find(Nvec > 1,1,'first');

dum = find(Nvec > 1,2,'last');

xlast = dum(1); 
