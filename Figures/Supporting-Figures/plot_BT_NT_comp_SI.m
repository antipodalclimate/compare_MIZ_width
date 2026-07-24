
comp_usable = (SSMI_BS_SH > 0.15 & AMSR_BS_SH > 0.15) & CDR_daily_SH < 0.8; 


is_NT = comp_usable & (CDR_daily_SH == SSMI_NT_SH) & (SSMI_NT_SH ~= SSMI_BS_SH); 
is_BS = comp_usable & (CDR_daily_SH == SSMI_BS_SH) & (SSMI_NT_SH ~= SSMI_BS_SH); 


% 
n_comp = sum(comp_usable & ~isnan(CDR_daily_SH) & (SSMI_NT_SH ~= SSMI_BS_SH),'all','omitmissing');
n_NT = sum(is_NT ,'all','omitmissing');
n_BS = sum(is_BS,'all','omitmissing');

fprintf('Total points is %2.2f million, with %2.2f from NT (%2.0f%%) and %2.2f (%2.0f%%)) from BS \n',n_comp/1e6,n_NT/1e6,100*n_NT/n_comp,n_BS/1e6,100*n_BS/n_comp);
fprintf('Median BT-derived SIC is %2.1f, compared to %2.1f from NT \n',100*mean(CDR_daily_SH(is_BS)),100*mean(CDR_daily_SH(is_NT)))


%% 


del_BS = SSMI_BS_SH - AMSR_BS_SH;
del_BS(~comp_usable) = nan; 

del_NT2 = AMSR_BS_SH - AMSR_NT_SH; 
del_NT2(~comp_usable) = nan; 

low_del = abs(del_BS) <= 0.15; 
high_del = abs(del_BS) > 0.15; 

low_del_NT2 = abs(del_NT2) <= 0.15; 
high_del_NT2  = abs(del_NT2) > 0.15; 


fprintf('%% of points with BT-sensor difference above 15%% is %2.2f \n',100*sum(high_del,'all','omitmissing')./(sum(low_del + high_del,'all','omitmissing'))); 
fprintf('The median CDR for these points is %2.2f \n',100*mean(CDR_daily_SH(high_del),'omitnan'))

fprintf('%% of points with NT2-BT difference above 15%% is %2.2f \n',100*sum(high_del_NT2,'all','omitmissing')./(sum(low_del_NT2 + high_del_NT2,'all','omitmissing'))); 
fprintf('The median CDR for these points is %2.2f \n',100*mean(CDR_daily_SH(high_del_NT2),'omitnan'))


%%

% figure(1)

Ax{1} = subplot('position',[.05 .05 .5 .9]);

worldmap([-90 -55],[-180 180]);

pcolorm(lat_SH,lon_SH,sum(is_BS,3)./sum(is_BS+is_NT,3))
make_HR_coastlines([.6 .6 .6]);
set(gca,'clim',[0 .5])
colormap(Ax{1},brewermap(20,'Reds'))
colorbar
title('Fraction of CDR-SIC from BT','interpreter','latex')
%% 

Ax{2} = subplot('position',[.55 .55 .225 .4]);

worldmap([-90 -55],[-180 180]);

pcolorm(lat_SH,lon_SH,mean(abs(del_BS),3,'omitnan'))
make_HR_coastlines([.6 .6 .6]);
set(gca,'clim',[0 1])
colorbar
title('Fraction of CDR-SIC from BT','interpreter','latex')
title('MAE(CDR-BT,AMSR2-BT)','interpreter','latex')
colormap(brewermap(40,'Reds'))


Ax{3} = subplot('position',[.55 .05 .225 .4]);

worldmap([-90 -55],[-180 180]);

pcolorm(lat_SH,lon_SH,sum(high_del,3,'omitnan')./(sum(low_del,3,'omitnan') + sum(high_del,3,'omitnan')))
make_HR_coastlines([.6 .6 .6]);
set(gca,'clim',[0 1])

colorbar
title('Fraction of CDR-SIC from BT','interpreter','latex')
title('Frac with bias above 0.15','interpreter','latex')

%%
Ax{4} = subplot('position',[.775 .55 .225 .4]);

worldmap([-90 -55],[-180 180]);

pcolorm(lat_SH,lon_SH,mean(abs(del_NT2),3,'omitnan'))
make_HR_coastlines([.6 .6 .6]);
set(gca,'clim',[0 .5])
colormap(Ax{2},brewermap(40,'Blues'))
colorbar
title('Fraction of CDR-SIC from BT','interpreter','latex')
title('MAE(CDR-BT,AMSR2-BT)','interpreter','latex')


Ax{5} = subplot('position',[.775 .05 .225 .4]);

worldmap([-90 -55],[-180 180]);

pcolorm(lat_SH,lon_SH,sum(high_del_NT2,3,'omitnan')./(sum(low_del_NT2,3,'omitnan') + sum(high_del_NT2,3,'omitnan')))
make_HR_coastlines([.6 .6 .6]);
set(gca,'clim',[0 1])

colorbar
title('Frac NT2-BT above 0.15','interpreter','latex')

%%


pos = [6.5 4];
set(gcf,'windowstyle','normal','position',[0 0 pos],'paperposition',[0 0 pos],'papersize',pos,'units','inches','paperunits','inches');
set(gcf,'windowstyle','normal','position',[0 0 pos],'paperposition',[0 0 pos],'papersize',pos,'units','inches','paperunits','inches');
print([OPTS.plot_save_str 'SI-comp-PM-algos'],'-dpdf','-r600');
