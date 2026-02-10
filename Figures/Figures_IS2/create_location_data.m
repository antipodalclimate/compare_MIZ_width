%% We want to analyse the location of each of the points we are comparing. 

load(fullfile(OPTS.code_folder,'Data','KDTrees','KDTree_100km.mat'),'lat_SH','lon_SH','KDTree','area_SH');


%%


fprintf(' Locating Grid Indices - ');
% Search for the grid values for each point via their lat/lon


ALL_posloc = knnsearch(KDTree,[latvals(keepvals==1) lonvals(keepvals==1)],'K',1);

MIZ_posloc = knnsearch(KDTree,[latvals(keepvals==1 & num_MIZ > 0) lonvals(keepvals==1 & num_MIZ > 0)],'K',1);

used_posloc = knnsearch(KDTree,[latvals(usable) lonvals(usable)],'K',1);


fprintf('Done \n')

num_stencils_all = accumarray(ALL_posloc,1 + 0*length(latvals),[numel(KDTree.X(:,1)) 1],@sum); 
num_stencils_all = reshape(num_stencils_all(end-numel(lat_SH(:))+1:end),size(lat_SH));

num_stencils_MIZ = accumarray(MIZ_posloc,1 + 0*length(latvals),[numel(KDTree.X(:,1)) 1],@sum); 
num_stencils_MIZ = reshape(num_stencils_MIZ(end-numel(lat_SH(:))+1:end),size(lat_SH));

num_stencils_usable = accumarray(used_posloc,1 + 0*length(latvals),[numel(KDTree.X(:,1)) 1],@sum); 
num_stencils_usable = reshape(num_stencils_usable(end-numel(lat_SH(:))+1:end),size(lat_SH));

%%
close all

clim_num = [1 5]

Ax{1} = subplot('position',[0 .05 .25 .85]);

worldmap([-90 -55],[-180 180]);
pcolorm(lat_SH,lon_SH,log10(num_stencils_all))
make_HR_coastlines([.6 .6 .6]);
% set(gca,'clim',climmer)
% set(gca,'clim',[-.2 .2])
% colorbar('position',[.55 .55 .025 .35]);
title('All Stencils','interpreter','latex')
set(gca,'clim',clim_num)

Ax{2} = subplot('position',[.3 .05 .25 .85]);
worldmap([-90 -55],[-180 180]);
pcolorm(lat_SH,lon_SH,log10(num_stencils_MIZ))
make_HR_coastlines([.6 .6 .6]);
% set(gca,'clim',climmer)
% set(gca,'clim',[-.2 .2])
% colorbar('position',[.575 .15 .025 .7]);
title('With an MIZ point','interpreter','latex')
set(gca,'clim',clim_num)


Ax{3} = subplot('position',[.6 .05 .25 .85]);
worldmap([-90 -55],[-180 180]);
pcolorm(lat_SH,lon_SH,log10(num_stencils_usable))
make_HR_coastlines([.6 .6 .6]);
set(gca,'clim',clim_num)

title('Crossing the MIZ','interpreter','latex')

cmapper = brewermap(5,'Dark2');
colormap(cmapper);



colorbar('position',[.925 .15 .025 .7]);

letter = {'(a)','(b)','(c)','(d)','(e)','(F)','(g)','(e)','(c)'};

delete(findall(gcf,'Tag','legtag'))

for i = 1:length(Ax)
    set(Ax{i},'fontname','helvetica','fontsize',9,'xminortick','on','yminortick','on')
    posy = get(Ax{i},'position');
    annotation('textbox',[posy(1) posy(2)+posy(4)-.025 .025 .025], ...
        'String',letter{i},'LineStyle','none','FontName','Helvetica', ...
        'FontSize',8,'Tag','legtag');

end

pos = [6.5 2.5];
set(gcf,'windowstyle','normal','position',[0 0 pos],'paperposition',[0 0 pos],'papersize',pos,'units','inches','paperunits','inches');
set(gcf,'windowstyle','normal','position',[0 0 pos],'paperposition',[0 0 pos],'papersize',pos,'units','inches','paperunits','inches');
print([OPTS.plot_save_str 'location-data'],'-dpdf','-r600');
