clear

addpath('/Users/chorvat/Library/CloudStorage/Dropbox-Brown/Christopher Horvat/Research Projects/Plot-Tools/')
addpath('/Users/chorvat/Library/CloudStorage/Dropbox-Brown/Christopher Horvat/Research Projects/Plot-Tools/NE_Coastlines/')

S1_fold_AT = '/Users/chorvat/Library/CloudStorage/Dropbox-Brown/Christopher Horvat/Research Projects/Active/Data/Sentinel-1/Tavri_Classified/';

S1_fold = '/Users/chorvat/Library/CloudStorage/Dropbox-Brown/Christopher Horvat/Research Projects/Active/Data/Sentinel-1/';

S1_filenames = dir([S1_fold '*.nc']);
S1_file = S1_filenames(2).name;



% S1_file = 'S1A_EW_GRDM_1SSH_20190225T014506_20190225T014611_026079_02E89C_56F1.nc';
% S1_file = 'S1A_EW_GRDM_1SSH_20190226T022533_20190226T022638_026094_02E936_7558.nc';

lat_S1 = ncread([S1_fold S1_file],'lat'); 
lon_S1 = ncread([S1_fold S1_file],'lon'); 

S1_lat_span = [min(lat_S1(:)) max(lat_S1(:))];
S1_lon_span = [min(lon_S1(:)) max(lon_S1(:))];


%%
Amp = ncread([S1_fold S1_file],'Amplitude_HH'); 
Int = ncread([S1_fold S1_file],'Intensity_HH'); 

%%
class_KT = imread([S1_fold_AT 'ThurSAR__corrected_mask_with_latlon.tiff']);
lat_KT = flipud(fliplr(class_KT(:,:,2))); 
lon_KT = flipud(fliplr(class_KT(:,:,3)));
class_KT = flipud(fliplr(class_KT(:,:,1))); 
amp_KT = flipud(fliplr(imread([S1_fold_AT 'S1A_EW_GRDM_1SSH_20190225_new_manual_binary.tif'])));


%% IS-2 track
mov_window = [12500 12500];


IS2_fold = '/Users/chorvat/Library/CloudStorage/Dropbox-Brown/Christopher Horvat/Research Projects/Active/ICESAT-2/MIZ-Definitions/Fig_Example_Track/';
IS2_filenames = dir([IS2_fold '*.h5']); 
% 
IS2_file = [IS2_fold '/' IS2_filenames(1).name];
% IS2_file = 'ATL07-02_20190225005459_08950201_005_01.h5'; 
beamname = '/gt3l';
% addpath(IS2_fold);
IS2_obj = preprocess_single_track(IS2_file,beamname);
AT_stats = get_AT_variables(IS2_obj);

track_span = [0 150];
end_span = find(IS2_obj.dist/1000 > max(track_span),1); 
lat_span = sort([IS2_obj.lat(1) IS2_obj.lat(end_span)]); 
lon_span = sort([IS2_obj.lon(1) IS2_obj.lon(end_span)]); 


%% Find locations of IS2 on the SAR data
inds = lat_S1(:) < lat_span(2) & lat_S1(:) > lat_span(1) & lon_S1(:) < lon_span(2) & lon_S1(:) > lon_span(1);
lat_X = lat_S1(inds);  
lon_X = lon_S1(inds);  
M = createns([lat_X,lon_X]); 

% Nearest neighbor of the IS2 lat/lon on the S1 lat/lon. Might take a
% while. 
ID = knnsearch(M,[IS2_obj.lat IS2_obj.lon],'K',1); 

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

%% Look at AT-classified image

% Take the indices in the image which are in the right range
inds_KT = lat_KT(:) < lat_span(2) & lat_KT(:) > lat_span(1) ... 
    & lon_KT(:) < lon_span(2) & lon_KT(:) > lon_span(1);

usable_class_KT = class_KT(inds_KT); 

lat_X_KT = lat_KT(inds_KT);  
lon_X_KT = lon_KT(inds_KT);  
M_KT = createns([lat_X_KT,lon_X_KT]); 




%%
% Nearest neighbor of the IS2 lat/lon on the S1 lat/lon. Might take a
% while. 
[ID_KT,D] = knnsearch(M_KT,[IS2_obj.lat IS2_obj.lon],'K',1); 

found_lats = lat_X_KT(ID_KT);
found_lons = lon_X_KT(ID_KT); 


class_AT_AT = usable_class_KT(ID_KT); 


S1_LIF = movmean(class_AT_AT,mov_window,'omitnan','samplepoints',IS2_obj.dist);


%%


% Remove outlier values
AT_stats.height_adj(AT_stats.height_adj > 1.5) = nan; 

height_smooth = movsum(AT_stats.height_adj(IS2_obj.is_ice).*IS2_obj.seg_len(IS2_obj.is_ice),mov_window,'omitnan','samplepoints',IS2_obj.dist(IS2_obj.is_ice)) ...
    ./ movsum(IS2_obj.seg_len(IS2_obj.is_ice),mov_window,'omitnan','samplepoints',IS2_obj.dist(IS2_obj.is_ice));

height_var = movsum((AT_stats.height_adj(IS2_obj.is_ice) - height_smooth).^2 .* IS2_obj.seg_len(IS2_obj.is_ice),mov_window,'omitnan','samplepoints',IS2_obj.dist(IS2_obj.is_ice)) ...
    ./ movsum(IS2_obj.seg_len(IS2_obj.is_ice),mov_window,'omitnan','samplepoints',IS2_obj.dist(IS2_obj.is_ice));

height_std = sqrt(height_var); 

%
xvals = IS2_obj.dist(IS2_obj.is_ice)/1000; 
latice = IS2_obj.lat(IS2_obj.is_ice); 
lonice = IS2_obj.lon(IS2_obj.is_ice); 

[X_WAF,X_PM,X_S1,X_AT,X_S1_KT,X_AT_adj] = deal(0); 

ind_WAF = find(AT_stats.WAF(IS2_obj.is_ice) < 0.075,1);
if ~isempty(ind_WAF)
X_WAF= xvals(ind_WAF);
end
ind_PM = find(AT_stats.SIC(IS2_obj.is_ice) > .8,1);
if ~isempty(ind_PM)
X_PM = xvals(ind_PM);
end
ind_S1 = find(amp_smooth(IS2_obj.is_ice) > .8,1);
if ~isempty(ind_S1)
X_S1 = xvals(ind_S1);
end

ind_S1_KT = find(S1_LIF(IS2_obj.is_ice) > .8,1);
if ~isempty(ind_S1_KT)
X_S1_KT = xvals(ind_S1_KT);
end

ind_AT = find(AT_stats.LIF(IS2_obj.is_ice) > .8,1);
if ~isempty(ind_AT)
X_AT = xvals(ind_AT);
end

ind_AT_adj = find(AT_stats.LIF_adj(IS2_obj.is_ice) > .8,1);
if ~isempty(ind_AT_adj)
X_AT = xvals(ind_AT_adj);
end

%
close all
addpath('Plot-Tools')
horvat_colors; 
Ax{1} = subplot('position',[.05 .455 .45 .5]);

latlim = sort([0.5 .25] + sort([IS2_obj.lat(1) IS2_obj.lat(end_span)]));

lonlim = sort([-1.75 2.25] + sort([IS2_obj.lon(1) IS2_obj.lon(end_span)]));

% if lonlim(2) > max(lon_S1(:))
%     lonlim(2) = max(lon_S1(:))
% end

worldmap((latlim),(lonlim)); 
setm(gca,'PlabelLocation',-70,'PlineLocation',-70,'PLabelRound',0);
setm(gca,'MLabelLocation',[-84 -82],'MLineLocation',[-84 -82],'MLabelRound',0,'MLabelParallel','south');
setm(gca,'grid','on','GLineWidth',0.5,'GLineStyle','--','GColor','w','fontsize',8);

% setm(gca,'maplatlimit',latlim,'parallellabel','off','meridianlabel','off');
% setm(gca,'plinevisible','off','mlinevisible','off'); 
% rgbmat = zeros([size(Amp) 3]);
% rgbmat(:,:,1) = Amp; 
% rgbmat(:,:,2) = Int; 

skipper = 4; 

inds_plot = lat_S1(:) < latlim(2) & lat_S1(:) > latlim(1) & lon_S1(:) < lonlim(2) & lon_S1(:) > lonlim(1);
Amp_plot = nan*Amp2; 
Amp_plot(inds_plot) = Amp2(inds_plot); 
Amp_plot(Amp_plot > 1000) = 1000; 
Amp_plot = Amp_plot / 1000; 

[a,b] = ind2sub(size(Amp),find(inds_plot));
l1 = min(a):skipper:max(a); 
l2 = min(b):skipper:max(b); 

pcolorm(lat_S1(l1,l2),lon_S1(l1,l2),Amp_plot(l1,l2));
set(gca,'clim',[0.5 1])
colormap((brewermap(9,'-greys')))
plotm(latice,lonice,'--r','linewidth',1)

s= scatterm(latice(1),lonice(1),20,'k','s','linewidth',3);
scatterm(latice(ind_WAF),lonice(ind_WAF),100,clabs(1,:),'s','linewidth',1.5);
scatterm(latice(ind_AT),lonice(ind_AT),100,clabs(2,:),'s','linewidth',1.5);
scatterm(latice(ind_PM),lonice(ind_PM),100,clabs(3,:),'s','linewidth',1.5);

% scatterm(latice(ind_S1),lonice(ind_S1),100,clabs(4,:),'s','linewidth',1.5);
scatterm(latice(ind_S1_KT),lonice(ind_S1_KT),100,clabs(4,:),'s','linewidth',1.5);


% add_coastlines; 
%
Ax{2} = subplot('position',[.1 .15 .8 .2],'replace');
plot(xvals,100*AT_stats.height_adj(IS2_obj.is_ice),'--','linewidth',0.05,'color',[.7 .7 .7])
hold on
plot(xvals,100*height_smooth,'k','linewidth',1)
plot(xvals,100*height_smooth + 100*height_std,'--k','linewidth',1)
plot(xvals,100*height_smooth - 100*height_std,'--k','linewidth',1)

ylim([-.2 .5]*100)
xlim([0 100]); 

xline(X_WAF,'color',clabs(1,:),'linewidth',3)
xline(X_AT,'color',clabs(2,:),'linewidth',3)
xline(X_PM,'color',clabs(3,:),'linewidth',3)
% xline(X_S1,'color',clabs(4,:),'linewidth',3)
xline(X_S1_KT,'color',clabs(4,:),'linewidth',3)

grid on; box on; 
ylabel('cm','interpreter','latex');
xlabel('Distance from 1st ice segment','interpreter','latex');

%
Ax{3} = subplot('position',[.525 .5 .35 .45],'replace');


plot(xvals,AT_stats.WAF(IS2_obj.is_ice),'linewidth',1,'color',clabs(1,:))
hold on
plot(xvals,AT_stats.LIF_adj(IS2_obj.is_ice),'linewidth',1,'color',clabs(2,:))
% plot(xvals,AT_stats.LIF(IS2_obj.is_ice),'linewidth',1,'color',clabs(2,:))
plot(xvals,AT_stats.SIC(IS2_obj.is_ice),'linewidth',1,'color',clabs(3,:))

% plot(xvals,amp_smooth(IS2_obj.is_ice),'linewidth',1,'color',clabs(4,:))
plot(xvals,S1_LIF(IS2_obj.is_ice),'linewidth',1,'color',clabs(4,:))


xlim(track_span); 
grid on; box on; 
xlabel('Distance from 1st ice segment','interpreter','latex');

% yyaxis right
% set(gca,'ycolor','k')
% plot(xvals,100*height_smooth,'k','linewidth',1)
xlim([0 100]); 
grid on; box on; 
ylabel('cm','interpreter','latex');
% plot(xvals,100*height_smooth + 100*height_std,'--k','linewidth',1)
% plot(xvals,100*height_smooth - 100*height_std,'--k','linewidth',1)
xline(X_WAF,'color',clabs(1,:),'linewidth',3)
xline(X_AT,'color',clabs(2,:),'linewidth',3)
xline(X_PM,'color',clabs(3,:),'linewidth',3)
% xline(X_S1_,'color',clabs(4,:),'linewidth',3)
xline(X_S1_KT,'color',clabs(4,:),'linewidth',3)
h = legend('AT-WAF','AT-LIF','PM-SIC','SAR-SIC','location','southeast');
set(h,'ItemTokenSize',[20 20]);
%

letter = {'(A)','(C)','(B)','(d)','(e)','(f)','(g)','(e)','(c)'};

delete(findall(gcf,'Tag','legtag'))

for i = 1:length(Ax)
    set(Ax{i},'fontname','helvetica','fontsize',9,'xminortick','on','yminortick','on')
    posy = get(Ax{i},'position');
    annotation('textbox',[posy(1)-.04 posy(2)+posy(4)+.04 .025 .025], ...
        'String',letter{i},'LineStyle','none','FontName','Helvetica', ...
        'FontSize',8,'Tag','legtag');
    
end


pos = [6.5 3.25]; 
set(gcf,'windowstyle','normal','position',[0 0 pos],'paperposition',[0 0 pos],'papersize',pos,'units','inches','paperunits','inches');
set(gcf,'windowstyle','normal','position',[0 0 pos],'paperposition',[0 0 pos],'papersize',pos,'units','inches','paperunits','inches');
print('/Users/chorvat/Library/CloudStorage/Dropbox-Brown/Christopher Horvat/Apps/Overleaf/IS2-Waves-PM/Figures/SAR-Comp','-dpdf','-r600');
print('/Users/chorvat/Library/CloudStorage/Dropbox-Brown/Christopher Horvat/Apps/Overleaf/2024-NASA-ROSES-PM/Proposal/Figures/SAR-Comp','-dpdf','-r600');

