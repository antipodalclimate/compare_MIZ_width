% compare sentinel and IS2

addpath('~/Dropbox (Brown)/Research Projects/Plot-Tools/')
addpath('~/Dropbox (Brown)/Research Projects/Plot-Tools/NE_Coastlines/')

S1_fold = '~/Dropbox-Brown/Christopher Horvat/Research Projects/Active/Data/Sentinel-1/Tavri_Classified/';

S1_files = dir([S1_fold '.tif']); 




% S1_file = 'S1B_EW_GRDM_1SSH_20210303T024209_20210303T024309_025844_03150B_B667_class_nowaves.tif';


class_KT = imread([S1_fold S1_file]);


%%
lat_KT = flipud(fliplr(class_KT(:,:,2))); 
lon_KT = flipud(fliplr(class_KT(:,:,3)));
class_KT = flipud(fliplr(class_KT(:,:,1))); 
amp_KT = flipud(fliplr(imread([S1_fold_AT 'S1A_EW_GRDM_1SSH_20190225_new_manual_binary.tif'])));

