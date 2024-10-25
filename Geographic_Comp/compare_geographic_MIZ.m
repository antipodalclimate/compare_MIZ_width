
clear
close all

AMSR_loc = '/Users/chorvat/Brown Dropbox/Christopher Horvat/Research Projects/Active/Data/SIC-Data/AMSR2-NT/AMSR2_SIC_daily.mat'; 
SSMI_loc = '/Users/chorvat/Brown Dropbox/Christopher Horvat/Research Projects/Active/Data/SIC-Data/NSIDC-CDR/Daily/NSIDC-CDR_daily.mat'; 

load(AMSR_loc,'AMSR_datenum','AMSR_NT_SH','AMSR_BS_SH');
load(SSMI_loc,'CDR_SIC_SH','lat_SH','lon_SH','CDR_datenum','BS_SIC_SH','NT_SIC_SH');

%%

[mutual_datenum,ia,ib] = intersect(AMSR_datenum,cell2mat(CDR_datenum));

% NASATEAM2 AMSR values
AMSR_NT_SH = AMSR_NT_SH(:,:,ia);
% CDR values
CDR_SIC_SH = CDR_SIC_SH(:,:,ib);
SSMI_BS_SH = BS_SIC_SH(:,:,ib); 
SSMI_NT_SH = NT_SIC_SH(:,:,ib); 

% BOOTSTRAP AMSR values
AMSR_BS_SH =  AMSR_BS_SH(:,:,ia); 

%%

AMSR_NT_SH(AMSR_NT_SH > 1) = nan; 
AMSR_BS_SH(AMSR_BS_SH > 1) = nan; 

MIZ_AMSR_NT = AMSR_NT_SH > 0.15 & AMSR_NT_SH < 0.8; 
MIZ_AMSR_BS = AMSR_BS_SH > 0.15 & AMSR_BS_SH < 0.8; 
MIZ_CDR = CDR_SIC_SH > 0.15 & CDR_SIC_SH < 0.8; 
MIZ_SSMI_NT = SSMI_NT_SH > 0.15 & SSMI_NT_SH < 0.8; 
MIZ_SSMI_BS = SSMI_BS_SH > 0.15 & SSMI_BS_SH < 0.8; 

ICE_AMSR_NT = AMSR_NT_SH > 0.15; 
ICE_AMSR_BS = AMSR_BS_SH > 0.15; 
ICE_CDR = CDR_SIC_SH > 0.15; 
ICE_SSMI_NT = SSMI_NT_SH > 0.15; 
ICE_SSMI_BS = SSMI_BS_SH > 0.15; 

prefac = 25*25/1e6; 

AMIZ_AMSR_BS = prefac*squeeze(sum(MIZ_AMSR_BS,[1 2],'omitnan'));
AMIZ_AMSR_NT = prefac*squeeze(sum(MIZ_AMSR_NT,[1 2],'omitnan'));
AMIZ_CDR = prefac*squeeze(sum(MIZ_CDR,[1 2],'omitnan')); 
AMIZ_SSMI_NT = prefac*squeeze(sum(MIZ_SSMI_NT,[1 2],'omitnan'));
AMIZ_SSMI_BS = prefac*squeeze(sum(MIZ_SSMI_BS,[1 2],'omitnan'));


SIA_AMSR_NT = prefac*squeeze(sum(AMSR_NT_SH,[1 2],'omitnan'));
SIA_AMSR_BS = prefac*squeeze(sum(AMSR_BS_SH,[1 2],'omitnan'));
SIA_CDR = prefac*squeeze(sum(CDR_SIC_SH,[1 2],'omitnan'));
SIA_SSMI_NT = prefac*squeeze(sum(SSMI_NT_SH,[1 2],'omitnan'));
SIA_SSMI_BS = prefac*squeeze(sum(SSMI_BS_SH,[1 2],'omitnan'));

SIE_AMSR_NT = prefac*squeeze(sum(ICE_AMSR_NT,[1 2],'omitnan'));
SIE_AMSR_BS = prefac*squeeze(sum(ICE_AMSR_BS,[1 2],'omitnan'));
SIE_CDR = prefac*squeeze(sum(ICE_CDR,[1 2],'omitnan')); 
SIE_SSMI_NT = prefac*squeeze(sum(ICE_SSMI_NT,[1 2],'omitnan'));
SIE_SSMI_BS = prefac*squeeze(sum(ICE_SSMI_BS,[1 2],'omitnan'));

% Five products
% AMSR2 - NT
% AMSR2 - BS
% SSMI - CDR
% SSMI - BS
% SSMI - NT

% Interesected in official product differences
% AMSR2 - SSMI-CDR
dAMIZ = AMIZ_AMSR_NT - AMIZ_CDR; 

% Then interested in same algorithms, different sensors
dAMIZ_BS = AMIZ_AMSR_BS - AMIZ_SSMI_BS; 
dAMIZ_NT = AMIZ_AMSR_NT - AMIZ_SSMI_NT; 

% Then interested in same sensor, different algorithms
dAMIZ_AMSR = AMIZ_AMSR_NT - AMIZ_AMSR_BS; 
dAMIZ_SSMI = AMIZ_SSMI_NT - AMIZ_SSMI_BS; 

% Now for shits, different sensor, different algos
dAMIZ_cross = dAMIZ_NT + dAMIZ_SSMI; 

dSIE = SIE_AMSR_NT - SIE_CDR; 

% Now SIA differences
dSIA = SIA_AMSR_NT - SIA_CDR; 
dSIA_NT = SIA_AMSR_NT - SIA_SSMI_NT;

xax = datetime(datestr(mutual_datenum)); 

close 

xlimmer = [datetime('20210101','InputFormat','yyyyMMdd') datetime('20211231','InputFormat','yyyyMMdd')];

subplot(411)

plot(xax,AMIZ_AMSR_NT,'r','linewidth',1);
hold on
plot(xax,AMIZ_AMSR_BS,'--r','linewidth',1);
plot(xax,AMIZ_CDR,'k','linewidth',2);
plot(xax,AMIZ_SSMI_NT,'b','linewidth',1);
plot(xax,AMIZ_SSMI_BS,'--b','linewidth',1);

xlim(xlimmer)

subplot(423)
plot(xax,SIE_AMSR_NT,'r','linewidth',1);
hold on
plot(xax,SIE_AMSR_BS,'--r','linewidth',1);
plot(xax,SIE_CDR,'k','linewidth',2);
plot(xax,SIE_SSMI_NT,'b','linewidth',1);
plot(xax,SIE_SSMI_BS,'--b','linewidth',1);

xlim(xlimmer)

subplot(424)
plot(xax,SIA_AMSR_NT,'r','linewidth',1);
hold on
plot(xax,SIA_AMSR_BS,'--r','linewidth',1);
plot(xax,SIA_CDR,'k','linewidth',2);
plot(xax,SIA_SSMI_NT,'b','linewidth',1);
plot(xax,SIA_SSMI_BS,'--b','linewidth',1);


xlim(xlimmer)


subplot(413)
plot(xax,dAMIZ,'-k','linewidth',2); % Official product differences 
hold on
plot(xax,dAMIZ_BS,'--r','linewidth',2); % Same algo, different sensors
plot(xax,dAMIZ_NT,'--b','linewidth',2); % Same algo, different sensors
plot(xax,dAMIZ_AMSR,'-.r','linewidth',1); % Same sensor, different algo
plot(xax,dAMIZ_SSMI,'-.b','linewidth',1); % Same sensor, different algo

plot(xax,dAMIZ_cross,'--','color','m','linewidth',1)

xlim(xlimmer)

subplot(414)
plot(xax,dSIA,'-k','linewidth',2); % Official product differences 
hold on
plot(xax,dSIA_NT,'b','linewidth',1)
plot(xax,dAMIZ,'--k','linewidth',1)

xlim(xlimmer)

pos = [9 10]; 
set(gcf,'windowstyle','normal','position',[0 0 pos],'paperposition',[0 0 pos],'papersize',pos,'units','inches','paperunits','inches');
set(gcf,'windowstyle','normal','position',[0 0 pos],'paperposition',[0 0 pos],'papersize',pos,'units','inches','paperunits','inches');
