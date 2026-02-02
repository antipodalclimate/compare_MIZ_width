% create_MIZ_width_panel

usable_beams = track_ID(intersections);


usable_waves = wave_flag(intersections) > 0 & WAF_width_dist(intersections) > 0 & MIZ_width_dist_CDR(intersections) > 0 & MIZ_width_dist_amsr(intersections) > 0; 

usable_nowaves = MIZ_width_dist_CDR(intersections) > 0 & MIZ_width_dist_amsr(intersections) > 0; 

beams_all = usable_beams(usable_waves | usable_nowaves);
beams_waves = usable_beams(usable_waves);
beams_nowaves = usable_beams(usable_nowaves);

%%

Wbins = linspace(0,250,25);
ratbins = logspace(-2,1,20);

wCDR = MIZ_width_CDR(usable_beams);
wAMSR = MIZ_width_amsr(usable_beams);
wWAF = WAF_width(usable_beams); 

wCDR_all = MIZ_width_dist_CDR(beams_all)/1000;
wAMSR_all = MIZ_width_dist_amsr(beams_all)/1000;
wWAF_all = WAF_width_dist(beams_all)/1000; 

wCDR_wav = MIZ_width_dist_CDR(beams_waves)/1000;
wAMSR_wav = MIZ_width_dist_amsr(beams_waves)/1000;
wWAF_wav = WAF_width_dist(beams_waves)/1000; 

wCDR_nowav = MIZ_width_dist_CDR(beams_nowaves)/1000;
wAMSR_nowav = MIZ_width_dist_amsr(beams_nowaves)/1000;
wWAF_nowav = WAF_width_dist(beams_nowaves)/1000; 


subplot('position',[.1 .1 .8 .15])

histogram(wCDR_wav,Wbins,'FaceColor',[.8 .2 .2],'FaceAlpha',0.5);
hold on

histogram(wAMSR_wav,Wbins,'FaceColor',[.2 .2 .8],'FaceAlpha',.5)
grid on; box on; 
% 
histogram(wWAF_wav,Wbins,'FaceColor',[.2 .8 .2],'FaceAlpha',.5)



legend('CDR','AMSR2','WAF')
hold off
xlim([0 max(Wbins)])
xlabel('km','interpreter','latex')
title('MIZ Width','interpreter','latex')

% %%
% subplot('position',[.55 .05 .35 .2])
% 
% 
% histogram(wCDR_wav,Wbins,'FaceColor',[.8 .2 .2],'FaceAlpha',0.5);
% hold on
% 
% histogram(wAMSR_wav,Wbins,'FaceColor',[.2 .2 .8],'FaceAlpha',.5)
% grid on; box on; 
% 
% histogram(wWAF_wav,Wbins,'FaceColor',[.2 .8 .2],'FaceAlpha',.5)


%%
