load_MIZ_waves; 


usable_all = (Nvals > 10) & usable_all; 
usable_all = usable_all & SICvals > 0.1 & LIFvals > 0.1;
usable_all = usable_all; %  & (timeval > 7 & timeval < 10); 
usable_all = usable_all & (isstrong == 1); 
usable_all = usable_all & npoints >= 0; 

%%
usable_tracks = (track_npoints >= 5) & track_is_strong == 1 & track_nMIZ > 0;  

wave_MIZ_widths = WAF_width(usable_tracks);
SSMI_MIZ_widths = MIZ_width(usable_tracks);
AMSR_MIZ_widths = MIZ_width_amsr(usable_tracks);

%%
close

subplot('position',[.1 .55 .35 .35])
histogram(SSMI_MIZ_widths(wave_MIZ_widths == 0),[0:100],'FaceColor',[.8 .2 .2],'FaceAlpha',.5); 
hold on
histogram(SSMI_MIZ_widths(wave_MIZ_widths > 0),[0:100],'FaceColor',[.2 .2 .8],'FaceAlpha',.5); 
hold off
legend('No Waves','Waves')
grid on; box on; 
xlim([0 25])
title('SSMI/S');

subplot('position',[.55 .55 .35 .35])


histogram(AMSR_MIZ_widths(wave_MIZ_widths == 0),[0:100],'FaceColor',[.8 .2 .2],'FaceAlpha',.5); 
hold on
histogram(AMSR_MIZ_widths(wave_MIZ_widths > 0),[0:100],'FaceColor',[.2 .2 .8],'FaceAlpha',.5); 
hold off
legend('No Waves','Waves')
grid on; box on; 
xlim([0 25])
title('AMSR2');

subplot('position',[.1 .1 .35 .35])

wavefac = (AMSR_MIZ_widths)./SSMI_MIZ_widths; 

% wavefac(wavefac < 0) = []; 
wavefac(isinf(wavefac)) = []; 

histogram(log10(wavefac)); 