close all

figure; 
pos = [4 3]; 
set(gcf,'windowstyle','normal','position',[0 0 pos],'paperposition',[0 0 pos],'papersize',pos,'units','inches','paperunits','inches');
set(gcf,'windowstyle','normal','position',[0 0 pos],'paperposition',[0 0 pos],'papersize',pos,'units','inches','paperunits','inches');

set(gca,'fontname','times','fontsize',8,'xminortick','on','yminortick','on')

subplot('position',[.125 .15 .775 .775])
% jbfill(Bincent,Hup',Hdn',[.4 .4 .4],[1 1 1],1,.3)

bias_bins = [0:.02:1];
wave_bins = [0:.01:1];

allbias = biasvals_LIF > -1; 
highbias = biasvals_LIF > 0.2; 

sumcolor = [.8 .4 .4];
wincolor = [.4 .4 .8];
shocolor = [.8 .8 .4];


subplot('position',[.125 .15 .775 .775])

histogram(LIFvals,bias_bins,'FaceColor',wincolor,'FaceAlpha',0.5,'Normalization','pdf','edgecolor','none');
hold on
histogram(LIFvals(highbias),bias_bins,'FaceColor',sumcolor,'FaceAlpha',0.5,'Normalization','pdf','edgecolor',sumcolor);
grid on; box on; 
xlabel('LIF','interpreter','latex')

legend('All','High Bias')
set(gca,'yticklabel','')
saveas(gcf,'LIF-highbias.pdf')
