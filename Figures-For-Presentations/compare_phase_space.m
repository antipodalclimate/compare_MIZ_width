close all

figure; 
pos = [4 3]; 
set(gcf,'windowstyle','normal','position',[0 0 pos],'paperposition',[0 0 pos],'papersize',pos,'units','inches','paperunits','inches');
set(gcf,'windowstyle','normal','position',[0 0 pos],'paperposition',[0 0 pos],'papersize',pos,'units','inches','paperunits','inches');

set(gca,'fontname','times','fontsize',8,'xminortick','on','yminortick','on')

subplot('position',[.125 .15 .775 .775])
% jbfill(Bincent,Hup',Hdn',[.4 .4 .4],[1 1 1],1,.3)

bias_bins = [-.5:.02:.5];
wave_bins = [0:.01:1];

wave_use = WAFvals > 0.075 & Dvals < 200; 
nowave_use = WAFvals < 0.075 & Dvals < 200; 
thick_use = Hvals > 0.15 & Dvals < 200; 
thin_use = Hvals < 0.1 & Dvals < 200; 

sumcolor = [.8 .4 .4];
wincolor = [.4 .4 .8];
shocolor = [.8 .8 .4];

histogram(-biasvals(wave_use),bias_bins,'FaceColor',sumcolor,'FaceAlpha',0.5,'Normalization','pdf','edgecolor','none');
hold on

histogram(-biasvals(nowave_use),bias_bins,'FaceColor',wincolor,'FaceAlpha',0.5,'Normalization','pdf','edgecolor','none');
grid on; box on; 
legend({'Wavy','Not Wavy'})
xlabel('CDR Bias','interpreter','latex')
set(gca,'yticklabel','')
saveas(gcf,'wave-influence.pdf')

clf

histogram(-biasvals(thick_use),bias_bins,'FaceColor',sumcolor,'FaceAlpha',0.5,'Normalization','pdf','edgecolor','none');
hold on

histogram(-biasvals(thin_use),bias_bins,'FaceColor',wincolor,'FaceAlpha',0.5,'Normalization','pdf','edgecolor','none');
grid on; box on; 
legend({'Thick','Thin'})
saveas(gcf,'thick-influence.pdf')
xlabel('CDR Bias','interpreter','latex')
set(gca,'yticklabel','')


%%


yyaxis right

plot(Bincent,bias_LIFvec,'-','color',[.2 .8 .8],'linewidth',2); 
set(gca,'ycolor',[.2 .8 .8])
ylabel('LIF bias','interpreter','latex')
xlim(xlimmer)

saveas(gcf,'H-vs-bias.pdf'); 


clf

subplot('position',[.125 .15 .775 .775])
% jbfill(Bincent,Hup',Hdn',[.4 .4 .4],[1 1 1],1,.3)

jbfill(Bincent,WAFup',WAFdn',[.4 .4 .4],[1 1 1],1,.3);
hold on
r1 = plot(Bincent,WAFvec,'k','linewidth',2);

plot(Bincent,WAFup,'--k');
plot(Bincent,WAFdn,'--k');
xlim(xlimmer)
% ylim([0 .5]);
set(gca,'ylim',[0 max(max(get(gca,'ylim')),0.25)])
ylabel('Wave Affected Fraction','interpreter','latex')
r2 = yline(.075,'--','color',[.4 .4 .4],'linewidth',1)
xline(0,'color',[.2 .2 .2],'linewidth',1)
grid on; box on;

xline(0,'label','CDR-defined MIZ','interpreter','latex','fontsize',8,'LabelOrientation','horizontal','LabelHorizontalAlignment','left');
xline(0,'label','CDR-defined CIZ','interpreter','latex','fontsize',8,'LabelOrientation','horizontal');
xline(-50,'--','color',[.8 .4 .4])


xline(50,'--','color',[.8 .4 .4])


yyaxis right

plot(Bincent,bias_LIFvec,'-','color',[.2 .8 .8],'linewidth',2); 
set(gca,'ycolor',[.2 .8 .8])
ylabel('LIF bias','interpreter','latex')
xlim(xlimmer)
saveas(gcf,'WAF-vs-bias.pdf'); 
