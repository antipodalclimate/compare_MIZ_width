close all

figure; 
pos = [6.5 2]; 
set(gcf,'windowstyle','normal','position',[0 0 pos],'paperposition',[0 0 pos],'papersize',pos,'units','inches','paperunits','inches');
set(gcf,'windowstyle','normal','position',[0 0 pos],'paperposition',[0 0 pos],'papersize',pos,'units','inches','paperunits','inches');

set(gca,'fontname','times','fontsize',8,'xminortick','on','yminortick','on')

subplot('position',[.125 .15 .775 .775])
% jbfill(Bincent,Hup',Hdn',[.4 .4 .4],[1 1 1],1,.3)

jbfill(Bincent,Hup',Hdn',[.4 .4 .4],[1 1 1],1,.3);
hold on
r1 = plot(Bincent,Hvec,'k','linewidth',2);

plot(Bincent,Hup,'--k');
plot(Bincent,Hdn,'--k');
xlim(xlimmer)
% ylim([0 .5]);
set(gca,'ylim',[0 max(max(get(gca,'ylim')),0.5)])
ylabel('Ice Height','interpreter','latex')
% r2 = yline(.075,'--','color',[.8 .4 .4],'linewidth',1)
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

saveas(gcf,'H-vs-bias-thick.pdf'); 

%%
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
