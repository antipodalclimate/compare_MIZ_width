function [p1,p2,p3,p4] = param_plot(Bincent,N,LIF_med,LIF_up,LIF_dn,AMSR_med,AMSR_up,AMSR_dn,CDR_med,CDR_up,CDR_dn,xlimmer,ylimmer)

yyaxis right
set(gca,'ycolor','k','yticklabel','')
p4 = jbfill(Bincent,N/sum(N),0*N/sum(N),[.8 .8 .8],[0 0 0],1,.4);
xlim(xlimmer);

yyaxis left
hold on
set(gca,'ycolor','k')

jbfill(Bincent,LIF_up',LIF_dn',[.4 .4 .8],[1 1 1],1,.2);
hold on
jbfill(Bincent,AMSR_up',AMSR_dn',[.8 .2 .2],[1 1 1],1,.3);
hold on
p3 = jbfill(Bincent,(CDR_up-CDR_med)',(CDR_dn - CDR_med)',[217 95 2]/256,[1 1 1],1,.3);
hold on
p1 = plot(Bincent,AMSR_med,'color',[.8 .2 .2],'linewidth',2);
hold on
p2 = plot(Bincent,LIF_med,'-','color',[.4 .4 .8],'linewidth',2);

xlim(xlimmer);
ylim(ylimmer)
drawnow

yline(0,'k','linewidth',1')

grid on; box on;



end