function param_plot(Bincent,N,LIF_med,LIF_up,LIF_dn,AMSR_med,AMSR_up,AMSR_dn,xlimmer,ylimmer)


hold on


jbfill(Bincent,LIF_up',LIF_dn',[.4 .4 .8],[1 1 1],1,.2);
hold on
jbfill(Bincent,AMSR_up',AMSR_dn',[.8 .2 .2],[1 1 1],1,.3);
hold on
plot(Bincent,AMSR_med,'color',[.8 .2 .2],'linewidth',2)
hold on
plot(Bincent,LIF_med,'color',[.4 .4 .8],'linewidth',2)

xlim(xlimmer);
ylim(ylimmer)
drawnow

yline(0,'k','linewidth',1')
ylabel('$\Delta$ from CDR','interpreter','latex')

grid on; box on;

yyaxis right
set(gca,'ycolor','k','yticklabel','')
jbfill(Bincent,N/sum(N),0*N/sum(N),[.8 .8 .8],[0 0 0],1,.8);
xlim(xlimmer);

end