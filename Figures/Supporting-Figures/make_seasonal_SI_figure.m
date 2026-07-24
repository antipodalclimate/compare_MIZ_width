% Create figure for month-by-month parametric plot

close all

%%

% Want to plot the parametric mapping between things. 

%% 
xlimmer_C = [.4 1];
xlimmer_W = [0 .5];
xlimmer_H = [-.5 1.5];

ylimmer = [-.2 .4];
clim_del = [-.3 .3];

mostr = {'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'};

for i = 1:12


    %%
    subplot(4,3,i)
    
    clear *_by_*
    
    param_use = SICvals_CDR <= 1 & Nsegvals >= 0 & monthvals == i;
    get_param_plot_vectors; 

    yyaxis right
    set(gca,'ycolor','k','yticklabel','')
    jbfill(Bincent_C,nC/sum(nC),0*nC/sum(nC),[.8 .8 .8],[0 0 0],1,.4);
    xlim(xlimmer);

    yyaxis left
    hold on
    set(gca,'ycolor','k')

    jbfill(Bincent_C,bias_AMSR_by_C_up',bias_AMSR_by_C_dn',[.8 .2 .2],[1 1 1],1,.3);
    hold on 
    plot(Bincent_C,bias_AMSR_by_C,'color',[.8 .2 .2],'linewidth',2);

    xlim(xlimmer_C);
    ylim(ylimmer)
    drawnow

    yline(0,'k','linewidth',1')

    grid on; box on;

    hold on
    yyaxis left
    hold on
    xline(.8,'linewidth',1,'label','MIZ')
    plot(Bincent_C(Bincent_C > 0.15),fitted(Bincent_C(Bincent_C > 0.15)),'--b','linewidth',1);
    plot(Bincent_C,1 - Bincent_C,'--r','linewidth',1)
    plot(Bincent_C,-Bincent_C,'--r','linewidth',1)
    xlabel('CDR SIC','interpreter','latex')
    ylabel('$\Delta$ from CDR','interpreter','latex')

    title(mostr{i})

end


allAxesInFigure = findall(gcf,'type','axes');
letter = {'(c)','(b)','(a)','(d)','(e)','(f)','(g)','(e)','(c)'};

for i = 1:length(allAxesInFigure)
    
 posy = get(allAxesInFigure(i),'position');

    set(allAxesInFigure(i),'fontname','times','fontsize',8,'xminortick','on','yminortick','on')
    
    % annotation('textbox',[posy(1) - .05 posy(2)+posy(4)-.015 .025 .025], ...
    %     'String',letter{i},'LineStyle','none','FontName','Helvetica', ...
    %     'FontSize',8,'Tag','legtag');

end

pos = [6.5 6.5]; 
set(gcf,'windowstyle','normal','position',[0 0 pos],'paperposition',[0 0 pos],'papersize',pos,'units','inches','paperunits','inches');
set(gcf,'windowstyle','normal','position',[0 0 pos],'paperposition',[0 0 pos],'papersize',pos,'units','inches','paperunits','inches');
print([OPTS.plot_save_str 'SI_Seasonal_1D'],'-dpdf','-r600');

