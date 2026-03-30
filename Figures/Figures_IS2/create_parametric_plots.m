

% Want to plot the parametric mapping between things. 
clear *_by_*

close all

param_use = SICvals_CDR <= 1 & Nsegvals >= 0;

get_param_plot_vectors; 


%% Here we want to examine the biases using all SIC values3
clf

xlimmer_C = [.4 1];
xlimmer_W = [0 .5];
xlimmer_H = [-.5 1.5];

ylimmer = [-.2 .4];
clim_del = [-.3 .3];

% Biases as a function of CDR SIC 
Ax{1} = subplot('position',[.1 .15 .275 .7]);

[p1,p2,p3,p4] = param_plot(Bincent_C,nC,bias_LIF_by_C,bias_LIF_by_C_up,bias_LIF_by_C_dn,bias_AMSR_by_C,bias_AMSR_by_C_up,bias_AMSR_by_C_dn, ...
    CDR_by_C,CDR_by_C_up,CDR_by_C_dn,xlimmer_C,ylimmer);

hold on
yyaxis left 
hold on
xline(.8,'linewidth',1,'label','MIZ')
plot(Bincent_C(Bincent_C > 0.15),fitted(Bincent_C(Bincent_C > 0.15)),'--b','linewidth',1);
plot(Bincent_C,1 - Bincent_C,'--r','linewidth',1)
plot(Bincent_C,-Bincent_C,'--r','linewidth',1)
xlabel('CDR SIC','interpreter','latex')
ylabel('$\Delta$ from CDR','interpreter','latex')
legend([p1 p2 p3 p4],'AMSR offset','LIF offset','CDR range','Count','orientation','horizontal','location',[.3 .9 .4 .045])

Ax{2} = subplot('position',[.4 .15 .275 .7]);

param_plot(Bincent_H,nH,bias_LIF_by_H,bias_LIF_by_H_up,bias_LIF_by_H_dn,bias_AMSR_by_H,bias_AMSR_by_H_up,bias_AMSR_by_H_dn, ...
    CDR_by_H,CDR_by_H_up,CDR_by_H_dn,xlimmer_H,ylimmer);
hold on
xlabel('IS2 Height','interpreter','latex')
yyaxis left 
ylabel('')
set(gca,'yticklabel','')


Ax{3} = subplot('position',[.7 .15 .25 .7]);
param_plot(Bincent_W,nW,bias_LIF_by_W,bias_LIF_by_W_up,bias_LIF_by_W_dn,bias_AMSR_by_W,bias_AMSR_by_W_up,bias_AMSR_by_W_dn, ...
    CDR_by_W,CDR_by_W_up,CDR_by_W_dn,xlimmer_W,ylimmer);
hold on
xlabel('IS2 WAF','interpreter','latex')
yyaxis left 
ylabel('')
set(gca,'yticklabel','')
yyaxis right
ylabel('Count','interpreter','latex')
xline(.075,'label','Horvat et al (2020) WMIZ','color',[.4 .4 .8],'linewidth',1,'interpreter','latex','fontsize',8,'LabelOrientation','horizontal');



allAxesInFigure = findall(gcf,'type','axes');
letter = {'(c)','(b)','(a)','(d)','(e)','(f)','(g)','(e)','(c)'};

for i = 1:length(allAxesInFigure)
    
 posy = get(allAxesInFigure(i),'position');

    set(allAxesInFigure(i),'fontname','times','fontsize',8,'xminortick','on','yminortick','on')
    
    annotation('textbox',[posy(1) posy(2)+posy(4)-.015 .025 .025], ...
        'String',letter{i},'LineStyle','none','FontName','Helvetica', ...
        'FontSize',8,'Tag','legtag');

end

pos = [6.5 2.5]; 
set(gcf,'windowstyle','normal','position',[0 0 pos],'paperposition',[0 0 pos],'papersize',pos,'units','inches','paperunits','inches');
set(gcf,'windowstyle','normal','position',[0 0 pos],'paperposition',[0 0 pos],'papersize',pos,'units','inches','paperunits','inches');
print([OPTS.plot_save_str '1D_Parametric_Comp'],'-dpdf','-r600');


%% 
close all

Ax{1} = subplot('position',[.075 .55 .1375 .35]);

pcolor(Bincent_H,Bincent_C, bias_LIF_by_HC);
set(gca,'xlim',xlimmer_H);
set(gca,'ylim',xlimmer_C)
xlabel('IS2 Height','interpreter','latex')
ylabel('CDR SIC','interpreter','latex')
set(gca,'clim',[-.5 .5])
shading flat
grid on; box on; 
title('LIF offset','interpreter','latex')

Ax{2} = subplot('position',[.225 .55 .1375 .35]);

pcolor(Bincent_H,Bincent_C, bias_AMSR_by_HC);
set(gca,'xlim',xlimmer_H);
set(gca,'ylim',xlimmer_C)
xlabel('IS2 Height','interpreter','latex')
set(gca,'clim',[-.5 .5])
set(gca,'yticklabel','')
shading flat
grid on; box on; 
title('AMSR2-NT2 offset','interpreter','latex')
colorbar('position',[.375 .55 .0125 .35])


Ax{3} = subplot('position',[.425 .55 .1375 .35]);
pcolor(Bincent_H,Bincent_C, 100*dark_by_HC);
set(gca,'xlim',xlimmer_H);
set(gca,'ylim',xlimmer_C)
set(gca,'yticklabel','')
% xlabel('IS2 Height','interpreter','latex')
% ylabel('CDR SIC','interpreter','latex')
colormap(brewermap(25,'Blues'))
set(gca,'clim',[0 25])
shading flat
grid on; box on; 
title('Dark Lead \%','interpreter','latex')
xlabel('IS2 Height','interpreter','latex')

Ax{4} = subplot('position',[.575 .55 .1375 .35]);
pcolor(Bincent_H,Bincent_C, 100*spec_by_HC);
set(gca,'xlim',xlimmer_H);
set(gca,'ylim',xlimmer_C)
set(gca,'yticklabel','')
% xlabel('IS2 Height','interpreter','latex')
% ylabel('CDR SIC','interpreter','latex')
colormap(brewermap(25,'Blues'))
set(gca,'clim',[0 25])
shading flat
grid on; box on; 
title('Specular Lead \%','interpreter','latex')
xlabel('IS2 Height','interpreter','latex')
colorbar('position',[.725 .55 .0125 .35])

set(Ax{1},'colormap',brewermap(25,'PuOr'))
set(Ax{2},'colormap',brewermap(25,'PuOr'))


Ax{5} = subplot('position',[.7875 .55 .1375 .35]);
pcolor(Bincent_H,Bincent_C, nseg_by_HC);
set(gca,'xlim',xlimmer_H);
set(gca,'ylim',xlimmer_C)
xlabel('IS2 Height','interpreter','latex')
% ylabel('CDR SIC','interpreter','latex')
set(gca,'yticklabel','')
set(gca,'clim',[0 2000])
shading flat
grid on; box on; 
title('Segment Number','interpreter','latex')

colormap(Ax{5},brewermap(20,'Spectral'))
colorbar('position',[.9375 .55 .0125 .35])



%%

clear *_by_*
param_use = (Hvals > 0) & (Hvals < 0.05) & (Nsegvals <= 1000);
get_param_plot_vectors; 


%%
Ax{6} = subplot('position',[.075 .1 .4 .325]);
cla
[p1,p2,p3,p4] = param_plot(Bincent_C,nC,bias_LIF_by_C,bias_LIF_by_C_up,bias_LIF_by_C_dn,bias_AMSR_by_C,bias_AMSR_by_C_up,bias_AMSR_by_C_dn, ...
    CDR_by_C,CDR_by_C_up,CDR_by_C_dn,xlimmer_C,ylimmer);

ylim([-.2 .4])
xlabel('CDR-SIC','interpreter','latex')
ylabel('$\Delta$ SIC','interpreter','latex')
title('Thin, low-photon ice: 0 $<$ H $<$ 5cm, N $\leq$ 1000','interpreter','latex')
%%

clear *_by_*
param_use = ~((Hvals > 0) & (Hvals < 0.05)) & Nsegvals >= 1000;

get_param_plot_vectors; 


%%

Ax{7} = subplot('position',[.55 .1 .4 .325]);
cla
[p1,p2,p3,p4] = param_plot(Bincent_C,nC,bias_LIF_by_C,bias_LIF_by_C_up,bias_LIF_by_C_dn,bias_AMSR_by_C,bias_AMSR_by_C_up,bias_AMSR_by_C_dn, ...
    CDR_by_C,CDR_by_C_up,CDR_by_C_dn,xlimmer_C,ylimmer);

xlabel('CDR-SIC','interpreter','latex')
title('Non-thin ice: H $\leq$ 0 or H $\geq$ 5cm, N $\geq$ 1000','interpreter','latex')
ylim([-.2 .4])

%%
allAxesInFigure = findall(gcf,'type','axes');
letter = {'(a)','(b)','(c)','(d)','(e)','(f)','(g)','(e)','(c)'};

for i = 1:length(Ax)
    
 posy = get(Ax{i},'position');

    set(Ax{i},'fontname','times','fontsize',8,'xminortick','on','yminortick','on')
    
    annotation('textbox',[posy(1) posy(2)+posy(4)-.015 .025 .025], ...
        'String',letter{i},'LineStyle','none','FontName','Helvetica', ...
        'FontSize',8,'Tag','legtag');
end


pos = [6.5 4]; 
set(gcf,'windowstyle','normal','position',[0 0 pos],'paperposition',[0 0 pos],'papersize',pos,'units','inches','paperunits','inches');
set(gcf,'windowstyle','normal','position',[0 0 pos],'paperposition',[0 0 pos],'papersize',pos,'units','inches','paperunits','inches');
print([OPTS.plot_save_str '2D_Parametric_Comp'],'-dpdf','-r600');


%%
% Ax{1} = subplot('position',[.075 .1 .225 .35]);
% 
% pcolor(Bincent_H,Bincent_C, bias_LIF_by_HC);
% set(gca,'xlim',xlimmer_H);
% set(gca,'ylim',xlimmer_C)
% xlabel('IS2 Height','interpreter','latex')
% ylabel('CDR SIC','interpreter','latex')
% set(gca,'clim',[-.5 .5])
% shading flat
% grid on; box on; 
% colorbar('position',[.325 .1 .0125 .35])
% title('LIF bias','interpreter','latex')


% 
% Ax{5} = subplot('position',[.075 .55 .225 .35]);
% 
% pcolor(Bincent_H,Bincent_C, spec_by_HC);
% set(gca,'xlim',xlimmer_H);
% set(gca,'ylim',xlimmer_C)
% xlabel('IS2 Height','interpreter','latex')
% ylabel('CDR SIC','interpreter','latex')
% set(gca,'clim',[0 .1])
% shading flat
% grid on; box on; 
% colorbar('position',[.325 .1 .0125 .35])
% title('Median Segment Count','interpreter','latex')
%%



%%

%%




%%

% allAxesInFigure = findall(gcf,'type','axes');
% letter = {'(c)','(b)','(a)','(d)','(e)','(f)','(g)','(e)','(c)'};
% 
% for i = 1:length(allAxesInFigure)
% 
%  posy = get(allAxesInFigure(i),'position');
% 
%     set(allAxesInFigure(i),'fontname','times','fontsize',8,'xminortick','on','yminortick','on')
% 
%     annotation('textbox',[posy(1) posy(2)+posy(4)-.015 .025 .025], ...
%         'String',letter{i},'LineStyle','none','FontName','Helvetica', ...
%         'FontSize',8,'Tag','legtag');
% 
% end
% 
% pos = [6.5 4.5]; 
% set(gcf,'windowstyle','normal','position',[0 0 pos],'paperposition',[0 0 pos],'papersize',pos,'units','inches','paperunits','inches');
% set(gcf,'windowstyle','normal','position',[0 0 pos],'paperposition',[0 0 pos],'papersize',pos,'units','inches','paperunits','inches');
% print([OPTS.plot_save_str '2D_Parametric_Comp'],'-dpdf','-r600');

%% Now look only at the overestimated values



%% try to remove the biased values

% param_use = (Hvals < 0.1 | Hvals > 0); 
% get_param_plot_vectors; 
% 
% ylimmer_bias = [-.15 .4];
% ylimmer_C = [0.15 1];
% 
% Ax{4} = subplot('position',[.1 .1 .35 .35]);
% 
% [p1,p2,p3,p4] = param_plot(Bincent_C,nC,bias_LIF_by_C,bias_LIF_by_C_up,bias_LIF_by_C_dn,bias_AMSR_by_C,bias_AMSR_by_C_up,bias_AMSR_by_C_dn, ...
%     CDR_by_C,CDR_by_C_up,CDR_by_C_dn,xlimmer_C,ylimmer);
% 
% hold on
% yyaxis left 
% hold on
% xline(.8,'linewidth',1,'label','MIZ')
% plot(Bincent_C(Bincent_C > 0.15),fitted(Bincent_C(Bincent_C > 0.15)),'--b','linewidth',1);
% plot(Bincent_C,1 - Bincent_C,'--r','linewidth',1)
% plot(Bincent_C,-Bincent_C,'--r','linewidth',1)
% xlabel('CDR SIC','interpreter','latex')
% ylabel('$\Delta$ from CDR','interpreter','latex')
% % legend([p1 p2 p3 p4],'AMSR offset','LIF offset','CDR range','Count','orientation','horizontal','location',[.3 .9 .4 .05])
% 
% %% 
% param_use = param_use_overbiased; 
% get_param_plot_vectors; 
% 
% 
% Ax{5} = subplot('position',[.55 .1 .35 .35]);
% 
% [p1,p2,p3,p4] = param_plot(Bincent_C,nC,bias_LIF_by_C,bias_LIF_by_C_up,bias_LIF_by_C_dn,bias_AMSR_by_C,bias_AMSR_by_C_up,bias_AMSR_by_C_dn, ...
%     CDR_by_C,CDR_by_C_up,CDR_by_C_dn,xlimmer_C,ylimmer);
% 
% hold on
% yyaxis left 
% hold on
% xline(.8,'linewidth',1,'label','MIZ')
% plot(Bincent_C(Bincent_C > 0.15),fitted(Bincent_C(Bincent_C > 0.15)),'--b','linewidth',1);
% plot(Bincent_C,1 - Bincent_C,'--r','linewidth',1)
% plot(Bincent_C,-Bincent_C,'--r','linewidth',1)
% xlabel('CDR SIC','interpreter','latex')
% ylabel('$\Delta$ from CDR','interpreter','latex')
% % 
% 
% 
% %% For AMSR-only figure
% 
% 
% Ax{4} = subplot('position',[.1 .1 .275 .35]);
% 
% pcolor(Bincent_W,Bincent_C, bias_AMSR_by_WC);
% set(gca,'xlim',xlimmer_W);
% set(gca,'ylim',xlimmer_C)
% xlabel('IS2 WAF','interpreter','latex')
% ylabel('CDR SIC','interpreter','latex')
% colormap(brewermap(25,'PuOr'))
% set(gca,'clim',[-.3 .3])
% shading flat
% title('LIF-CDR Offsets','interpreter','latex')
% grid on; box on; 
% 
% 
% Ax{5} = subplot('position',[.4 .1 .275 .35]);
% 
% pcolor(Bincent_W,Bincent_H, bias_AMSR_by_WH);
% set(gca,'xlim',xlimmer_W);
% set(gca,'ylim',xlimmer_H)
% xlabel('IS2 WAF','interpreter','latex')
% ylabel('IS2 Height','interpreter','latex')
% colormap(brewermap(25,'PuOr'))
% set(gca,'clim',[-.3 .3])
% shading flat
% grid on; box on; 
% % title('LIF-CDR Offsets','interpreter','latex')
% 
% Ax{6} = subplot('position',[.7 .1 .275 .35]);
% 
% pcolor(Bincent_H,Bincent_C, bias_AMSR_by_HC);
% set(gca,'xlim',xlimmer_H);
% set(gca,'ylim',xlimmer_C)
% xlabel('IS2 Height','interpreter','latex')
% ylabel('CDR SIC','interpreter','latex')
% colormap(brewermap(25,'PuOr'))
% set(gca,'clim',[-.3 .3])
% shading flat
% grid on; box on; 
% % title('LIF-CDR Offsets','interpreter','latex')


