
close all

% Want to plot the parametric mapping between things. 



param_use = SICvals <= 1 & Hvals >= 0; 

nbins = 100*floor(round(sqrt(sum(param_use)))/250) + 1;



Hbins = linspace(-6,6,nbins);
Bincent_H = 0.5*(Hbins(1:end-1) + Hbins(2:end));

Wbins = linspace(0,2,nbins);
Wbins = Wbins - 0.5*Wbins(2); 
Bincent_W = 0.5*(Wbins(1:end-1) + Wbins(2:end));

Cbins = linspace(0,1,nbins);
Bincent_C = 0.5*(Cbins(1:end-1) + Cbins(2:end));


[nH,~,mapper_H] = histcounts(Hvals(param_use),Hbins);
[nW,~,mapper_W] = histcounts(WAFvals(param_use),Wbins);
[nC,~,mapper_C] = histcounts(SICvals(param_use),Cbins);

mapper_C(SICvals(param_use) > 1) = length(Cbins)-1;

ncutoff = nbins;


xlimmer_C = [Cbins(find(nC>ncutoff,1)) 1];
xlimmer_H = [Hbins(find(nH>ncutoff,1)) Hbins(find(nH>ncutoff,1,'last'))];
xlimmer_W = [Wbins(find(nW>ncutoff,1)) Wbins(find(nW>ncutoff,1,'last'))];


ylimmer = [-.15 .4];
%%

LIF_by_C = accumarray(mapper_C,LIFvals(param_use),[length(Cbins)-1 1],@nanmedian);
LIF_by_C_up = accumarray(mapper_C,LIFvals(param_use),[length(Cbins)-1 1],upval);
LIF_by_C_dn = accumarray(mapper_C,LIFvals(param_use),[length(Cbins)-1 1],dnval);


SIC_by_C = accumarray(mapper_C,SICvals(param_use),[length(Cbins)-1 1],@nanmedian);
SIC_by_C_up = accumarray(mapper_C,SICvals(param_use),[length(Cbins)-1 1],upval);
SIC_by_C_dn = accumarray(mapper_C,SICvals(param_use),[length(Cbins)-1 1],dnval);


AMSR_by_C = accumarray(mapper_C,SICvals_amsr(param_use),[length(Cbins)-1 1],@nanmedian);
AMSR_by_C_up = accumarray(mapper_C,SICvals_amsr(param_use),[length(Cbins)-1 1],upval);
AMSR_by_C_dn = accumarray(mapper_C,SICvals_amsr(param_use),[length(Cbins)-1 1],dnval);


bias_by_C = accumarray(mapper_C,biasvals(param_use),[length(Cbins)-1 1],@nanmedian);
bias_by_C_up = accumarray(mapper_C,biasvals(param_use),[length(Cbins)-1 1],upval);
bias_by_C_dn = accumarray(mapper_C,biasvals(param_use),[length(Cbins)-1 1],dnval);

bias_LIF_by_C = accumarray(mapper_C,biasvals_LIF(param_use),[length(Cbins)-1 1],@nanmedian);
bias_LIF_by_C_up = accumarray(mapper_C,biasvals_LIF(param_use),[length(Cbins)-1 1],upval);
bias_LIF_by_C_dn = accumarray(mapper_C,biasvals_LIF(param_use),[length(Cbins)-1 1],dnval);


%%

LIF_by_H = accumarray(mapper_H,LIFvals(param_use),[length(Hbins)-1 1],@nanmedian);
LIF_by_H_up = accumarray(mapper_H,LIFvals(param_use),[length(Hbins)-1 1],upval);
LIF_by_H_dn = accumarray(mapper_H,LIFvals(param_use),[length(Hbins)-1 1],dnval);


SIC_by_H = accumarray(mapper_H,SICvals(param_use),[length(Hbins)-1 1],@nanmedian);
SIC_by_H_up = accumarray(mapper_H,SICvals(param_use),[length(Hbins)-1 1],upval);
SIC_by_H_dn = accumarray(mapper_H,SICvals(param_use),[length(Hbins)-1 1],dnval);


AMSR_by_H = accumarray(mapper_H,SICvals_amsr(param_use),[length(Hbins)-1 1],@nanmedian);
AMSR_by_H_up = accumarray(mapper_H,SICvals_amsr(param_use),[length(Hbins)-1 1],upval);
AMSR_by_H_dn = accumarray(mapper_H,SICvals_amsr(param_use),[length(Hbins)-1 1],dnval);


bias_by_H = accumarray(mapper_H,biasvals(param_use),[length(Hbins)-1 1],@nanmedian);
bias_by_H_up = accumarray(mapper_H,biasvals(param_use),[length(Hbins)-1 1],upval);
bias_by_H_dn = accumarray(mapper_H,biasvals(param_use),[length(Hbins)-1 1],dnval);

bias_LIF_by_H = accumarray(mapper_H,biasvals_LIF(param_use),[length(Hbins)-1 1],@nanmedian);
bias_LIF_by_H_up = accumarray(mapper_H,biasvals_LIF(param_use),[length(Hbins)-1 1],upval);
bias_LIF_by_H_dn = accumarray(mapper_H,biasvals_LIF(param_use),[length(Hbins)-1 1],dnval);

%%
LIF_by_W = accumarray(mapper_W,LIFvals(param_use),[length(Wbins)-1 1],@nanmedian);
LIF_by_W_up = accumarray(mapper_W,LIFvals(param_use),[length(Wbins)-1 1],upval);
LIF_by_W_dn = accumarray(mapper_W,LIFvals(param_use),[length(Wbins)-1 1],dnval);


SIC_by_W = accumarray(mapper_W,SICvals(param_use),[length(Wbins)-1 1],@nanmedian);
SIC_by_W_up = accumarray(mapper_W,SICvals(param_use),[length(Wbins)-1 1],upval);
SIC_by_W_dn = accumarray(mapper_W,SICvals(param_use),[length(Wbins)-1 1],dnval);


AMSR_by_W = accumarray(mapper_W,SICvals_amsr(param_use),[length(Wbins)-1 1],@nanmedian);
AMSR_by_W_up = accumarray(mapper_W,SICvals_amsr(param_use),[length(Wbins)-1 1],upval);
AMSR_by_W_dn = accumarray(mapper_W,SICvals_amsr(param_use),[length(Wbins)-1 1],dnval);


bias_by_W = accumarray(mapper_W,biasvals(param_use),[length(Wbins)-1 1],@nanmedian);
bias_by_W_up = accumarray(mapper_W,biasvals(param_use),[length(Wbins)-1 1],upval);
bias_by_W_dn = accumarray(mapper_W,biasvals(param_use),[length(Wbins)-1 1],dnval);

bias_LIF_by_W = accumarray(mapper_W,biasvals_LIF(param_use),[length(Wbins)-1 1],@nanmedian);
bias_LIF_by_W_up = accumarray(mapper_W,biasvals_LIF(param_use),[length(Wbins)-1 1],upval);
bias_LIF_by_W_dn = accumarray(mapper_W,biasvals_LIF(param_use),[length(Wbins)-1 1],dnval);

%%

subplot(1,3,1)


jbfill(Bincent_C,bias_LIF_by_C_up',bias_LIF_by_C_dn',[.4 .4 .8],[1 1 1],1,.2);
hold on
jbfill(Bincent_C,bias_by_C_up',bias_by_C_dn',[.8 .2 .2],[1 1 1],1,.3);
hold on
plot(Bincent_C,bias_by_C,'color',[.8 .2 .2],'linewidth',2)
hold on
plot(Bincent_C,bias_LIF_by_C,'color',[.4 .4 .8],'linewidth',2)
xline(.8,'linewidth',1,'label','MIZ')
% yline(.8,'linewidth',1,'label','MIZ')
% line([0 1],[0 1],'Color',[.3 .3 .3],'linestyle','--')
% title('SIC Comparison','interpreter','latex');

xlim(xlimmer_C);
ylim(ylimmer)
% ylim([.15 1]);

plot(Bincent_C(Bincent_C > 0.15),fitted(Bincent_C(Bincent_C > 0.15)),'--b','linewidth',1);
plot(Bincent_C,1 - Bincent_C,'--r','linewidth',1)
plot(Bincent_C,-Bincent_C,'--r','linewidth',1)
yline(0,'k','linewidth',1')
ylabel('$\Delta$ from CDR','interpreter','latex')

grid on; box on;

yyaxis right
set(gca,'ycolor','k','yticklabel','')
bar(Bincent_C,nC/sum(nC),1,'FaceColor',[.8 .8 .8],'FaceAlpha',.5,'EdgeColor','none');
xlim(xlimmer_C);

xlabel('CDR SIC')




%%

subplot(1,3,2)


% jbfill(Bincent_H,LIF_by_H_up',LIF_by_H_dn',[.4 .4 .8],[1 1 1],1,.2);
% hold on
% jbfill(Bincent_H,SIC_by_H_up',SIC_by_H_dn',[.4 .4 .4],[1 1 1],1,.3);
% hold on
% jbfill(Bincent_H,AMSR_by_H_up',AMSR_by_H_dn',[.8 .2 .2],[1 1 1],1,.3);
% hold on
% plot(Bincent_H,SIC_by_H,'-k','linewidth',2)
% plot(Bincent_H,LIF_by_H,'color',[.4 .4 .8],'linewidth',2)
% 
% 
% plot(Bincent_H,AMSR_by_H,'color',[.8 .2 .2],'linewidth',2)
% % xline(.8,'linewidth',1,'label','MIZ')
% % yline(.8,'linewidth',1,'label','MIZ')
% % line([0 1],[0 1],'Color',[.3 .3 .3],'linestyle','--')
% % title('SIC Comparison','interpreter','latex');
% 
% xlim(xlimmer_H);
% % ylim([.15 1]);
% 
% % ylabel('AMSR2-NT2 SIC','interpreter','latex');
% % xlabel('AMSR2-BS SIC','interpreter','latex')
% grid on; box on;
% 
% yyaxis right
% set(gca,'ycolor','k','yticklabel','')
% bar(Bincent_H,nH/sum(nH),1,'FaceColor',[.4 .4 .8],'FaceAlpha',.5,'EdgeColor','none');
% xlim(xlimmer_H);
% ylim([0 1])
% 
% 
% subplot(2,2,2)

%jbfill(Bincent_C,plot_up',plot_dn',[.4 .4 .4],[1 1 1],1,.3);
% hold on

jbfill(Bincent_H,bias_LIF_by_H_up',bias_LIF_by_H_dn',[.4 .4 .8],[1 1 1],1,.2);
hold on
jbfill(Bincent_H,bias_by_H_up',bias_by_H_dn',[.8 .2 .2],[1 1 1],1,.3);
hold on
plot(Bincent_H,bias_by_H,'color',[.8 .2 .2],'linewidth',2)
hold on
plot(Bincent_H,bias_LIF_by_H,'color',[.4 .4 .8],'linewidth',2)
yline(0,'k','linewidth',1')
% line([0 1],[0 1],'Color',[.3 .3 .3],'linestyle','--')
% title('SIC Comparison','interpreter','latex');

xlim(xlimmer_H);
ylim(ylimmer)
grid on; box on;

yyaxis right
set(gca,'ycolor','k','yticklabel','')
bar(Bincent_H,nH/sum(nH),1,'FaceColor',[.8 .8 .8],'FaceAlpha',.5,'EdgeColor','none');
xlim(xlimmer_H);
xlabel('IS2 Surface Height')

%%


subplot(1,3,3)

%jbfill(Bincent_C,plot_up',plot_dn',[.4 .4 .4],[1 1 1],1,.3);
% hold on

jbfill(Bincent_W,bias_LIF_by_W_up',bias_LIF_by_W_dn',[.4 .4 .8],[1 1 1],1,.2);
hold on
jbfill(Bincent_W,bias_by_W_up',bias_by_W_dn',[.8 .2 .2],[1 1 1],1,.3);
hold on
plot(Bincent_W,bias_by_W,'color',[.8 .2 .2],'linewidth',2)
hold on
plot(Bincent_W,bias_LIF_by_W,'color',[.4 .4 .8],'linewidth',2)
xline(.075,'linewidth',1,'label','WMIZ')
% yline(.8,'linewidth',1,'label','MIZ')
% line([0 1],[0 1],'Color',[.3 .3 .3],'linestyle','--')
% title('SIC Comparison','interpreter','latex');

xlim(xlimmer_W);
ylim(ylimmer)
yline(0,'k','linewidth',1')

% ylim([.15 1]);

% ylabel('AMSR2-NT2 SIC','interpreter','latex');
% xlabel('AMSR2-BS SIC','interpreter','latex')
grid on; box on;

yyaxis right
set(gca,'ycolor','k','yticklabel','')
bar(Bincent_W,nW/sum(nW),1,'FaceColor',[.8 .8 .8],'FaceAlpha',.5,'EdgeColor','none');
xlim(xlimmer_W);
xlabel('IS2 WAF')
ylim([0 .05])

%%



allAxesInFigure = findall(gcf,'type','axes');
letter = {'(c)','(b)','(a)','(d)','(e)','(f)','(g)','(e)','(c)'};

for i = 1:length(allAxesInFigure)
    
 posy = get(allAxesInFigure(i),'position');

    set(allAxesInFigure(i),'fontname','times','fontsize',8,'xminortick','on','yminortick','on')
    
    annotation('textbox',[posy(1) - .05 posy(2)+posy(4)-.015 .025 .025], ...
        'String',letter{i},'LineStyle','none','FontName','Helvetica', ...
        'FontSize',8,'Tag','legtag');

end

pos = [6.5 3]; 
set(gcf,'windowstyle','normal','position',[0 0 pos],'paperposition',[0 0 pos],'papersize',pos,'units','inches','paperunits','inches');
set(gcf,'windowstyle','normal','position',[0 0 pos],'paperposition',[0 0 pos],'papersize',pos,'units','inches','paperunits','inches');
print([OPTS.plot_save_str 'Parametric_Comp'],'-dpdf','-r600');