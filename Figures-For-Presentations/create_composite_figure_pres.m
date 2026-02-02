close all

figure; 
pos = [4.5 3]; 
set(gcf,'windowstyle','normal','position',[0 0 pos],'paperposition',[0 0 pos],'papersize',pos,'units','inches','paperunits','inches');
set(gcf,'windowstyle','normal','position',[0 0 pos],'paperposition',[0 0 pos],'papersize',pos,'units','inches','paperunits','inches');

set(gca,'fontname','times','fontsize',8,'xminortick','on','yminortick','on')

subplot('position',[.125 .15 .775 .775])



hold on
p2 = plot(Bincent,LIFvec,'color',[.4 .4 .8],'linewidth',2); 
plot(Bincent,LIFup,'--b'); 
plot(Bincent,LIFdn,'--b'); 
jbfill(Bincent,LIFup',LIFdn',[.2 .2 .8],[1 1 1],1,.4);


xlim(xlimmer)
grid on; box on; 
ylim([-.1 1])
yline(0,'-k'); 
ylabel('Ice Fraction','Interpreter','latex')
xline(0,'color',[.2 .2 .2],'linewidth',1)
xlabel('Kilometers Relative to PM-MIZ','interpreter','latex')


yyaxis right
p4 = plot(Bincent,Nvec,'color',[.8 .4 .8],'LineWidth',2);
set(gca,'ycolor','k')
ylabel('Number of Samples','interpreter','latex')
xlim(xlimmer)

saveas(gcf,'LIF-only.pdf'); 
saveas(gcf,'LIF-only.jpeg'); 

%%

yyaxis left; 
hold on 

p1 = plot(Bincent,SICvec,'k','linewidth',2); 

jbfill(Bincent,SICup',SICdn',[.4 .4 .4],[1 1 1],1,.3);
hold on

if IS2_DATA.v6

% jbfill(Bincent,SICup_amsr',SICdn_amsr',[.8 .2 .2],[1 1 1],1,.3);
jbfill(Bincent,bias_LIFup',bias_LIFdn',[.2 .8 .8],[1 1 1],1,.2);

end

hold on

p7 = plot(Bincent,bias_LIFvec,'--','color',[.2 .8 .8],'linewidth',2); 

% legend([p1 p2 p4 p7],{'CDR','LIF','Number','Bias'})

saveas(gcf,'LIF-to-CDR.pdf'); 
saveas(gcf,'LIF-to-CDR.jpeg'); 

%%

if IS2_DATA.v6

p3 = plot(Bincent,SICvec_amsr,'color',[.8 .2 .2],'linewidth',2); 
p6 = plot(Bincent,biasvec,'--','color',[.8 .2 .2],'linewidth',2); 

end

% legend([p1 p2 p4 p7 p3 p6],{'CDR','LIF','Number','Bias','NT2','NT2-bias'})


saveas(gcf,'LIF-to-CDR-AMSR.pdf'); 
saveas(gcf,'LIF-to-CDR-AMSR.jpeg'); 


%%


% jbfill(Bincent,bias_LIFup',bias_LIFdn',[.2 .8 .8],[1 1 1],1,.2);
% hold on
% 


% p5 = plot(Bincent,LIF_spec_vec,'--','color',[.8 .4 .2],'linewidth',2); 
% p6 = plot(Bincent,LIF_dark_vec,'--','color',[.2 .4 .8],'linewidth',2); 


% xlimmer = [-500 500];
% 
% % xlimmer = [Bincent(xfirst) Bincent(xlast)]; 
% % xlimmer = min(abs(xlimmer))*[-1 1];
% 
% 
% %  set(gca,'xticklabel','')
% % fitted = @(c) -1.639*(c-.6122).^2 + 0.2316; % This is weighted fit
% fitted = @(c) -1.604*(c-.6138).^2 + 0.229; % This is naive fit. 
% 
% p5 = plot(Bincent,fitted(SICvec),'--b','linewidth',1);
% 
% 
% 
% if IS2_DATA.v6
% 
% legend([p1 p2 p3 p4 p5 p6 p7],{'CDR','LIF','AMSR2','Number','Fit','AMSR2 Bias','LIF Bias'},'location','best')
% 
% else
% 
%     legend([p1 p2 p4 p5 p7],{'CDR','LIF','Number','Spec','LIF Bias'},'location','best')
% 
% end

%% SECTION TITLE