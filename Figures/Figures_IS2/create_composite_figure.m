close all

[binct,~,binval] = histcounts(Dvals,Dbins);

disp('------')
fprintf('Fraction of tracks that are segmented here is %2.2f \n',100*nusable_segmented./nusable);
fprintf('Fraction of the tracks that is MIZ is %2.2f \n',100*sum(Dvals < 0)./numel(Dvals));
fprintf('Fraction of the tracks that is MIZ or wavy is %2.2f \n',100*sum(Dvals < 0 | WAFvals > wave_thresh)./numel(Dvals));

fprintf('Fraction of those that have bias is %2.2f \n',100*sum(WAFvals > wave_thresh)./numel(WAFvals));


% binval(binval == 0) = length(Dbins);

upval = @(x) prctile(x,75);
dnval = @(x) prctile(x,25);

Nvec = accumarray(binval,1,[length(Dbins)-1 1],@sum);
%%

Nsegvec = accumarray(binval,Nsegvals,[length(Dbins)-1 1],@mean);
Nsegvar = accumarray(binval,Nsegvals,[length(Dbins)-1 1],@std);

SICvec = accumarray(binval,SICvals,[length(Dbins)-1 1],@median);

SICup = accumarray(binval,SICvals,[length(Dbins) - 1 1],upval);
SICup(isinf(SICup)) = 1;
SICdn = accumarray(binval,SICvals,[length(Dbins) - 1 1],dnval);

if IS2_DATA.v6

    SICvec_amsr = accumarray(binval,SICvals_amsr,[length(Dbins) - 1 1],@median);
    SICup_amsr= accumarray(binval,SICvals_amsr,[length(Dbins) - 1 1],upval);
    SICdn_amsr = accumarray(binval,SICvals_amsr,[length(Dbins) - 1 1],dnval);

    biasvec = accumarray(binval,biasvals,[length(Dbins) - 1 1],@median);
    biasup= accumarray(binval,biasvals,[length(Dbins) - 1 1],upval);
    biasdn = accumarray(binval,biasvals,[length(Dbins) - 1 1],dnval);


end

bias_LIFvec = accumarray(binval,biasvals_LIF,[length(Dbins) - 1 1],@median);
bias_LIFup= accumarray(binval,biasvals_LIF,[length(Dbins) - 1 1],upval);
bias_LIFdn = accumarray(binval,biasvals_LIF,[length(Dbins) - 1 1],dnval);


WAFvec = accumarray(binval,WAFvals,[length(Dbins) - 1 1],@median);
WAFup = accumarray(binval,WAFvals,[length(Dbins) - 1 1],upval);
WAFdn = accumarray(binval,WAFvals,[length(Dbins) - 1 1],dnval);

Hvec = accumarray(binval,Hvals,[length(Dbins) - 1 1],@median);
Hup = accumarray(binval,Hvals,[length(Dbins) - 1 1],upval);
Hdn = accumarray(binval,Hvals,[length(Dbins) - 1 1],dnval);

LIFvec = accumarray(binval,LIFvals,[length(Dbins) - 1 1],@median);
LIFup = accumarray(binval,LIFvals,[length(Dbins) - 1 1],upval);
LIFdn = accumarray(binval,LIFvals,[length(Dbins) - 1 1],dnval);

LIF_spec_vec = accumarray(binval,LIF_spec_vals,[length(Dbins) - 1 1],@median);
LIF_spec_up = accumarray(binval,LIF_spec_vals,[length(Dbins) - 1 1],upval);
LIF_spec_dn = accumarray(binval,LIF_spec_vals,[length(Dbins) - 1 1],dnval);

LIF_dark_vec = accumarray(binval,LIF_dark_vals,[length(Dbins) - 1 1],@median);
LIF_dark_up = accumarray(binval,LIF_dark_vals,[length(Dbins) - 1 1],upval);
LIF_dark_dn = accumarray(binval,LIF_dark_vals,[length(Dbins) - 1 1],dnval);



Hvar = accumarray(binval,Hvals,[length(Dbins) - 1 1],@std);

Evec = accumarray(binval,Evals,[length(Dbins) - 1 1],@mean);
Evar = accumarray(binval,Evals,[length(Dbins) - 1 1],@std);

%
xfirst = find(Nvec > 1,1,'first');

dum = find(Nvec > 1,2,'last');

xlast = dum(1);

figure;

subplot('position',[.1 .55 .8 .4])


jbfill(Bincent_D,SICup',SICdn',[.4 .4 .4],[1 1 1],1,.3);
hold on

if IS2_DATA.v6

    jbfill(Bincent_D,SICup_amsr',SICdn_amsr',[.8 .2 .2],[1 1 1],1,.3);
    jbfill(Bincent_D,biasup',biasdn',[.8 .8 .8],[1 1 1],1,.3);

end


hold on
p1 = plot(Bincent_D,SICvec,'k','linewidth',2);

if IS2_DATA.v6

    p3 = plot(Bincent_D,SICvec_amsr,'color',[.8 .2 .2],'linewidth',2);
    p6 = plot(Bincent_D,biasvec,'--','color',[.8 .8 .8],'linewidth',2);

end

p7 = plot(Bincent_D,bias_LIFvec,'--','color',[.2 .8 .8],'linewidth',2);

jbfill(Bincent_D,bias_LIFup',bias_LIFdn',[.2 .8 .8],[1 1 1],1,.2);
hold on
jbfill(Bincent_D,LIFup',LIFdn',[.2 .2 .8],[1 1 1],1,.4);

hold on

p2 = plot(Bincent_D,LIFvec,'color',[.4 .4 .8],'linewidth',2);

% p5 = plot(Bincent_D,LIF_spec_vec,'--','color',[.8 .4 .2],'linewidth',2);
% p6 = plot(Bincent_D,LIF_dark_vec,'--','color',[.2 .4 .8],'linewidth',2);


xlimmer = [-500 500];

% xlimmer = [Bincent_D(xfirst) Bincent_D(xlast)];
% xlimmer = min(abs(xlimmer))*[-1 1];

plot(Bincent_D,LIFup,'--b');
plot(Bincent_D,LIFdn,'--b');
xlim(xlimmer)
grid on; box on;
ylim([-.1 1])
yline(0,'-k');
ylabel('Ice Fraction','Interpreter','latex')
xline(0,'color',[.2 .2 .2],'linewidth',1)
set(gca,'xticklabel','')
% fitted = @(c) -1.639*(c-.6122).^2 + 0.2316; % This is weighted fit
fitted = @(c) -1.604*(c-.6138).^2 + 0.229; % This is naive fit.

p5 = plot(Bincent_D,fitted(SICvec),'--b','linewidth',1);


yyaxis right
p4 = plot(Bincent_D,Nvec,'color',[.8 .4 .8],'LineWidth',2);
set(gca,'ycolor','k')
ylabel('Number of Samples','interpreter','latex')
xlim(xlimmer)


xline(-50,'--','color',[.8 .4 .4])
xline(50,'--','color',[.8 .4 .4])


if IS2_DATA.v6

    legend([p1 p2 p3 p4 p5 p6 p7],{'CDR','LIF','AMSR2-NT2','Stencil #','Fit','NT2-CDR','LIF-CDR'},'location','best')

else

    legend([p1 p2 p4 p5 p7],{'CDR','LIF','Number','Spec','LIF Bias'},'location','best')

end



%% 
subplot('position',[.1 .325 .8 .2])

% jbfill(Bincent_D,Hup',Hdn',[.4 .4 .4],[1 1 1],1,.3)

jbfill(Bincent_D,Hup',Hdn',[.4 .4 .4],[1 1 1],1,.3);
hold on
r1 = plot(Bincent_D,Hvec,'k','linewidth',2);

plot(Bincent_D,Hup,'--k');
plot(Bincent_D,Hdn,'--k');
xlim(xlimmer)
% ylim([0 .5]);
set(gca,'ylim',[0 max(max(get(gca,'ylim')),0.5)])
ylabel('Ice Height','interpreter','latex')
% r2 = yline(.075,'--','color',[.8 .4 .4],'linewidth',1)
xline(0,'color',[.2 .2 .2],'linewidth',1)
grid on; box on;
% xlabel('Kilometers Relative to PM-MIZ','interpreter','latex')

xline(0,'label','CDR-defined MIZ','interpreter','latex','fontsize',8,'LabelOrientation','horizontal','LabelHorizontalAlignment','left');
xline(0,'label','CDR-defined CIZ','interpreter','latex','fontsize',8,'LabelOrientation','horizontal');
xline(-50,'--','color',[.8 .4 .4])
xline(50,'--','color',[.8 .4 .4])
set(gca,'xticklabel','')

%%

subplot('position',[.1 .1 .8 .2])

% jbfill(Bincent_D,Hup',Hdn',[.4 .4 .4],[1 1 1],1,.3)

jbfill(Bincent_D,WAFup',WAFdn',[.4 .4 .4],[1 1 1],1,.3);
hold on
r1 = plot(Bincent_D,WAFvec,'k','linewidth',2);

plot(Bincent_D,WAFup,'--k');
plot(Bincent_D,WAFdn,'--k');
xlim(xlimmer)
% ylim([0 .5]);
set(gca,'ylim',[0 max(max(get(gca,'ylim')),0.5)])
ylabel('WAF','interpreter','latex')

yline(.075,'label','Horvat et al (2020) WMIZ','color',[.4 .4 .8],'linewidth',1,'interpreter','latex','fontsize',8,'LabelOrientation','horizontal');

xline(0,'color',[.2 .2 .2],'linewidth',1)
grid on; box on;
xlabel('Kilometers Relative to PM-MIZ','interpreter','latex')

xline(0,'label','CDR-defined MIZ','interpreter','latex','fontsize',8,'LabelOrientation','horizontal','LabelHorizontalAlignment','left');
xline(0,'label','CDR-defined CIZ','interpreter','latex','fontsize',8,'LabelOrientation','horizontal');
xline(-50,'--','color',[.8 .4 .4])
xline(50,'--','color',[.8 .4 .4])

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

pos = [6.5 4]; 
set(gcf,'windowstyle','normal','position',[0 0 pos],'paperposition',[0 0 pos],'papersize',pos,'units','inches','paperunits','inches');
set(gcf,'windowstyle','normal','position',[0 0 pos],'paperposition',[0 0 pos],'papersize',pos,'units','inches','paperunits','inches');
print([OPTS.plot_save_str 'MIZ-SIC-comp-all'],'-dpdf','-r600');
