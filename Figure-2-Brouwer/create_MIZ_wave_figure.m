
nusable = sum(usable_all);
nusable_segmented = sum(usable); 


if IS2_DATA.v6

usable = usable & ~isnan(SICvals_amsr); 

SICvals_amsr = SICvals_amsr(usable); 
biasvals = biasvals(usable);

end


Nvals = Nvals(usable); 
SICvals = SICvals(usable); 
LIFvals = LIFvals(usable); 
LIF_spec_vals = LIF_spec_vals(usable); 
LIF_dark_vals = LIF_dark_vals(usable); 
Dvals = Dvals(usable); 
Hvals = Hvals(usable); 
Evals = Evals(usable); 
WAFvals = WAFvals(usable);
wavytracks = wavytracks(usable); 

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


Nsegvec = accumarray(binval,Nvals,[length(Dbins)-1 1],@mean); 
Nsegvar = accumarray(binval,Nvals,[length(Dbins)-1 1],@std); 

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

subplot('position',[.1 .6 .8 .35])

jbfill(Bincent,SICup',SICdn',[.4 .4 .4],[1 1 1],1,.3);
hold on

if IS2_DATA.v6

jbfill(Bincent,SICup_amsr',SICdn_amsr',[.8 .2 .2],[1 1 1],1,.3);
jbfill(Bincent,biasup',biasdn',[.8 .8 .8],[1 1 1],1,.3);

end


jbfill(Bincent,LIFup',LIFdn',[.2 .2 .8],[1 1 1],1,.3);
hold on
p1 = plot(Bincent,SICvec,'k','linewidth',2); 

if IS2_DATA.v6

p3 = plot(Bincent,SICvec_amsr,'color',[.8 .2 .2],'linewidth',2); 
p6 = plot(Bincent,biasvec,'color',[.8 .8 .8],'linewidth',2); 

end

p2 = plot(Bincent,LIFvec,'color',[.4 .4 .8],'linewidth',2); 

p5 = plot(Bincent,LIF_spec_vec,'--','color',[.8 .4 .2],'linewidth',2); 
% p6 = plot(Bincent,LIF_dark_vec,'--','color',[.2 .4 .8],'linewidth',2); 


xlimmer = [-500 500];
 
% xlimmer = [Bincent(xfirst) Bincent(xlast)]; 
% xlimmer = min(abs(xlimmer))*[-1 1];

plot(Bincent,LIFup,'--b'); 
plot(Bincent,LIFdn,'--b'); 
xlim(xlimmer)
grid on; box on; 
ylim([0 1])
ylabel('Ice Fraction','Interpreter','latex')
xline(0,'color',[.2 .2 .2],'linewidth',1)

yyaxis right
p4 = plot(Bincent,Nvec,'color',[.8 .4 .8],'LineWidth',2);
set(gca,'ycolor','k')
ylabel('Number of Samples','interpreter','latex')
xlim(xlimmer)

if IS2_DATA.v6

legend([p1 p2 p3 p4 p5 p6],{'CDR','LIF','AMSR2','Number','Spec','Bias'},'location','best')

else

    legend([p1 p2 p4 p5],{'CDR','LIF','Number','Spec'},'location','best')

end

%%
% subplot('position',[.1 .1 .8 .25])
% 
% jbfill(Bincent,WAFup',WAFdn',[.4 .4 .4],[1 1 1],1,.3)
% hold on
% r1 = plot(Bincent,WAFvec,'k','linewidth',2); 
% 
% plot(Bincent,WAFup,'--k'); 
% plot(Bincent,WAFdn,'--k'); 
% xlim(xlimmer)
% ylim([0 .5]); 
% ylabel('Wave Affected Fraction','interpreter','latex')
% r2 = yline(.075,'--','color',[.8 .4 .4],'linewidth',1)
% xline(0,'color',[.2 .2 .2],'linewidth',1)
% grid on; box on;
% xlabel('Kilometers Relative to PM-MIZ','interpreter','latex')
% legend([r1 r2],{'Wave Affected Fraction','WMIZ Threshold'})

%% 
subplot('position',[.1 .35 .8 .2])


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
% r2 = yline(.075,'--','color',[.8 .4 .4],'linewidth',1)
xline(0,'color',[.2 .2 .2],'linewidth',1)
grid on; box on;
xlabel('Kilometers Relative to PM-MIZ','interpreter','latex')
% legend([r1 r2],{'Wave Affected Fraction','WMIZ Threshold'})
%% SECTION TITLE
subplot('position',[.1 .1 .8 .2])


% jbfill(Bincent,Hup',Hdn',[.4 .4 .4],[1 1 1],1,.3)

jbfill(Bincent,Hup',Hdn',[.4 .4 .4],[1 1 1],1,.3);
hold on
r1 = plot(Bincent,Hvec,'k','linewidth',2); 

plot(Bincent,Hup,'--k'); 
plot(Bincent,Hdn,'--k'); 
xlim(xlimmer)
% ylim([0 .5]); 
set(gca,'ylim',[0 max(max(get(gca,'ylim')),0.5)])
ylabel('Freeboard Height','interpreter','latex')
% r2 = yline(.075,'--','color',[.8 .4 .4],'linewidth',1)
xline(0,'color',[.2 .2 .2],'linewidth',1)
grid on; box on;
xlabel('Kilometers Relative to PM-MIZ','interpreter','latex')
%%
allAxesInFigure = findall(gcf,'type','axes');
letter = {'(b)','(a)','(c)','(d)','(e)','(f)','(g)','(e)','(c)'};

for i = 1:length(allAxesInFigure)
    
 posy = get(allAxesInFigure(i),'position');

    set(allAxesInFigure(i),'fontname','times','fontsize',8,'xminortick','on','yminortick','on')
    
    annotation('textbox',[posy(1) - .025 posy(2)+posy(4) + .035 .025 .025], ...
        'String',letter{i},'LineStyle','none','FontName','Helvetica', ...
        'FontSize',8,'Tag','legtag');

end