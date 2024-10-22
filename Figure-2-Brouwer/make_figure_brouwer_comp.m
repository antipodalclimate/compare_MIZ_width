function make_figure_brouwer_comp(OPTS)
% Load in the segmented statistics. Each is an array of stats which is
% indexed by
% nT - number of tracks
% nB - number of beams (this is frequently double the actual number of
% beams because we count Northward and Southward halves of a single track
% the same way

load(OPTS.output_str);

nT = size(MIZ_DATA.timer,1);
nB = size(MIZ_DATA.timer,2);

%%

for i = 1:nT
    for j = 1:2 % Both forward and reverse tracks, all strong beams
        
        beaminds = (nB/2)*floor(j/2)+1:(nB/2)*floor(j/2) + nB/2; 
        
        if ~isempty(MIZ_DATA.SIC{i,j})
            
            % Distance from the edge
            var1 = vertcat(MIZ_DATA.D_to_edge{i,beaminds});
            var2 = vertcat(MIZ_DATA.SIC{i,beaminds});
            
            Nvals = vertcat(MIZ_DATA.N{i,beaminds});
            
            
            % Usable values to do the correlation
            usable = var1 < 2e5 & Nvals > 100;
            
            if sum(usable) > 1
            
            cmat = corrcoef(var1(usable),var2(usable));
            
            CC(i,j) = cmat(1,2);
            
            end
            
        else
            CC(i,j) = nan;
        end
        
    end
end


%%
Nvals = vertcat(MIZ_DATA.N{:}); 
if IS2_DATA.v6

    SICvals_amsr = vertcat(MIZ_DATA.SIC_amsr{:}); 

end

SICvals = vertcat(MIZ_DATA.SIC{:}); 
LIFvals = vertcat(MIZ_DATA.LIF{:}); 
Hvals = vertcat(MIZ_DATA.H{:}); 
Evals = vertcat(MIZ_DATA.E{:}); 
WAFvals = vertcat(MIZ_DATA.WAF{:}); 

Dvals = (vertcat(MIZ_DATA.D_to_MIZ{:})/1000); 

Dbins = -12.5*1000:25:12.5*1000; 
Bincent = 0.5*(Dbins(1:end-1) + Dbins(2:end));
Bincent(end+1) = Bincent(end) + Bincent(2) - Bincent(1); 

usable = (Nvals > 5) & ~isnan(Dvals); 

Nvals = Nvals(usable); 

if IS2_DATA.v6

SICvals_amsr = SICvals_amsr(usable); 

end

SICvals = SICvals(usable); 
LIFvals = LIFvals(usable); 
Dvals = Dvals(usable); 
Hvals = Hvals(usable); 
Evals = Evals(usable); 
WAFvals = WAFvals(usable);

[binct,~,binval] = histcounts(Dvals,Dbins);

% binval(binval == 0) = length(Dbins); 

upval = @(x) prctile(x,75); 
dnval = @(x) prctile(x,25); 

Nvec = accumarray(binval,1,[length(Dbins) 1],@sum); 


Nsegvec = accumarray(binval,Nvals,[length(Dbins) 1],@nanmean); 
Nsegvar = accumarray(binval,Nvals,[length(Dbins) 1],@std); 

SICvec = accumarray(binval,SICvals,[length(Dbins) 1],@nanmedian); 
SICup = accumarray(binval,SICvals,[length(Dbins) 1],upval); 
SICdn = accumarray(binval,SICvals,[length(Dbins) 1],dnval); 

if IS2_DATA.v6

SICvec_amsr = accumarray(binval,SICvals_amsr,[length(Dbins) 1],@nanmedian); 
SICup_amsr= accumarray(binval,SICvals_amsr,[length(Dbins) 1],upval); 
SICdn_amsr = accumarray(binval,SICvals_amsr,[length(Dbins) 1],dnval); 

end

WAFvec = accumarray(binval,WAFvals,[length(Dbins) 1],@nanmean); 
WAFup = accumarray(binval,WAFvals,[length(Dbins) 1],upval); 
WAFdn = accumarray(binval,WAFvals,[length(Dbins) 1],dnval); 

Hvec = accumarray(binval,Hvals,[length(Dbins) 1],@nanmedian); 
Hup = accumarray(binval,Hvals,[length(Dbins) 1],upval); 
Hdn = accumarray(binval,Hvals,[length(Dbins) 1],dnval); 


LIFvec = accumarray(binval,LIFvals,[length(Dbins) 1],@nanmedian); 
LIFup = accumarray(binval,LIFvals,[length(Dbins) 1],upval); 
LIFdn = accumarray(binval,LIFvals,[length(Dbins) 1],dnval); 


Hvar = accumarray(binval,Hvals,[length(Dbins) 1],@nanstd); 

Evec = accumarray(binval,Evals,[length(Dbins) 1],@nanmean); 
Evar = accumarray(binval,Evals,[length(Dbins) 1],@nanstd); 

%
close

subplot('position',[.1 .4 .8 .5])

jbfill(Bincent,SICup',SICdn',[.4 .4 .4],[1 1 1],1,.3)
hold on

if IS2_DATA.v6

jbfill(Bincent,SICup_amsr',SICdn_amsr',[.8 .2 .2],[1 1 1],1,.3)

end


jbfill(Bincent,LIFup',LIFdn',[.2 .2 .8],[1 1 1],1,.3)
hold on
p1 = plot(Bincent,SICvec,'k','linewidth',2); 

if IS2_DATA.v6

p3 = plot(Bincent,SICvec_amsr,'color',[.8 .2 .2],'linewidth',2); 

end

p2 = plot(Bincent,LIFvec,'color',[.4 .4 .8],'linewidth',2); 


xfirst = find(Nvec > 10,1,'first');

dum = find(Nvec > 10,2,'last');

xlast = dum(1); 

xlimmer = [Bincent(xfirst) Bincent(xlast)]; 
xlimmer = min(abs(xlimmer))*[-1 1];

plot(Bincent,LIFup,'--b'); 
plot(Bincent,LIFdn,'--b'); 
xlim(xlimmer)
grid on; box on; 
ylim([0 1])
ylabel('Ice Fraction','Interpreter','latex')
xline(0,'color',[.2 .2 .2],'linewidth',1)

yyaxis right
p3 = plot(Bincent,Nvec,'color',[.8 .4 .4],'LineWidth',2);
set(gca,'ycolor','k')
ylabel('Number of Samples','interpreter','latex')
xlim(xlimmer)
legend([p1 p2 p3],{'PM-SIC','LIF','Number'},'location','best')

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
subplot('position',[.1 .1 .8 .25])


% jbfill(Bincent,Hup',Hdn',[.4 .4 .4],[1 1 1],1,.3)

jbfill(Bincent,WAFup',WAFdn',[.4 .4 .4],[1 1 1],1,.3)
hold on
r1 = plot(Bincent,WAFvec,'k','linewidth',2); 

plot(Bincent,WAFup,'--k'); 
plot(Bincent,WAFdn,'--k'); 
xlim(xlimmer)
% ylim([0 .5]); 
set(gca,'ylim',[0 max(get(gca,'ylim'))])
ylabel('Wave Affected Fraction','interpreter','latex')
% r2 = yline(.075,'--','color',[.8 .4 .4],'linewidth',1)
xline(0,'color',[.2 .2 .2],'linewidth',1)
grid on; box on;
xlabel('Kilometers Relative to PM-MIZ','interpreter','latex')
% legend([r1 r2],{'Wave Affected Fraction','WMIZ Threshold'})
%% SECTION TITLE
% DESCRIPTIVE TEXT

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

pos = [6.5 5]; 
set(gcf,'windowstyle','normal','position',[0 0 pos],'paperposition',[0 0 pos],'papersize',pos,'units','inches','paperunits','inches');
set(gcf,'windowstyle','normal','position',[0 0 pos],'paperposition',[0 0 pos],'papersize',pos,'units','inches','paperunits','inches');
% print('/Users/chorvat/Library/CloudStorage/Dropbox-Brown/Christopher Horvat/Apps/Overleaf/2024-NASA-ROSES-PM/Proposal/Figures/MIZ-SIC-comp','-dpdf','-r600');
%%