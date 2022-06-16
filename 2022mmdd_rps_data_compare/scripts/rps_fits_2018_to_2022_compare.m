function rps_fits_2018_to_2022_compare()



%%
fd0 = {};
load('z:/dev/rps/rps_beam_fits_cut_2018.mat')
fd.xpol = fd.aparam(:,2);
fd.phi = fd.phi_d;
fd.phi_err = fd.aerr(:,1);
fd0{1} = fd;

load('z:/dev/rps/rps_beam_fits_mirror_persch_cut.mat')
%load('z:/dev/rps/rps_beam_fits_cut.mat')
fd0{2} = fd;
fd = fd0;

load('z:/dev/rps/fpu_data_obs.mat')

figdir = 'C:\Users\James\Documents\GitHub\postings\2022mmdd_rps_data_compare\figs';


xpols = NaN(2,2640);
for yearind = 1:2
    fd0 = fd{yearind};
    for chind = 1:length(fd0.ch)
        ci = find(fd0.ch==chind);
        
        if ~isempty(ci)
            val = reshape(fd0.xpol(ci),[],1);
            err = reshape(fd0.phi_err(ci),[],1);
            err(isnan(val))=1e10;
            
            xpols(yearind,chind) = wmean(val,1./err,1);
        end
    end
end
% Use only common data:
indxpol = isnan(diff(xpols,1));




phis = NaN(2,2640);
for yearind = 1:2
    fd0 = fd{yearind};
    for chind = 1:length(fd0.ch)
        ci = find(fd0.ch==chind);
        
        if ~isempty(ci)
            ch_chi = atand(tand(p.theta(chind)+p.chi(chind)));
            
            if abs(ch_chi) < 5
                val = reshape(atand(tand(fd0.phi(ci))),[],1);
            else
                val = reshape(atand(tand(fd0.phi(ci)-90))+90,[],1);
            end
            err = reshape(fd0.phi_err(ci),[],1);
            err(isnan(val))=1e10;
            
            phis(yearind,chind) = wmean(val,1./err,1);
        end
    end
end
% Use only common data:
ind = isnan(diff(phis,1)) | abs(diff(phis,1))>1 | indxpol;
xpols(:,ind) = NaN;
phis(:,ind) = NaN;


% Calc the pair-diff pols
[phi_diffs, xpol_diffs] = deal(NaN(2,length(p_ind.a)));
for yearind = 1:2
    for chind = 1:length(p_ind.a)
        ia = p_ind.a(chind);
        ib = p_ind.b(chind);
        if ~isnan(phis(yearind,ia))
            
            Q = ([cosd(2*phis(yearind,ia))-xpols(yearind,ia)*cosd(2*phis(yearind,ia))]-...
                [cosd(2*phis(yearind,ib))-xpols(yearind,ib)*cosd(2*phis(yearind,ib))])./...
                (2+xpols(yearind,ia)+xpols(yearind,ib));
            U = ([sind(2*phis(yearind,ia))-xpols(yearind,ia)*sind(2*phis(yearind,ia))]-...
                [sind(2*phis(yearind,ib))-xpols(yearind,ib)*sind(2*phis(yearind,ib))])./...
            (2+xpols(yearind,ia)+xpols(yearind,ib));
        
            phi_diffs(yearind,chind) = atan2(U,Q)./2*180/pi;
            xpol_diffs(yearind,chind) = 1-sqrt(Q.^2+U.^2);
        
        end
    end
end
ind = p.mce(p_ind.a)==0;
phi_diffs(:,ind) = atand(tand(phi_diffs(:,ind)-90));

%%
fig = figure(1);
fig.Position(3:4) = [1500 500];


lims = 0.03;
res = 50;
pols = {'a','b'};
Nlims = 250;

for polind = 1:2
    clf;
    subplot(1,3,1)
    edges = (-1:2/res:1)*lims;
    N = histc(xpols(1,p_ind.(pols{polind})),edges);
    bar(edges,N,'histc')
    grid on
    xlabel('Xpol Efficiency')
    title('All RPS2018 Xpol Efficiency')
    xlim([-1 1]*lims)
    ylim([0 Nlims])
    
    subplot(1,3,2)
    edges = (-1:2/res:1)*lims;
    N = histc(xpols(2,p_ind.(pols{polind})),edges);
    bar(edges,N,'histc')
    grid on
    xlabel('Xpol Efficiency')
    title('All RPS2022 Xpol Efficiency')
    xlim([-1 1]*lims)
    ylim([0 Nlims])
    
    subplot(1,3,3)
    edges = (-1:2/res:1)*lims;
    N = histc(diff(xpols(:,p_ind.(pols{polind})),1),edges);
    bar(edges,N,'histc')
    grid on
    xlabel('Xpol Efficiency')
    title({'All RPS2018 Minus RPS2022',...
        sprintf('Mean: %1.3f STD: %1.3f',nanmean(diff(xpols(:,p_ind.(pols{polind})),1)),...
        nanstd(diff(xpols(:,p_ind.(pols{polind})),1)))})
    xlim([-1 1]*lims)
    ylim([0 Nlims])
    
    fname = sprintf('xpol_hist_%s.png',pols{polind});
    saveas(fig,fullfile(figdir,fname))
end
%%

if 0
    fig = figure(2);
    fig.Position(3:4) = [900 700];
    clf;
    
    plot_tiles(xpols(1,:),p,'fig',fig,'clim',[-1 1]*0.02,'title','2018 Xpol Effs');
    
    fig = figure(3);
    fig.Position(3:4) = [900 700];
    clf;
    plot_tiles(xpols(2,:),p,'fig',fig,'clim',[-1 1]*0.02,'title','2022 Xpol Effs');
    
    %
    fig = figure(4);
    fig.Position(3:4) = [900 700];
    clf;
    plot_tiles(diff(xpols,1),p,'fig',fig,'clim',[-1 1]*0.005,'title','Xpol Effs 2018-2022');
end

%%
fig = figure(5);
fig.Position(3:4) = [900,600];
clf; hold on;

lims = 0.03;
offs = 0;

h = scatterhist(xpols(1,p_ind.rgl100a),xpols(1,p_ind.rgl100b),'kernel','on');
xlim([-1 1]*lims)
ylim([-1 1]*lims)
grid on
title(' RPS2018 Pol A vs Pol B')
xlabel({'Xpol_{2018} [Degrees]','Pol A detectors'})
ylabel({'Xpol_{2018} [Degrees]','Pol B detectors'})

fname = 'xpol_scatterhist_2018.png';
saveas(fig,fullfile(figdir,fname))

fig = figure(6);
fig.Position(3:4) = [900,600];
clf; hold on;

lims = 0.03;
offs = -0;

h = scatterhist(xpols(2,p_ind.rgl100a),xpols(2,p_ind.rgl100b),'kernel','on');
xlim([-1 1]*lims)
ylim([-1 1]*lims)
grid on
title(' RPS2022 Pol A vs Pol B')
xlabel({'Xpol_{2022} [Degrees]','Pol A detectors'})
ylabel({'Xpol_{2022} [Degrees]','Pol B detectors'})

fname = 'xpol_scatterhist_2022.png';
saveas(fig,fullfile(figdir,fname))


fig = figure(7);
fig.Position(3:4) = [900,600];
clf; hold on;

h = scatterhist(xpols(1,p_ind.rgl100a)-xpols(2,p_ind.rgl100a),xpols(1,p_ind.rgl100b)-xpols(2,p_ind.rgl100b),'kernel','on');
ylim([-1 1]*0.03)
xlim([-1 1]*0.03)
grid on
title(' RPS2018 minus RPS2022 Pol A vs Pol B')
xlabel({'Xpol_{2018} - Xpol_{2022} [Degrees]','Pol A detectors'})
ylabel({'Xpol_{2018} - Xpol_{2022} [Degrees]','Pol B detectors'})

fname = 'xpol_scatterhist_diff.png';
saveas(fig,fullfile(figdir,fname))

%%

fig = figure(8);
fig.Position(3:4) = [1500 500];
clf;

for polind = 1:2
    clf;
    lims = 3;
    res = 50;
    Nlim = 130;
    offs = {-2.3, 87};
    subplot(1,3,1)
    edges = (-1:2/res:1)*lims+offs{polind};
    N = histc(phis(1,p_ind.(pols{polind})),edges);
    bar(edges,N,'histc')
    grid on
    xlabel('\phi [Degrees]')
    title('All RPS2018 Pol Angles')
    xlim([-1 1]*lims+offs{polind})
    ylim([0 Nlim])
    
    subplot(1,3,2)
    edges = (-1:2/res:1)*lims+offs{polind};
    N = histc(phis(2,p_ind.(pols{polind})),edges);
    bar(edges,N,'histc')
    grid on
    xlabel('\phi [Degrees]')
    title('All RPS2022 Pol Angles')
    xlim([-1 1]*lims+offs{polind})
    ylim([0 Nlim])
    
    
    lims = 1;
    res = 50;
    offs = -0;
    
    
    subplot(1,3,3)
    edges = (-1:2/res:1)*lims+offs;
    N = histc(diff(phis(:,p_ind.(pols{polind})),1),edges);
    bar(edges,N,'histc')
    grid on
    xlabel('\phi [Degrees]')
    title({'All RPS2018 Minus RPS2022',...
        sprintf('Mean: %1.3f STD: %1.3f',nanmean(diff(phis(:,p_ind.(pols{polind})),1)),...
        nanstd(diff(phis(:,p_ind.(pols{polind})),1)))})
    xlim([-1 1]*lims+offs)
    ylim([0 Nlim])
    
    fname = sprintf('phi_hist_%s.png',pols{polind});
    saveas(fig,fullfile(figdir,fname))
end

%%
fig = figure(9);
fig.Position(3:4) = [900,600];
clf; hold on;

lims = 1.5;
offs = -3;

h = scatterhist(phis(1,p_ind.rgl100a),phis(1,p_ind.rgl100b),'kernel','on');
xlim([-1 1]*lims-2.3)
ylim([-1 1]*lims+87)
grid on
title(' RPS2018 Pol A vs Pol B')
xlabel({'\phi_{2018} [Degrees]','Pol A detectors'})
ylabel({'\phi_{2018} [Degrees]','Pol B detectors'})

fname = 'phi_scatterhist_2018.png';
saveas(fig,fullfile(figdir,fname))

fig = figure(10);
fig.Position(3:4) = [900,600];
clf; hold on;

lims = 1.5;
offs = -3;

h = scatterhist(phis(2,p_ind.rgl100a),phis(2,p_ind.rgl100b),'kernel','on');
xlim([-1 1]*lims-2.3)
ylim([-1 1]*lims+87)
grid on
title(' RPS2022 Pol A vs Pol B')
xlabel({'\phi_{2022} [Degrees]','Pol A detectors'})
ylabel({'\phi_{2022} [Degrees]','Pol B detectors'})

fname = 'phi_scatterhist_2022.png';
saveas(fig,fullfile(figdir,fname))


fig = figure(11);
fig.Position(3:4) = [900,600];
clf; hold on;

h = scatterhist(phis(1,p_ind.rgl100a)-phis(2,p_ind.rgl100a),phis(1,p_ind.rgl100b)-phis(2,p_ind.rgl100b),'kernel','on');
ylim([-1 1]*1)
xlim([-1 1]*1)
grid on
title(' RPS2018 minus RPS2022 Pol A vs Pol B')
xlabel({'\phi_{2018} - \phi_{2022} [Degrees]','Pol A detectors'})
ylabel({'\phi_{2018} - \phi_{2022} [Degrees]','Pol B detectors'})

fname = 'phi_scatterhist_diff.png';
saveas(fig,fullfile(figdir,fname))

%% scatterhist of just 2018 vs 2022 phis

phis_mce = phis;
ind = p.mce==0 & ismember((1:2640)',p_ind.a);
phis_mce(:,ind) = phis_mce(:,ind)-90;
ind = p.mce==0 & ismember((1:2640)',p_ind.b);
phis_mce(:,ind) = atand(tand(phis_mce(:,ind)-90));
fig = figure(12);
fig.Position(3:4) = [900,600];
clf; hold on;

offs = {-2.3 87.5};
scaling = {2, 2};
polind = {p_ind.a, p_ind.b,...
    };
pols = {'A','B'};
for pltind = 1:2
    h = scatterhist(phis_mce(1,polind{pltind}),phis_mce(2,polind{pltind}),'kernel','on');
    %scatter(phis_mce(1,polind{pltind}),phis_mce(2,polind{pltind}),14,p.mce(polind{pltind}),'filled')
    hold(h(1),'on')
    plot(h(1),[-1 1]*scaling{pltind}+offs{pltind},[-1 1]*scaling{pltind}+offs{pltind},'k--')
    xlim([-1 1]*scaling{pltind}+offs{pltind})
    ylim([-1 1]*scaling{pltind}+offs{pltind})
    grid on
    title(sprintf('RPS2018 Vs RPS2022 Pol %s',pols{pltind}))
    xlabel({'\phi_{2018} [Degrees]',sprintf('Pol %s detectors',pols{pltind})})
    ylabel({'\phi_{2022} [Degrees]',sprintf('Pol %s detectors',pols{pltind})})
    
    fname = sprintf('phi_scatterhist_pol_%s.png',pols{pltind});
    saveas(fig,fullfile(figdir,fname))
end

%% scatterhist of just 2018 vs 2022 phis
% separated by 0/90 instead of designated pol A/B

fig = figure(14);
fig.Position(3:4) = [900,600];
clf; hold on;

offs = {-2.3 87.5};
scaling = {2, 2};
polind = {[intersect(find(p.mce~=0),p_ind.a); intersect(find(p.mce==0),p_ind.b)],...
    [intersect(find(p.mce~=0),p_ind.b); intersect(find(p.mce==0),p_ind.a)]...
    };
pols = {'0','90'};
for pltind = 1:2
    h = scatterhist(phis(1,polind{pltind}),phis(2,polind{pltind}),'kernel','on');
    %scatter(phis_mce(1,polind{pltind}),phis_mce(2,polind{pltind}),14,p.mce(polind{pltind}),'filled')
    hold(h(1),'on')
    plot(h(1),[-1 1]*scaling{pltind}+offs{pltind},[-1 1]*scaling{pltind}+offs{pltind},'k--')
    xlim([-1 1]*scaling{pltind}+offs{pltind})
    ylim([-1 1]*scaling{pltind}+offs{pltind})
    grid on
    title(sprintf('RPS2018 Vs RPS2022 Pol %s',pols{pltind}))
    xlabel({'\phi_{2018} [Degrees]',sprintf('Pol %s detectors',pols{pltind})})
    ylabel({'\phi_{2022} [Degrees]',sprintf('Pol %s detectors',pols{pltind})})
    
    fname = sprintf('phi_scatterhist_pol_%s.png',pols{pltind});
    saveas(fig,fullfile(figdir,fname))
end



%% scatterhist of just 2018 vs 2022 xpol

fig = figure(13);
fig.Position(3:4) = [900,600];
clf; hold on;

offs = {0 0 0 0};
scaling = {0.04, 0.04, 0.04, 0.04};
polind = {p_ind.a p_ind.b,...
    [intersect(find(p.mce~=0),p_ind.a); intersect(find(p.mce==0),p_ind.b)],...
    [intersect(find(p.mce~=0),p_ind.b); intersect(find(p.mce==0),p_ind.a)]...
    };
pols = {'A','B','0','90'};
for pltind = 1:4
    h = scatterhist(xpols(1,polind{pltind}),xpols(2,polind{pltind}),'kernel','on');
    %scatter(phis_mce(1,polind{pltind}),phis_mce(2,polind{pltind}),14,p.mce(polind{pltind}),'filled')
    hold(h(1),'on')
    plot(h(1),[-1 1]*scaling{pltind}+offs{pltind},[-1 1]*scaling{pltind}+offs{pltind},'k--')
    xlim([-1 1]*scaling{pltind}+offs{pltind})
    ylim([-1 1]*scaling{pltind}+offs{pltind})
    grid on
    title(sprintf('RPS2018 Vs RPS2022 Pol %s',pols{pltind}))
    xlabel({'\epsilon_{2018} [Degrees]',sprintf('Pol %s detectors',pols{pltind})})
    ylabel({'\epsilon_{2022} [Degrees]',sprintf('Pol %s detectors',pols{pltind})})
    
    fname = sprintf('xpol_scatterhist_pol_%s.png',pols{pltind});
    saveas(fig,fullfile(figdir,fname))
end

%% scatterhist of phi-diffs



offs = {-2.3, 0};
scaling = {2, 0.04};
vals = {phi_diffs,xpol_diffs};
valnames = {'phi','xpol'};
valtitles = {'\phi','xpol'};
pols = {'diff','diff'};
for pltind = 1:2
    fig = figure(16+pltind);
fig.Position(3:4) = [900,600];
clf; hold on;
    h = scatterhist(vals{pltind}(1,:),vals{pltind}(2,:),'kernel','on');
    %scatter(phis_mce(1,polind{pltind}),phis_mce(2,polind{pltind}),14,p.mce(polind{pltind}),'filled')
    hold(h(1),'on')
    plot(h(1),[-1 1]*scaling{pltind}+offs{pltind},[-1 1]*scaling{pltind}+offs{pltind},'k--')
    xlim([-1 1]*scaling{pltind}+offs{pltind})
    ylim([-1 1]*scaling{pltind}+offs{pltind})
    grid on
    title(sprintf('RPS2018 Vs RPS2022 %s-diff',valtitles{pltind}))
    xlabel({sprintf('%s_{2018} [Degrees]',valtitles{pltind}),sprintf('Pol %s detectors',pols{pltind})})
    ylabel({sprintf('%s_{2022} [Degrees]',valtitles{pltind}),sprintf('Pol %s detectors',pols{pltind})})
    
    fname = sprintf('%s_scatterhist_pol_%s.png',valnames{pltind},pols{pltind});
    %saveas(fig,fullfile(figdir,fname))
end




