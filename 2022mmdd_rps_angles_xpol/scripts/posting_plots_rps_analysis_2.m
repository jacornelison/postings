function posting_plots_rps_analysis_2()

%% Run this first
addpath('z:/dev/rps/')
addpath('z:/pipeline/beammap')
addpath('z:/pipeline/util')
addpath('z:/dev/diff_polarization/')
addpath('z:/dev')
figdir = fullfile('C:','Users','James','Documents','GitHub','postings','2022mmdd_rps_angles_xpol','figs','');

%% Then this
%clear all
close all
clc

p22 = load('z:/dev/rps/fpu_data_obs.mat'); % pointing from 2022
p18 = load('z:/dev/rps/fpu_data_B18.mat'); % RGL's from B18
p = p22.p;
p_ind = p18.p_ind;

% shortcuts for 0/90 orientiations
inda = p_ind.a;
indb = p_ind.b;

ind0 = [inda(ismember(inda,find(p.mce~=0))) ...
    indb(ismember(indb,find(p.mce==0)))];
ind90 = [indb(ismember(indb,find(p.mce~=0))) ...
    inda(ismember(inda,find(p.mce==0)))];

% Load type 9 data
%load('z:/dev/rps/type9_fit_dat_mirrorfitted_cut.mat')
fd_type9 = load('z:/dev/rps/rps_beam_fits_type9_21feb_rerun.mat');
fd_type9 = fd_type9.fd;
fd_type9 = rps_cut_fitdata(fd_type9,p,p_ind,false);
load('z:/dev/rps/sch_type9.mat')
fd_type9 = get_pair_params(fd_type9,ind0,ind90);
[fd_type9, ~,~,~,~] = get_pol_params_per_obs(fd_type9,p);


% Load 2018 data
load('z:/dev/rps/rps_beam_fits_cut_2018.mat')
fd.xpol = fd.aparam(:,2);
fd.phi = fd.phi_d;
fd.phi_err = fd.aerr(:,1);
fd.schnum = fd.sch;
fd.rowind = fd.row;
fd_2018 = fd;
fd_2018 = get_pair_params(fd_2018,ind0,ind90);
[fd_2018, phis_2018, phi_pair_2018, xpols_2018, poleff_pair_2018,~,~,~] = get_pol_params_per_obs(fd_2018,p);

%
clc
load('z:/dev/rps/rps_obs_info.mat')
%load('z:/dev/rps/rps_beam_fits_type5_withbparam.mat')
load('z:/dev/rps/rps_beam_fits_type5_21feb_rerun.mat');
fd = rps_cut_fitdata(fd,p,p_ind);%,1,figdir);
%fd = rps_cut_fitdata(fd,p,[]);%,1,figdir);
%
fd = get_pair_params(fd,ind0,ind90);
[fd, phis, phi_pair, xpols, poleff_pair,n1s,n2s,amps] = get_pol_params_per_obs(fd,p,scheds);

load('z:/pipeline/beammap/viridis_cm.mat')
load('z:/dev/rps/sch_type5.mat')
load('z:/dev/rps/pm.mat')


% Load no tilt data
FD = load('z:/dev/rps/rps_beam_fits_type5_notilt_cut.mat');
fd_notilt = FD.fd; clear FD;
fd_notilt = get_pair_params(fd_notilt,ind0,ind90);
[fd_notilt, phis_nt, phi_pair_nt, xpols_nt, poleff_pair_nt,~,~,~] = get_pol_params_per_obs(fd_notilt,p,scheds);


cmlines = colormap('lines');

[fd.theta, fd.r] = pol2cart(fd.x,fd.y);
fd.theta = fd.theta*180/pi;


% Expected pol params from FPU data
phis_exp = atand(tand(p.chi+p.chi_thetaref))';
eps_exp = p.epsilon';

[phi_pair_exp, poleff_pair_exp] = deal(NaN(1,2640));
[phi_pair_exp(ind0), poleff_pair_exp(ind0)] = calc_pair_diff_pol(phis_exp(ind0),phis_exp(ind90),eps_exp(ind0),eps_exp(ind90));


%% Fig 2.1 Angle vs. channel

set(groot,'defaultLineMarkerSize',12)
yrnames = {'2022','2018','exp'};
V = {phis, phis_2018, phis_exp};
V2 = {phi_pair, phi_pair_2018, phi_pair_exp};
for yearind = 1:length(V)


    fig = figure(21);
    fig.Position(3:4) = [900 500];
    clf;

    clc
    cmlines = colormap('lines');
    for dkind = 1:size(V{yearind},1)
        X = {ind0, ind90, ind0, ind0};
        Y = {V{yearind}(dkind,ind0), V{yearind}(dkind,ind90), ...
            V{yearind}(dkind,ind0)-V{yearind}(dkind,ind90)+90 ,V2{yearind}(dkind,ind0)};
        ylab = {'\phi_0','\phi_{90}','\phi_0-\phi_{90}+90','\phi_{pair}'};
        ylims = {[-5 1], [85 91], [-3 3],[-5 1]};

        for pltind = 1:length(X)


            subplot(length(X),1,pltind)
            hold on



            plot(X{pltind},Y{pltind},'.','Color',cmlines(dkind,:))
            grid on
            ylabel([ylab{pltind} ' [Degrees]'])

            for mceind = 0:3
                chind = max(find(p.mce==mceind));
                plot([1 1]*chind,ylims{pltind},'k--')
            end
            ylim(ylims{pltind})
            xlim([-50 2700])
        end

    end
    xlabel('Channel Number')
    sgtitle('Angles Vs. Channel')
    fname = sprintf('phi_vs_chan_090_%s.png',yrnames{yearind});
    saveas(fig,fullfile(figdir,fname))
end


%% Fig 2.1.2 Xpol vs. channel

yrnames = {'2022','2018','exp'};
V = {xpols, xpols_2018,eps_exp};
V2 = {poleff_pair,poleff_pair_2018,poleff_pair_exp};
for yearind = 1:length(V)


    fig = figure(22);
    fig.Position(3:4) = [900 500];
    clf;


    cmlines = colormap('lines');
    for dkind = 1:size(V{yearind},1)
        X = {ind0, ind90, ind0};
        Y = {V{yearind}(dkind,ind0), V{yearind}(dkind,ind90), V2{yearind}(dkind,ind0)};
        ylab = {'\epsilon_0','\epsilon_{90}','1-\epsilon_{pair}'};
        ylims = {[-1 1]*0.05, [-1 1]*0.05, [-0.1 2]*0.001};

        for pltind = 1:length(X)


            subplot(3,1,pltind)
            hold on

            plot(X{pltind},Y{pltind},'.','Color',cmlines(dkind,:))
            grid on
            ylabel([ylab{pltind} ''])

            for mceind = 0:3
                chind = max(find(p.mce==mceind));
                plot([1 1]*chind,ylims{pltind},'k--')
            end

            ylim(ylims{pltind})
            xlim([-50 2700])
        end

    end
    xlabel('Channel Number')
    sgtitle('Xpol Vs. Channel')
    fname = sprintf('xpol_vs_chan_090_%s.png',yrnames{yearind});
    saveas(fig,fullfile(figdir,fname))
end

%% Fig 2.1 N1 vs. channel

set(groot,'defaultLineMarkerSize',12)
yrnames = {'2022'};
V = {n1s};

for yearind = 1:length(V)
    fig = figure(21);
    fig.Position(3:4) = [900 500];
    clf;

    clc
    cmlines = colormap('lines');
    for dkind = 1:size(V{yearind},1)
        X = {ind0, ind90};
        Y = {V{yearind}(dkind,ind0), V{yearind}(dkind,ind90)};
        ylab = {'N1_0','N1_{90}'};
        ylims = {[-1 1]*0.05, [-1 1]*0.05};

        for pltind = 1:length(X)


            subplot(length(X),1,pltind)
            hold on



            plot(X{pltind},Y{pltind},'.','Color',cmlines(dkind,:))
            grid on
            ylabel([ylab{pltind} ''])

            for mceind = 0:3
                chind = max(find(p.mce==mceind));
                plot([1 1]*chind,ylims{pltind},'k--')
            end
            ylim(ylims{pltind})
            xlim([-50 2700])
        end

    end
    xlabel('Channel Number')
    sgtitle('N1 Vs. Channel')
    fname = sprintf('N1_vs_chan_090_%s.png',yrnames{yearind});
    saveas(fig,fullfile(figdir,fname))
end

%% Fig 2.1 N1 vs. channel

set(groot,'defaultLineMarkerSize',12)
yrnames = {'2022'};
V = {n2s};

for yearind = 1:length(V)
    fig = figure(21);
    fig.Position(3:4) = [900 500];
    clf;

    clc
    cmlines = colormap('lines');
    for dkind = 1:size(V{yearind},1)
        X = {ind0, ind90};
        Y = {V{yearind}(dkind,ind0), V{yearind}(dkind,ind90)};
        ylab = {'N2_0','N2_{90}'};
        ylims = {[-1 1]*0.05, [-1 1]*0.05};

        for pltind = 1:length(X)


            subplot(length(X),1,pltind)
            hold on



            plot(X{pltind},Y{pltind},'.','Color',cmlines(dkind,:))
            grid on
            ylabel([ylab{pltind} ''])

            for mceind = 0:3
                chind = max(find(p.mce==mceind));
                plot([1 1]*chind,ylims{pltind},'k--')
            end
            %ylim(ylims{pltind})
            xlim([-50 2700])
        end

    end
    xlabel('Channel Number')
    sgtitle('N2 Vs. Channel')
    fname = sprintf('N2_vs_chan_090_%s.png',yrnames{yearind});
    saveas(fig,fullfile(figdir,fname))
end

%% Fig 2.1 Amp vs. channel

set(groot,'defaultLineMarkerSize',12)
yrnames = {'2022'};
V = {amps};

for yearind = 1:length(V)
    fig = figure(21);
    fig.Position(3:4) = [900 500];
    clf;

    clc
    cmlines = colormap('lines');
    for dkind = 1:size(V{yearind},1)
        X = {ind0, ind90};
        Y = {V{yearind}(dkind,ind0), V{yearind}(dkind,ind90)};
        ylab = {'A_0 [ADU]','A_{90} [ADU]'};
        ylims = {[0 1]*250, [0 1]*250};

        for pltind = 1:length(X)


            subplot(length(X),1,pltind)
            hold on



            plot(X{pltind},Y{pltind},'.','Color',cmlines(dkind,:))
            grid on
            ylabel([ylab{pltind} ''])

            for mceind = 0:3
                chind = max(find(p.mce==mceind));
                plot([1 1]*chind,ylims{pltind},'k--')
            end
            ylim(ylims{pltind})
            xlim([-50 2700])
        end

    end
    xlabel('Channel Number')
    sgtitle('Amp Vs. Channel')
    fname = sprintf('amp_vs_chan_090_%s.png',yrnames{yearind});
    saveas(fig,fullfile(figdir,fname))
end

%% Angle VS. DK per tile


trow = [2 3 4 1 2 3 4 5 1 2 3 4 5 1 2 3 4 5 2 3 4];
tcol = [5 5 5 4 4 4 4 4 3 3 3 3 3 2 2 2 2 2 1 1 1];
rc2idx = reshape(1:25,5,5)';
fig = figure(1);
fig.Position(3:4) = [750 650];
clf;
t = tiledlayout(5,5);
t.Padding = 'tight';
t.TileSpacing = 'tight';
P = phi_pair-repmat(nanmedian(phi_pair,1),10,1);
[b bi] = sort(dks);
for tind = 1:20
    idx = p.tile==tind;
    nexttile(rc2idx(trow(tind),tcol(tind)))
    hold on
    plot(dks(bi),nanmean(P(bi,idx),2),'.')
    errorbar(dks(bi),nanmean(P(bi,idx),2),nanstd(P(bi,idx),[],2),'CapSize',0)
    grid on
    lims = [-1 1]*0.11;
    ylim(lims)
    xlim([-180 180])
    pbaspect([1 1 1])
    ttl = title(sprintf('Tile %i',tind));
    ttl.Position(2) = 0.075;
    if tind == 18
        xb = xlabel('DK [Deg]');
        xb.Position(2) = -0.075;
        ylabel({'\phi_{pair}-md(\phi) [Deg]'})
    end
end

%% Angle VS. obs
%load('z:/dev/rps/rps_phi_final_2022.mat')
P = phi_pair-repmat(nanmedian(phi_pair,1),10,1);
fig = figure(2);
fig.Position(3:4) = [450 420];
clf; hold on
errorbar(1:10,nanmean(P,2),nanstd(P,[],2),'.','CapSize',0,'LineWidth',1.25)
plot(1:10,nanmean(P,2),'.','MarkerSize',14)
%plot([-180 180],polyval(C,[-180 180]),'k--')

grid on
lims = [-1 1]*0.085;
ylim(lims)
xlim([0 10])
xlabel('Obs')
ylabel({'\phi_{pair}-md(\phi) [Deg]'})
title({'\phi_{pair} Vs. Obs','Median-subtracted per pair, avg''d per obs'})
pbaspect([1 1 1])
fname = 'phi_pair_vs_obs_medsub';
saveas(fig,fullfile(figdir,fname),'png')

%% Angle VS. DK overall
%load('z:/dev/rps/rps_phi_final_2022.mat')
P = phi_pair-repmat(nanmedian(phi_pair,1),10,1);
fig = figure(2);
fig.Position(3:4) = [450 420];
clf; hold on
C = polyfitweighted(dks,nanmean(P,2)',1,nanstd(P,[],2)');
errorbar(dks(bi),nanmean(P(bi,:),2),nanstd(P(bi,:),[],2),'CapSize',0)
plot(dks(bi),nanmean(P(bi,:),2),'.','MarkerSize',14)
%plot([-180 180],polyval(C,[-180 180]),'k--')

grid on
lims = [-1 1]*0.085;
ylim(lims)
xlim([-180 180])
xlabel('Dk [Deg]')
ylabel({'\phi_{pair}-md(\phi) [Deg]'})
title({'\phi_{pair} Vs. Dk','Median-subtracted per pair, avg''d per obs'})
pbaspect([1 1 1])
fname = 'phi_pair_vs_dk_medsub';
saveas(fig,fullfile(figdir,fname),'png')


%%

fig = figure(2);
fig.Position(3:4) = [450 420];
clf; hold on
edges = (-1:1/8:1)*0.06;
N = histc(nanmean(P,2),edges);
bar(edges,N,'histc')
grid on
title('Mean Med-Sub \phi_{pair} per-obs')
xlabel({'\phi_{pair}-md(\phi) [Deg]'})
ylabel('N')
pbaspect([1 1 1])
fname = 'phi_pair_avg_perobs_medsub';
saveas(fig,fullfile(figdir,fname),'png')

%% Fig. 2.1.3 Averaged tile plots

vals = {phi_pair poleff_pair, n1s, n2s, amps};
valname = {'phi','xpol','n1','n2','amp'};

labs = {'\phi_{pair} [Deg]','1-\epsilon_{pair}','N1','N2','A [ADU]'};
ttls = {'\phi_{pair} Tile Plot','\phi_{pair} Tile Plot - Median Subtracted';...
    '1-\epsilon_{pair} Tile Plot','1-\epsilon_{pair} Tile Plot - Median Subtracted';...
    'N1 Tile Plot','N1 Tile Plot - Median Subtracted';...
    'N2 Tile Plot','N2 Tile Plot - Median Subtracted';...
    'Amp Tile Plot','Amp Tile Plot - Median Subtraced';...
    };

pltname = {'','_medsub'};
lims = {[-3.25 -1], [-1 1]*0.25;...
    [0 1]*2e-4, [-1 1]*0.0001;...
    [-1 1]*8e-3, [-1 1]*8e-3;...
    [-1 1]*8e-3, [-1 1]*8e-3;...
    [50 200], [-1 1]*60;...
    };

clc
for valind = 1:length(vals)
    for pltind = 1:length(pltname)
        V = nanmean(vals{valind},1);


        fig = figure(22+pltind);
        fig.Position(3:4) = [900 800];
        if pltind == 2
            for tileind = 1:20
                ind = p.tile==tileind;
                V(ind) = V(ind)-nanmedian(V(ind));
            end
        end
        if ismember(valind, 1:2)
            V(ind90) = V(ind0);
        end

        plot_tiles(V,p,'fig',fig,'pair','mean','clim',lims{valind,pltind},'clab',labs{valind},'title',ttls{valind,pltind});
        colormap(cm)
        fname = sprintf('%s_tile_plot%s.png',valname{valind},pltname{pltind});
        saveas(fig,fullfile(figdir,fname))

    end
end


%% Compare Tilt Correction

if ~exist('tiltcorrfig','var') | ~ishandle(tiltcorrfig)
    count = 1;
    while ishandle(count)
        count = count+1;
    end
end
tiltcorrfig = figure(count);
tiltcorrfig.Position(3:4) = [750 400];
clf; hold on;
t = tiledlayout(1,2);
t.TileSpacing = 'none';
t.Padding = 'loose';

clc
%edges = (-1:0.075:1)*0.6;
edges = (-1:0.075:1)*0.2;

V1 = reshape(phi_pair_nt-nanmedian(phi_pair_nt,1),[],1);
V1 = V1;
V2 = reshape(phi_pair-nanmedian(phi_pair,1),[],1);
V2 = V2;
idx = ~isnan(V1.*V2) & V2~=0 & V1~=0;

nexttile(2)
V = V2(idx);
N = histc(V,edges);
Nmax = max(N);
b = bar(edges,N,'histc');
b.FaceColor = cmlines(1,:);
M = nanmean(V);
S = nanstd(V);
L = length(find(~isnan(V)));
title({'With Tilt Meter Correction',...
    sprintf('M: %0.3f, STD: %0.3f, N: %0i',M,S,L)})
pbaspect([1 1 1])
xlim([edges(1) edges(end)])
grid on
ylim([-0.1 Nmax*1.1])
ax = gca;
ax.YTickLabel = {};
xlabel('\phi_{pair} - median(\phi_{pair}) [Deg]')

nexttile(1)
V = V1(idx);
N = histc(V,edges);
b = bar(edges,N,'histc');
b.FaceColor = cmlines(1,:);
M = nanmean(V);
S = nanstd(V);
L = length(find(~isnan(V)));
title({'No Tilt Meter Correction',...
    sprintf('M: %0.3f, STD: %0.3f, N: %0i',M,S,L)})
pbaspect([1 1 1])
xlim([edges(1) edges(end)])
grid on
ylim([-0.1 Nmax*1.1])
xlabel('\phi_{pair} - median(\phi_{pair}) [Deg]')
ylabel('N')

fname = 'tilt_corr_hist';
saveas(tiltcorrfig,fullfile(figdir,fname),'png')

%% Find the maximum STD for 2022 data subsets

fd0 = fd;
%idx = fd.obsnum~=1;
%fd = structcut(fd,idx);
Niter = 1000;
nobs = length(find(~isnan(unique(fd.obsnum))));
[M, S, N] = deal(NaN(Niter,1));
idx = NaN(Niter,ceil(nobs/2));
for iterind = 1:Niter
    ind1 = randperm(nobs,ceil(nobs/2));
    ind2 = setdiff(1:nobs,ind1);

    Vdiff = nanmean(phi_pair(ind1,:),1)-nanmean(phi_pair(ind2,:),1);
    M(iterind) = nanmean(Vdiff);
    S(iterind) = nanstd(Vdiff);
    N(iterind) = length(find(~isnan(Vdiff)));
    idx(iterind,:) = ind1;
end
[~, mi] = max(abs(M.*sqrt(N)./S));
idx_max = idx(mi,:);
[~, mi] = min(abs(M.*sqrt(N)./S));
idx_min = idx(mi,:);

% Figure 3.2.1 Consistency Checks Part 2

clc
V0 = {phi_pair; poleff_pair};
V = {phi_pair,phi_pair,phi_pair_2018,phi_pair_exp;...
    poleff_pair,poleff_pair,poleff_pair_2018,poleff_pair_exp};
lims1 = {[-3.5 -1], [-3.5 -1], [-3.5 -1], [-3.5 -1];...
    [0 6e-4], [0 6e-4], [0 6e-4], [0 6e-4]};
lims2 = {[-1 1]*0.4,[-1 1]*0.4, [-1 1]*0.4, [-1 1]*2;
    [-1 1]*0.6e-3,[-1 1]*0.6e-3,[-1 1]*0.6e-3,[-1 1]*0.6e-3};
ttls = {'2022','2022','2018','B18 FPU Data'};
pltnames = {'2022_min','2022_max','2018','fpu'};
valnames = {'phi','xpol'};
valunits = {'\phi_{pair} [Degrees]','1-Poleff'};

for valind = 1%1:size(V,1)

    for pltind = 4%1:size(V,2)
        V1ttl = ttls{1};
        V2ttl = ttls{pltind};
        if pltind == 1
            ind1 = idx_min;
            ind2 = find(~ismember(1:10,ind1));
            V1 = V0{valind}(ind1,:);
            V2 = V{valind,pltind}(ind2,:);
            V1ttl = [V1ttl '\_SUB1'];
            V2ttl = [V2ttl '\_SUB2'];
        elseif pltind == 2
            ind1 = idx_max;
            ind2 = find(~ismember(1:10,ind1));
            V1 = V0{valind}(ind1,:);
            V2 = V{valind,pltind}(ind2,:);
            V1ttl = [V1ttl '\_SUB1'];
            V2ttl = [V2ttl '\_SUB2'];
        else
            V1 = V0{valind};
            V2 = V{valind,pltind};
        end
        V1 = nanmean(V1,1);
        V2 = nanmean(V2,1);

        fig = figure(320+pltind);
        fig.Position(3:4) = [900 500];
        clf;

        subplot(1,2,1)
        hold on
        plot(lims1{valind,pltind},lims1{valind,pltind},'k--')
        scatter(V1,V2,14,cmlines(1,:),'filled')
        xlim(lims1{valind,pltind})
        ylim(lims1{valind,pltind})
        grid on
        xlabel(sprintf('%s %s',V1ttl,valunits{valind}))
        ylabel(sprintf('%s %s',V2ttl,valunits{valind}))

        subplot(1,2,2)
        hold on
        edges = lims2{valind,pltind}(1):diff(lims2{valind,pltind})/30:lims2{valind,pltind}(2);
        N = histc(V1-V2,edges);
        b = bar(edges,N,'histc');
        b.FaceColor = cmlines(1,:);
        grid on
        xlim(lims2{valind,pltind})
        Nchans = length(find(~isnan(V1-V2)));
        M = nanmean(V1-V2);
        S = nanstd(V1-V2);
        title({...
            sprintf('%s minus %s',V1ttl,V2ttl),...
            sprintf('M: %0.4f | S: %0.4f | N: %03i | EOM: %0.4f',M,S,Nchans,S./sqrt(Nchans))...
            });

        fname = sprintf('consistplot_%s_2022_vs_%s.png',valnames{valind},pltnames{pltind});
        %saveas(fig,fullfile(figdir,fname))

    end
end
fd = fd0;

%% Type 9 schedules

% Load type 9 data
load('z:/dev/rps/sch_type9.mat')

clc
% Collate each schedule into phi_pair
unqsch = unique(fd_type9.schnum);
[phi_pair_type9, poleff_pair_type9, nrots, dks, times] = deal(NaN(length(unqsch),1));
for schind = 1:length(unqsch)
    inda = find(fd_type9.schnum==unqsch(schind) & fd_type9.ch==696);
    indb = find(fd_type9.schnum==unqsch(schind) & fd_type9.ch==697);

    if ~isempty(inda) & ~isempty(indb)
        [phi_pair_type9(schind), poleff_pair_type9(schind)] = calc_pair_diff_pol(fd_type9.phi(inda),fd_type9.phi(indb),fd_type9.xpol(inda),fd_type9.xpol(indb));
        nrots(schind) = fd_type9.nrots(inda);
        dks(schind) = fd_type9.dk_cen(inda);
        times(schind) = fd_type9.time(inda);
    end

end


%% We need some extra stuff for the next section
% Grab time, obs number, and phi_pair in the structure

% Tilt info. This takes a while... like 10 minutes.
load('z:/dev/rps/rps_tilt_data_2022.mat')

%%

fd = get_tilt_params(fd,lj_data);
fd = get_all_other_params(fd);

%%
fd_type9 = get_tilt_params(fd_type9,lj_data);
fd_type9 = get_all_other_params(fd_type9);

%% phi vs. other stuff
% Need to run the above before we run this.

%%
%fd = get_pair_params(fd,ind0,ind90);
%fd_type9 = get_pair_params(fd_type9,ind0,ind90);


if ~exist('yvxfig','var') | ~ishandle(yvxfig)
    count = 1;
    while ishandle(count)
        count = count+1;
    end
end
yvxfig = figure(count);
yvxfig.Position(3:4) = [900 400];
clf


fd.az_cen = wrapTo360(fd.az_cen);
fd.tod = fd.t-floor(fd.t);
fd_type9.az_cen = wrapTo360(fd_type9.az_cen);
fd_type9.tod = fd_type9.t-floor(fd_type9.t);

clc
medsubnames = {'','_medsub'};
for medind = 1:length(medsubnames)


    Y = {'phi', 'phi','phi_pair',...
        'xpol','xpol','poleff',...
        'n1','n1','n1_pair',...
        'n2','n2','n2_pair',...
        'Amp','Amp','amp_pair'};

    yidx = {ismember(fd.ch,ind0), ismember(fd.ch,ind90), ismember(fd.ch,ind0),...
        ismember(fd.ch,ind0), ismember(fd.ch,ind90), ismember(fd.ch,ind0),...
        ismember(fd.ch,ind0), ismember(fd.ch,ind90), ismember(fd.ch,ind0),...
        ismember(fd.ch,ind0), ismember(fd.ch,ind90), ismember(fd.ch,ind0),...
        ismember(fd.ch,ind0), ismember(fd.ch,ind90), ismember(fd.ch,ind0)};

    yidx_type9 = {ismember(fd_type9.ch,ind0), ismember(fd_type9.ch,ind90), ismember(fd_type9.ch,ind0),...
        ismember(fd_type9.ch,ind0), ismember(fd_type9.ch,ind90), ismember(fd_type9.ch,ind0),...
        ismember(fd_type9.ch,ind0), ismember(fd_type9.ch,ind90), ismember(fd_type9.ch,ind0),...
        ismember(fd_type9.ch,ind0), ismember(fd_type9.ch,ind90), ismember(fd_type9.ch,ind0),...
        ismember(fd_type9.ch,ind0), ismember(fd_type9.ch,ind90), ismember(fd_type9.ch,ind0)};

    ynames = {'phi_0','phi_90','phi_p',...
        'xpol_0','xpol_90','xpol_p',...
        'n1_0','n1_90','n1_p',...
        'n2_0','n2_90','n2_p',...
        'amp_0','amp_90','amp_p',...
        };

    Ylabs = {'\phi_0 [Degrees]','\phi_{90} [Degrees]','\phi_{pair} [Degrees]',...
        '\epsilon_0','\epsilon_{90}','1-Pol Eff.',...
        'n1_0','n1_{90}','n1',...
        'n2_0','n2_{90}','n2',...
        'A_0','A_{90}','A',...
        };

    %yttls = {'Pol 0','Pol 90','Pair-Diff'};
    if medind == 1
        ylims = {[-4.5 0], [-4.5 0]+90, [-4.5 0],...
            [-1 1]*0.02,[-1 1]*0.02, [-1 10]*1e-4,...
            [-1 1]*0.07,[-1 1]*0.07,[-1 1]*0.07,...
            [-1 1]*0.07,[-1 1]*0.07,[-1 1]*0.07,...
            [0 1]*250,[0 1]*250,[0 1]*250,...
            };
    else
        ylims = {[-1 1]*1, [-1 1]*1, [-1 1]*0.5,...
            [-1 1]*0.01,[-1 1]*0.01, [-1 1]*6e-4,...
            [-1 1]*0.07,[-1 1]*0.07,[-1 1]*0.07...
            [-1 1]*0.07,[-1 1]*0.07,[-1 1]*0.07...
            [-1 1]*30,[-1 1]*30,[-1 1]*30,...
            };
    end

    X = {'t','obsnum','az_cen','el_cen','dk_cen','az_cen_sun','el_cen_sun',...
        'x','y','xm','ym','r','theta','thetam','tod',...
        'tilt_out','tilt_out_std','tilt_temp','peak_diff'};
    Xlabs = {'Time [dd-mmm]','Observation Number','Azimuth [Degrees]','Elevation [Degrees]',...
        'Deck [Degrees]','Sun Azimuth [Degrees]','Sun Elevation [Degrees]',...
        'Inst.-Fixed X [Degrees]','Inst.-Fixed Y [Degrees]',...
        'Mirror-Fixed X [Degrees]','Mirror-Fixed Y [Degrees]','Inst.-Fixed r [Degrees]',...
        'Inst.-Fixed \theta [Degrees]','Mirror-Fixed \theta [Degrees]','Time-of-Day [Days]',...
        'Tilt Out [Degrees]','Tilt Out STD [Degrees]','Tilt Temp [deg-C]','Mod Curve Peak-to-Peak Diff'};


    for yind = 9:3:15%1:size(Y,2)

        % type 5 medsub
        if medind == 2
            [medvals] = deal(NaN(1,2640));
            for chind = 1:length(medvals)
                ind = fd.ch==chind;
                medvals(chind) = nanmedian(fd.(Y{yind})(ind));
            end
        else
            medvals = zeros(1,2640);
        end

        % Type 9 medsub
        if medind == 2
            [medvals_t9] = deal(NaN(1,2640));
            for chind = 1:length(medvals_t9)
                ind = fd_type9.ch==chind;
                medvals_t9(chind) = nanmedian(fd_type9.(Y{yind})(ind));
            end
        else
            medvals_t9 = zeros(1,2640);
        end


        for xind = 1:length(X)
            clf; hold on;
            if xind == 1
                x = mjd2datenum(fd.(X{xind})(yidx{yind}));
            else
                x = fd.(X{xind})(yidx{yind});
            end
            y = fd.(Y{yind})(yidx{yind})-medvals(fd.ch(yidx{yind}));
            scatter(x,y,14,fd.obsnum(yidx{yind})','filled')

            % Type 9 plot
            if xind == 1
                x = mjd2datenum(fd_type9.(X{xind})(yidx_type9{yind}));
                y = fd_type9.(Y{yind})(yidx_type9{yind})-medvals(fd_type9.ch(yidx_type9{yind}));
                plot(x,y,'x','Color',cmlines(2,:),'MarkerSize',14,'LineWidth',2)
            elseif xind ~= 2
                x = fd_type9.(X{xind})(yidx_type9{yind});
                y = fd_type9.(Y{yind})(yidx_type9{yind})-medvals(fd_type9.ch(yidx_type9{yind}));
                plot(x,y,'x','Color',cmlines(2,:),'MarkerSize',14,'LineWidth',2)
            end

            ylim(ylims{yind})

            colormap(cm)
            grid on
            xlabel(Xlabs{xind})
            if xind==1
                datetick('x','dd-mmm','keeplimits')
            end
            legend({'Type 5 Standard','Type 9 (Ch: 696/697 only)'})
            ylabel(Ylabs{yind})
            fname = sprintf('scatter_%s_vs_%s%s.png',ynames{yind},X{xind},medsubnames{medind});
            saveas(yvxfig,fullfile(figdir,fname))
        end

    end

end

% From these plots we discovered two things:
% -- Anticorrelation between A/B angles
% -- Anticorrelation/Correlation between N1/N2 for 0/90 dets
%% Plot psi0 vs psi90
clc

fig = figure(1);
fig.Position(3:4) = [600 500];
clf; hold on;

V1 = atand(tand(phis(:,ind0)-nanmedian(phis(:,ind0),1)));
V2 = atand(tand(phis(:,ind90)-nanmedian(phis(:,ind90),1)));
legttls = titles;
%ttls{end+1} = 'Slope=-1';
lims = [-1 1]*0.6;
plot(lims,-lims,'--','Color',[1 1 1]*0.75,'LineWidth',1)
clear z
for obsind = 1:length(dks)
    idx = V1(obsind,:)~=0 & V2(obsind,:)~=0;
    z(obsind) = plot(V1(obsind,idx)',V2(obsind,idx)','.','MarkerSize',14,'Color',cmlines(1,:));
end
xlim(lims)
ylim(lims)
grid on
xlabel({'\phi_{d,0}-Md(\phi_{d,0}) [Degrees]',''})
ylabel({'','\phi_{d,90}-Md(\phi_{d,90}) [Degrees]'})
%leg = legend(z,legttls);
%title(leg,'DK''s:')
title('Pol 90 Vs Pol 0')
pbaspect([1 1 1])

fname = 'pol90_vs_pol0_type5.png';
saveas(fig,fullfile(figdir,fname),'png')

%% Now plot the covariance vs DK
[C D] = deal([]);
for obsind = 1:length(dks)
    idx = V1(obsind,:)~=0 & V2(obsind,:)~=0;
    C0 = nancov(V1(obsind,:),V2(obsind,:));
    D0 = corrcov(C0);
    C(end+1) = C0(1,2);
    D(end+1) = D0(1,2);
end

clear z
fig = figure(2);
fig.Position(3:4) = [600 300];
clf; hold on;
[s si] = sort(dks);
plot(dks(si),C(si),'Color',cmlines(1,:))
z(1) = plot(dks(si),C(si),'.','Color',cmlines(1,:),'MarkerSize',14);
grid on
xlabel('DK [Degrees]')
ylabel('\phi_d Covariance [Deg^2]')

% What about type 9's?
clc
unqsch = unique(fd_type9.schnum);
[phis9] = deal(NaN(length(unqsch),2)); % row: obs, col: 0/90
[dks9] = deal(NaN(size(unqsch)));
for obsind = 1:length(unqsch)
    schind = find(fd_type9.schnum==unqsch(obsind));
    if length(schind)==2
        dk0 = mean(fd_type9.dk_cen(schind));
        if ~isnan(dk0)
            dks9(obsind) = dk0;
            chind = fd_type9.ch(schind)==696; % 696 is phi0
            phis9(obsind,1) = fd_type9.phi_medsub(schind(chind));
            phis9(obsind,2) = fd_type9.phi_medsub(schind(~chind));
            dks9(obsind) = mean(fd_type9.dk_cen(schind));
        end
    end
end

idx = inrange(dks9,-10,10);
fig = figure(3);
fig.Position(3:4) = [600 500];
clf; hold on;
plot(phis9(idx,1),phis9(idx,2),'.')
lims = [-1 1]*0.6;
plot(lims,-lims,'--','Color',[1 1 1]*0.75,'LineWidth',1)
xlim(lims)
ylim(lims)
colormap(cm)
grid on
title('Type 9 Channel 696 vs. 676 ')
xlabel({'\phi_{d,0}-Md(\phi_{d,0}) [Degrees]','Ch 696'})
ylabel({'Ch 697','\phi_{d,90}-Md(\phi_{d,90}) [Degrees]'})
pbaspect([1 1 1])
fname = 'pol90_vs_pol0_type9.png';
saveas(fig,fullfile(figdir,fname),'png')


C = cov(phis9(idx,1),phis9(idx,2));
fig = figure(2);
z(2) = plot(0,C(1,2),'x','Color',cmlines(2,:));
grid on
ylim([-0.02 0])
legend(z,{'Standard Type 5','Type 9 (696/697 only)'},'Location','southeast')
fname = 'cov_vs_dk.png';
title('A/B Covariance vs. DK')
saveas(fig,fullfile(figdir,fname),'png')



clear z
fig = figure(4);
fig.Position(3:4) = [600 300];
clf; hold on;
[s si] = sort(dks);
plot(dks(si),D(si),'Color',cmlines(1,:))
z(1) = plot(dks(si),D(si),'.','Color',cmlines(1,:),'MarkerSize',14);
grid on
xlabel('DK [Degrees]')
ylabel('\phi_d Correlation')

% What about type 9's?
clc
unqsch = unique(fd_type9.schnum);
[phis9] = deal(NaN(length(unqsch),2)); % row: obs, col: 0/90
[dks9] = deal(NaN(size(unqsch)));
for obsind = 1:length(unqsch)
    schind = find(fd_type9.schnum==unqsch(obsind));
    if length(schind)==2
        dk0 = mean(fd_type9.dk_cen(schind));
        if ~isnan(dk0)
            dks9(obsind) = dk0;
            chind = fd_type9.ch(schind)==696; % 696 is phi0
            phis9(obsind,1) = fd_type9.phi_medsub(schind(chind));
            phis9(obsind,2) = fd_type9.phi_medsub(schind(~chind));
            dks9(obsind) = mean(fd_type9.dk_cen(schind));
        end
    end
end

idx = inrange(dks9,-10,10);
C = cov(phis9(idx,1),phis9(idx,2));
D = corrcov(C);
fig = figure(4);
z(2) = plot(0,D(1,2),'x','Color',cmlines(2,:));
grid on
ylim([-1.1 0.1])
legend(z,{'Standard Type 5','Type 9 (696/697 only)'},'Location','southeast')
fname = 'corr_vs_dk.png';
title('A/B Correlation vs. DK')
saveas(fig,fullfile(figdir,fname),'png')


%% Create sims for correlated amp-scaling noise

% Grab the mod curves for A dets
idx = find(ismember(fd.ch,ind0) & fd.obsnum == 7);

[B, B90, R] = deal([]);
fig = figure(2);
clf; hold on;

isused = false(size(idx));
for idxind = 1:length(idx)
    fdi = idx(idxind);
    fdi90 = fd.pair_idx(fdi);
    if fd.nrots(fdi) == 13 & ~isnan(fdi90)
        b = fd.bparam{fdi}(6:end)./fd.Amp(fdi);
        b90 = fd.bparam{fdi90}(6:end)./fd.Amp(fdi90);
        if ~isnan(sum(b+b90))
        B(end+1,:) = b;
        B90(end+1,:) = b90;
        R(end+1,:) = fd.rot{fdi}; 
        plot(R(end,:),B(end,:))
        isused(idxind) = true;
        end
    end
end

grid on


%% Fits to mock data with correlated amp-scaling noise.

clc
mxfev = 100000;
mxiter = 100000;
options = optimset('TolFun',1e-10,'MaxIter',mxiter,'MaxFunEvals',mxfev,'Display','off');

Niter = 2000;
sigN1 = 0.01*0;
sigN2 = 0.01*0;
sigA = 0.1*0;
sigB = 0.1*0;
sigAmp = 0.013;
sigE = 0.001*0;
dk = 0;

tic;
clear Pa Pb Ras
fig = figure(3);
clf; hold on;
for iterind = 1:Niter
    rot0 = (-180:15:180)+dk;
    rot = rot0+0*normrnd(0,0.001,size(rot0));
    phia = normrnd(0,sigA,1,1);
    phib = normrnd(90,sigB,1,1);
    N1 = normrnd(0,sigN1,1,1);
    N2 = normrnd(0,sigN2,1,1);
    epsa = normrnd(0.01,sigE,1,1);
    epsb = normrnd(0.01,sigE,1,1);
    Anoise = normrnd(0,sigAmp,size(rot));

    parm = [phia epsa N1 N2 1];
    Ra = rps_get_mod_model(parm,rot);
    Ra = Ra+(Ra).*Anoise;
    Ras(iterind,:) = Ra;
    parm = [phib epsb N1 N2 1];
    Rb = rps_get_mod_model(parm,rot);
    Rb = Rb+(Rb).*Anoise;

    lb = [-10 -0.5 -10 -10 0];
    ub = [169 0.5 10 10 1e3];
    guess = [0 0 0 0 max(Ra)/2];
    parm_a = lsqcurvefit(@rps_get_mod_model,guess,rot0,Ra,lb,ub,options);
    Pa(iterind) = parm_a(1);

    %lb = [-20+90 -0.5 -10 -10 0];
    %ub = [20+90 0.5 10 10 1e3];
    guess = [90 0 0 0 max(Rb)/2];
    parm_b = lsqcurvefit(@rps_get_mod_model,guess,rot0,Rb,lb,ub,options);
    Pb(iterind) = parm_b(1);

    plot(rot0,Ra)

end
toc;

%
fig = figure(1);
fig.Position(3:4) = [600 500];
clf; hold on;

lims = [-1 1]*0.6;
plot(lims,-lims,'--','Color',[1 1 1]*0.75,'LineWidth',1)
plot(Pa,Pb-90,'.','MarkerSize',14,'Color',cmlines(1,:));
xlim(lims)
ylim(lims)
grid on
pbaspect([1 1 1])
xlabel({'\phi_{d,0} [Degrees]',''})
ylabel({'','\phi_{d,90}-90^\circ [Degrees]'})
title(sprintf('\\phi 0 Vs. 90 of N=%i noise-correlated sims',Niter))
C = (cov(Pa,Pb));
D = corrcov(C);
text(0.1,0.4,{sprintf('0/90 Cov: %1.3f [Deg^2]',C(1,2))...
    sprintf('0/90 Corr: %1.3f',D(1,2))})

fname = 'pol90_vs_pol0_sim.png';
saveas(fig,fullfile(figdir,fname),'png')


%% Grab the mod curves for A/B dets 

[Amodel Adata] = deal({});
fig = figure(1);
fig.Position(3:4) = [600 500];
clf; hold on;

for obsind = 8:9%1:length(dks)
idx = find(ismember(fd.ch,ind0) & fd.obsnum == obsind);

[B, B90, R, Bmod,Bmod90] = deal([]);
% fig = figure(2);
% clf; hold on;

isused = false(size(idx));
for idxind = 1:length(idx)
    fdi = idx(idxind);
    fdi90 = fd.pair_idx(fdi);
    if fd.nrots(fdi) == 13 && ~isnan(fdi90)
        b = fd.bparam{fdi}(6:end)./fd.Amp(fdi);
        b90 = fd.bparam{fdi90}(6:end)./fd.Amp(fdi90);
        if ~isnan(sum(b+b90))
        B(end+1,:) = b;
        B90(end+1,:) = b90;
        R(end+1,:) = fd.rot{fdi}; 
        Bmod(end+1,:) = rps_get_mod_model([fd.phi(fdi),fd.xpol(fdi),fd.n1(fdi),fd.n2(fdi),fd.Amp(fdi)],fd.rot{fdi}+fd.phi_s(fdi))./fd.Amp(fdi);
        Bmod90(end+1,:) = rps_get_mod_model([fd.phi(fdi90),fd.xpol(fdi90),fd.n1(fdi90),fd.n2(fdi90),fd.Amp(fdi90)],fd.rot{fdi90}+fd.phi_s(fdi90))./fd.Amp(fdi90);
        %plot(R(end,:),B(end,:))
        %grid on
        isused(idxind) = true;
        end
    end
end

% Plot fracdiff residuals A vs. B per zeta
Bres = (Bmod-B)./Bmod;
Bres90 = (Bmod90-B90)./Bmod90;

idx = [2 3 5 6 8 9 11 12];
idx = [1 4 7 10];
%idx = 1:13;
lims = [-1 1]*0.06;
plot(lims,lims,'--','Color',[1 1 1]*0.75,'LineWidth',1)
plot(Bres(:,idx),Bres90(:,idx),'.','MarkerSize',14,'Color',cmlines(1,:));
xlim(lims)
ylim(lims)
grid on
end

pbaspect([1 1 1])
xlabel({'Pol0',''})
ylabel({'','Pol90'})
title('Mod Curve fractional residuals')

fname = 'pol90_vs_pol0_mod_curve_res.png';
saveas(fig,fullfile(figdir,fname),'png')


%% Plot modcurve residuals!?

if ~exist('modcurvefig','var') | ~ishandle(modcurvefig)
    count = 1;
    while ishandle(count)
        count = count+1;
    end
end
modcurvefig = figure(count);
clf; hold on;
t = tiledlayout(2,1);
t.TileSpacing = 'tight';
t.Padding = 'tight';

clc
tic;
for obsind = 1:length(dks)
    clf;
    obschans = find(fd.obsnum==obsind);
    for chidx = 1:length(obschans)
        chind = obschans(chidx);
        if inrange(fd.phi(chind),-10,10)
            nexttile(1)
            hold on
        else
            nexttile(2)
            hold on
        end

        Adata = fd.bparam{chind}(6:end);
        Amodel = rps_get_mod_model([fd.phi(chind),fd.xpol(chind),fd.n1(chind),fd.n2(chind),fd.Amp(chind)],fd.rot{chind}+fd.phi_s(chind));

        idx = ~isnan(fd.rot{chind}) & ~isnan(Adata) & ~isnan(Amodel);
        plot(fd.rot{chind}(idx),(Adata(idx)-Amodel(idx))./fd.Amp(chind))

    end

    nexttile(1)
    grid on;
    ylabel({'Pol 0','Fractional Residual'})
    ylim([-1 1]*0.1)
    nexttile(2)
    grid on
    ylabel({'Pol 90','Fractional Residual'})
    xlabel('\zeta_{grav} [Degrees]')
    ylim([-1 1]*0.1)

    fname = sprintf('modcurveres_dk%s',titles{obsind});
    saveas(modcurvefig,fullfile(figdir,fname),'png')
end
toc


%% Angle jacks ... drop this for now.
% Gotta come up with a better method for this.

mxfev = 100000;
mxiter = 100000;
options = optimset('TolFun',1e-10,'MaxIter',mxiter,'MaxFunEvals',mxfev,'Display','off');

clc
[jack1, jack2] = deal(NaN(length(fd.ch),5));
tic;
for chind = 1:length(fd.ch)
    %if ~isnan(fd.obsnum(chind)) & ismember(dks(fd.obsnum(chind)),[0,90])
    ch = fd.ch(chind);
    A = fd.bparam{chind}(6:end);
    rot_adj = fd.rot{chind}+fd.phi_s(chind);
    %rot_adj = fd.rot{chind}-dks(fd.obsnum(chind))+90;
    %angguess = atand(tand(p.chi(ch)+p.chi_thetaref(ch)));
    parm0 = [fd.phi(chind),fd.xpol(chind), fd.n1(chind), fd.n2(chind),fd.Amp(chind)];
    lb = [-20 -1e-6 -1e-6 -1e-6 -1e-5] + parm0;
    ub = [20 1e-6 1e-6 1e-6 1e-5]+parm0;
    %lb = [-20 -0.5 -10 -10 -10]+parm0;
    %ub = [20 0.5 10 10 10]+parm0;
    guess = parm0;

    if fd.nrots(chind) == 13
        idx = sign(cosd(rot_adj+45))>0;
        jack1(chind,:) = lsqcurvefit(@rps_get_mod_model,guess,rot_adj(idx),A(idx),lb,ub,options);
        jack2(chind,:) = lsqcurvefit(@rps_get_mod_model,guess,rot_adj(~idx),A(~idx),lb,ub,options);
    end
    %end
end
toc
scatter(1:length(fd.ch),jack1(:,1)-jack2(:,1),12,fd.dk_cen,'filled')

%% Angle jacks with type 9's ... drop this for now.
% Gotta come up with a better method for this.

mxfev = 100000;
mxiter = 100000;
options = optimset('TolFun',1e-10,'MaxIter',mxiter,'MaxFunEvals',mxfev,'Display','off');

clc
[jack1, jack2] = deal(NaN(length(fd_type9.ch),5));
phi = [];
tic;
for chind = 1:length(fd_type9.ch)
    %if ~isnan(fd_type9.obsnum(chind)) & ismember(dks(fd_type9.obsnum(chind)),[0,90])
    ch = fd_type9.ch(chind);
    A = fd_type9.bparam{chind}(6:end);
    rot_adj = fd_type9.rot{chind}+fd_type9.phi_s(chind);
    %rot_adj = fd_type9.rot{chind}-dks(fd_type9.obsnum(chind))+90;
    %angguess = atand(tand(p.chi(ch)+p.chi_thetaref(ch)));
    parm0 = [fd_type9.phi(chind),fd_type9.xpol(chind), fd_type9.n1(chind), fd_type9.n2(chind),fd_type9.Amp(chind)];
    lb = [-20 -1e-6 -1e-6 -1e-6 -1e-5] + parm0;
    ub = [20 1e-6 1e-6 1e-6 1e-5]+parm0;
    %lb = [-20 -0.5 -10 -10 -10]+parm0;
    %ub = [20 0.5 10 10 10]+parm0;
    guess = parm0;
    phi(end+1) = fd_type9.phi(chind);
    if fd_type9.nrots(chind) == 37

        idx = sign(cosd(rot_adj+45))>0;
        idx2 = ~idx;
        %idx = 2:3:length(rot_adj);
        %idx2 = 3:3:length(rot_adj);
        jack1(chind,:) = lsqcurvefit(@rps_get_mod_model,guess,rot_adj(idx),A(idx),lb,ub,options);
        jack2(chind,:) = lsqcurvefit(@rps_get_mod_model,guess,rot_adj(idx2),A(idx2),lb,ub,options);
    end
    %end
end
toc

clf; hold on;
%scatter(1:length(fd_type9.ch),jack1(:,1)-jack2(:,1),12,fd_type9.dk_cen,'filled')
scatter(1:length(fd_type9.ch),phi'-jack1(:,1),12,'filled','Color',cmlines(1,:))
scatter(1:length(fd_type9.ch),phi'-jack2(:,1),12,'filled','Color',cmlines(2,:))
%idx = fd_type9.nrots == 37;
%scatter(1:length(fd_type9.ch)(idx),fd_type9.phi(idx),12,fd_type9.dk_cen,'filled')


%% Calculate phi_s_prime

clc
source = struct;
source.distance = 195.5;

mirror = struct;
mirror.height = 1.4592;
p0 = rmfield(p,'expt');

[fd.xp,fd.yp,fd.phi_sp] = deal(NaN(size(fd.ch)));
for chind = 1:length(fd.ch)
    source.azimuth = fd.src_az(chind);
    source.height = source.distance*tand(fd.src_el(chind));
    mirror.tilt = fd.mirr_tilt(chind);
    mirror.roll = fd.mirr_roll(chind);
    [fd.xp(chind),fd.yp(chind),fd.phi_sp(chind)] = beam_map_pointing_model(fd.az_cen(chind),fd.el_cen(chind),fd.dk_cen(chind),...
        model,'bicep3',mirror,source,structcut(p0,fd.ch(chind)));
end

%% Cut stuff: Plot B18 rgls and dets that pass RPS cuts

A = NaN(size(p.gcp));
A(p_ind.rgl) = 0;
A(~isnan(nanmean(phi_pair,1))) = 1;

fig = figure(13453);
fig.Position(3:4) = [900 800];
ttl = 'B18 RGL pairs that pass RPS cuts';
plot_tiles(A,p,'pair','sum','fig',13453,'title',ttl);
fname = 'B18_and_RPS_cut_plot';
saveas(fig,fullfile(figdir,fname),'png')

%% Med sub parameters vs other parameters plots

V = {phis,xpols,n1s,n2s,amps};
Vttl = {'\phi','\epsilon','N1','N2','A'};
lims = {0.4, 0.01, 0.05, 0.05, 30};
pltttl = {'0','90'};
idx = {ind0, ind90};

clc
for pltind = 1:length(idx)
    fig = figure(5093+pltind);
    fig.Position(3:4) = [800 800];
    clf;
    t = tiledlayout(5,5);
    t.TileSpacing = 'compact';
    t.Padding = 'tight';
    rc2idx = reshape(1:25,5,5);
    for xind = 1:(length(V)-0)
        for yind = (xind+0):length(V)
            nexttile(rc2idx(xind,yind))
            x = V{xind}(:,idx{pltind});
            x = x-nanmedian(x,1);
            y = V{yind}(:,idx{pltind});
            y = y-nanmedian(y,1);
            obs = repmat((1:10),1,size(y,2));
            ax = gca;
            if xind == yind
                edges = (-1:0.075:1)*lims{xind};
                hold on;

                N = histc(reshape(x,[],1),edges);
                b = bar(edges,N,'histc');

                title(Vttl{xind})
            else
                scatter(reshape(x,[],1),reshape(y,[],1),14,reshape(obs,[],1),'filled')
                colormap(cm)

                ylim([-1 1]*lims{yind})
                if xind == 1
                    ylabel(Vttl{yind})
                else
                    ax.YTickLabels = {};
                end

            end
            

            if yind == 5
                xlabel(Vttl{xind})
            else
                ax.XTickLabels = {};
            end
            xlim([-1 1]*lims{xind})

            grid on
            pbaspect([1 1 1])
        end
    end
    sgtitle({'Pol Parameter Corner Plots',sprintf('\\phi_{%s} Detectors',pltttl{pltind})})
    fname = sprintf('pol_param_corner_plots_%s',pltttl{pltind});
    saveas(fig,fullfile(figdir,fname),'png')
end

%% Look at type 9 no homing

starttime = 59606+7/24+05/24/60;
endtime = 59606+11/24+30/24/60;

fig = figure(40958);
clf;
idx = find(inrange(fd_type9.t,starttime,endtime));
std(fd_type9.phi_pair(idx(1:2:end)))


%% Check N chans used per year vs how many we have from the RPS.
% Run on odyssey.
Nchans = NaN(3,2);
yrs = {'2016','2017','2018'};
k = ParameterRead('aux_data/chi/phi_dummy_bicep3_20150101.csv');
usedchans = zeros(2,2640);

for yrind = 1:3
% load the channel flags used by B2018
clear chflags
chflags = get_default_chflags([],yrs{yrind});

clear get_array_info
[p1, p1_ind] = get_array_info(20160505+10000*(yrind-1),'obs','obs','obs','obs',chflags);

flags = struct();
flags.filebase{1} = 'chi/phi_dummy';
flags.par_name = {'hasmeas'};
flags.low = [0.5];
flags.high = [1.5];

chflags = structcat(chflags,flags);

clear get_array_info
[p2, p2_ind] = get_array_info(20160505+10000*(yrind-1),'obs','obs','obs','obs',chflags);

Nchans(yrind,:) = [length(p1_ind.rgl100a), length(p2_ind.rgl100a)];
usedchans(1,p1_ind.rgl100) = 1;
usedchans(2,p2_ind.rgl100) = 1;

end



if 1
    clc
    fprintf('\tNchans\t\n')
    fprintf('Year\tCMB\tRPS\t%%\n')
    for yrind = 1:3

        fprintf('%s\t%03i\t%03i\t%1.0f\n',yrs{yrind},Nchans(yrind,1),Nchans(yrind,2),Nchans(yrind,2)./Nchans(yrind,1)*100)
    end


    fprintf('\n\nTotal unique pairs used\n')
    fprintf('CMB\tRPS\t%%\n')
    fprintf('%i\t%i\t%1.0f\n',sum(usedchans(1,:),2)/2,sum(usedchans(2,:),2)/2,(sum(usedchans(2,:),2)./sum(usedchans(1,:),2))*100)
end


%% End of main function
function [fd, phis, phi_pair, xpol, poleffs,n1s,n2s,amps] = get_pol_params_per_obs(fd,p,obscell)

if ~exist('obscell','var')
    obscell = num2cell(unique(fd.schnum));
end


len = length(obscell);
fd.obsnum = NaN(size(fd.ch));
for chind = 1:length(fd.ch)
    % Obs Number
    for obsind = 1:len
        if ismember(fd.schnum(chind),obscell{obsind})
            fd.obsnum(chind) = obsind;
        end
    end
end

[phis, phi_pair, xpol, poleffs,n1s,n2s,amps] = deal(NaN(len,length(p.gcp)));
for obsind = 1:len
    for chind = 1:length(p.gcp)
        idx = find(fd.ch==chind & fd.obsnum==obsind);
        if any(idx)
            phis(obsind,chind) = nanmean(fd.phi(idx));
            xpol(obsind,chind) = nanmean(fd.xpol(idx));
            n1s(obsind,chind) = nanmean(fd.n1(idx));
            n2s(obsind,chind) = nanmean(fd.n2(idx));
            amps(obsind,chind) = nanmean(fd.Amp(idx));
            if any(fd.pair_idx(idx))
                phi_pair(obsind,chind) = nanmean(fd.phi_pair(idx));
                poleffs(obsind,chind) = nanmean(fd.poleff(idx));
            end
        end
    end
end


function fd = get_pair_params(fd,ind0,ind90)

checkflds = {'n1','n2','Amp'};
for fldind = 1:length(checkflds)
    if ~isfield(fd,checkflds{fldind})
        fd.(checkflds{fldind}) = NaN(size(fd.ch));
    end
end

[fd.phi_pair, fd.poleff, fd.pair_idx] = deal(NaN(size(fd.ch)));
for chind = 1:length(fd.ch)
    % Phi/poleff pair
    isind0 = find(ismember(ind0,fd.ch(chind)));
    if ~isempty(isind0)
        has90 = find(fd.schnum == fd.schnum(chind) & fd.rowind==fd.rowind(chind) & fd.ch == ind90(isind0));
        if ~isempty(has90)
            fd.pair_idx(chind) = has90(1);
            [fd.phi_pair(chind) fd.poleff(chind)] = calc_pair_diff_pol(fd.phi(chind),...
                nanmean(fd.phi(has90)),fd.xpol(chind),nanmean(fd.xpol(has90)));
        end
    end

idx = find(~isnan(fd.pair_idx));
[fd.n1_pair, fd.n2_pair,fd.amp_pair] = deal(NaN(size(fd.ch))); 
fd.n1_pair(idx) = (fd.n1(fd.pair_idx(idx))+fd.n1(idx))/2;
fd.n2_pair(idx) = (fd.n2(fd.pair_idx(idx))+fd.n2(idx))/2;
fd.amp_pair(idx) = (fd.Amp(fd.pair_idx(idx))+fd.Amp(idx))/2;

end

function fd = get_tilt_params(fd,lj_data)

for chind = 1:length(fd.ch)
    % Tilt info
    if 1
        tiltcal = [0.284157 -0.0118];
        ind = lj_data.time >= fd.t(chind) & lj_data.time<=fd.t2(chind);
        fd.tilt_out(chind) = nanmean(polyval(tiltcal,lj_data.AIN0(ind)));
        fd.tilt_out_std(chind) = nanstd(polyval(tiltcal,lj_data.AIN0(ind)));
        fd.tilt_temp(chind) = nanmean(lj_data.AIN2(ind)*100);
    end
end
%fd.tilt_out = interp1(lj_data.time',lj_data.AIN0',mean([fd.t1;fd.t2],1));
%fd.tilt_out = polyval(tiltcal,fd.tilt_out);



function fd = get_all_other_params(fd)

% Peak difference in residuals
fd.peak_diff = NaN(size(fd.ch));
for chind = 1:length(fd.ch)

    B = fd.bparam{chind}(6:end);
    R = fd.rot{chind};
    [m midx] = nanmax(B);
    Rmax = R(midx);
    signidx = 1;
    if midx > length(B)/2
        signidx = -1;
    end

    midx2 = find(inrange(R,Rmax+signidx*180-5,Rmax+signidx*180+5));
    if ~isempty(midx2)
        peak_diff = (B(min(midx,midx2))-B(max(midx,midx2)))/max(B([midx, midx2]));
        if abs(peak_diff)<0.1
            fd.peak_diff(chind) = peak_diff;
        end
    end
end

% Mirror coords
[fd.theta, fd.r] = cart2pol(fd.x,fd.y);
fd.theta = wrapTo360(fd.theta*180/pi);
fd.thetam = wrapTo360(fd.theta-fd.dk_cen);
[fd.xm, fd.ym] = pol2cart(fd.thetam*pi/180,fd.r);

% Grab the moon-sun angles
% MJD, , ,raapp,decapp
fname = fullfile('z:/dev/sun_check_2022Aug12.txt');
f = fopen(fname);
hd_sun = textscan(f,'%f%s%s%f%f','delimiter',',','HeaderLines',60);
hd_sun{1} = hd_sun{1}-2400000.5;
fclose(f);
%fd.az_cen_sun = interp1(hd_sun{1},hd_sun{4},fd.t_cen);
fd.az_cen_sun = wrapTo180(interp1(hd_sun{1},unwrap(hd_sun{4}*pi/180)*180/pi,fd.t));
fd.el_cen_sun = interp1(hd_sun{1},hd_sun{5},fd.t);

%
clc
% MJD, , ,raapp,decapp
fname = fullfile('z:/dev/moon_check_2022Aug12.txt');
f = fopen(fname);
hd_moon = textscan(f,'%f%s%s%f%f','delimiter',',','HeaderLines',61);
hd_moon{1} = hd_moon{1}-2400000.5;
fclose(f);
%fd.az_cen_moon = interp1(hd_moon{1},hd_moon{4},fd.t_cen);
fd.az_cen_moon = interp1(hd_moon{1},unwrap(hd_moon{4}*pi/180)*180/pi,fd.t);
fd.el_cen_moon = interp1(hd_moon{1},hd_moon{5},fd.t);


