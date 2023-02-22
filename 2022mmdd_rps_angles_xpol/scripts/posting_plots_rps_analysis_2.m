function posting_plots_rps_analysis_2()

%% Run this first
addpath('z:/dev/rps/')
addpath('z:/pipeline/beammap')
addpath('z:/pipeline/util')
addpath('z:/dev/diff_polarization/')
addpath('z:/dev')
%% Then this
clear all
close all
clc

p22 = load('z:/dev/rps/fpu_data_obs.mat'); % pointing from 2022
p18 = load('z:/dev/rps/fpu_data_B18.mat'); % RGL's from B18
p = p22.p;
p_ind = p18.p_ind;

% Load no tilt data
load('z:/dev/rps/rps_beam_fits_type5_notilt_cut.mat')
fd_notilt = fd;


% Load type 9 data
%load('z:/dev/rps/type9_fit_dat_mirrorfitted_cut.mat')
load('z:/dev/rps/rps_beam_fits_type9_21feb_rerun.mat');
fd_type9 = rps_cut_fitdata(fd,p,p_ind,false);
load('z:/dev/rps/sch_type9.mat')


% Load 2018 data
load('z:/dev/rps/rps_beam_fits_cut_2018.mat')
fd.xpol = fd.aparam(:,2);
fd.phi = fd.phi_d;
fd.phi_err = fd.aerr(:,1);
fd.schnum = fd.sch;
fd_2018 = fd;


load('z:/dev/rps/rps_obs_info.mat')
%load('z:/dev/rps/rps_beam_fits_type5_withbparam.mat')
load('z:/dev/rps/rps_beam_fits_type5_21feb_rerun.mat');
fd = rps_cut_fitdata(fd,p,p_ind,0);
load('z:/pipeline/beammap/viridis_cm.mat')
%load('z:/dev/sims/coaddopt_6600.mat')
load('z:/dev/rps/sch_type5.mat')
load('z:/dev/rps/pm.mat')
figdir = fullfile('C:','Users','James','Documents','GitHub','postings','2022mmdd_rps_angles_xpol','figs','');

cmlines = colormap('lines');

[fd.theta, fd.r] = pol2cart(fd.x,fd.y);
fd.theta = fd.theta*180/pi;

%fd.phi_corr = mirror_diff_pol_calc(fd.phi,fd.r,fd.theta,2200,2200,fd.dk_cen,fd.mirr_tilt,fd.mirr_roll);
%fd.phi_corr = reshape(fd.phi_corr,size(fd.ch));

len = length(scheds);
for chind = 1:length(fd.ch)
    % Obs Number
    for schind = 1:len
        if ismember(fd.schnum(chind),scheds{schind})
            fd.obsnum(chind) = schind;
        end
    end
end

% If 1 use chflags from B18
if 1
inda = p_ind.rgla;
indb = p_ind.rglb;

ind0 = [p_ind.rgl100a(ismember(p_ind.rgl100a,find(p.mce~=0))) ...
    p_ind.rgl100b(ismember(p_ind.rgl100b,find(p.mce==0)))];
ind90 = [p_ind.rgl100b(ismember(p_ind.rgl100b,find(p.mce~=0))) ...
    p_ind.rgl100a(ismember(p_ind.rgl100a,find(p.mce==0)))];

else
inda = p_ind.a;
indb = p_ind.b;

ind0 = [p_ind.a(ismember(p_ind.a,find(p.mce~=0))) ...
    p_ind.b(ismember(p_ind.b,find(p.mce==0)))];
ind90 = [p_ind.b(ismember(p_ind.b,find(p.mce~=0))) ...
    p_ind.a(ismember(p_ind.a,find(p.mce==0)))];
end


%% Then this 
% Grab the phis per obs
clc

[phis, xpols, phis_err, xs, ys] = deal(NaN(len,2640));
for schind = 1:len
    for chind = 1:2640
        ci = find(fd.ch==chind & ismember(fd.schnum,scheds{schind}));

        if ~isempty(ci)

            %ch_chi = atand(tand(p.theta(chind)+p.chi(chind)));

            val = reshape(fd.phi(ci),[],1);
            %val = reshape(atand(tand(fd.phi(ci)-fd.phi_s(ci)+fd.dk_cen(ci)-90)),[],1);
            xp = reshape(fd.xpol(ci),[],1);
            err = reshape(fd.phi_err(ci),[],1);

            %err(isnan(val))=1e10;
            xs(schind,chind) = mean(fd.x(ci));
            ys(schind,chind) = mean(fd.y(ci));
            phis(schind,chind) = mean(val);%wmean(val,1./err,1);
            xpols(schind,chind) = mean(xp);%wmean(xp,1./err,1);
            phis_err(schind,chind) = nanmin(err);

        end
    end
end

% Calculate pair-diff angles
% Loop over channels to account for MCE0;
[phi_pair, poleff_pair, phi_pair_err,Qpair,Upair] = deal(NaN(len,2640));
for schind = 1:len
    for chind = 1:length(ind0)

        pha = phis(schind,ind0(chind));
        phb = phis(schind,ind90(chind));
        ea = xpols(schind,ind0(chind));
        eb = xpols(schind,ind90(chind));

        phi_pair_err(schind,ind0(chind)) = max(phis_err(schind,[ind0(chind) ind90(chind)]));

            [phi_pair(schind,ind0(chind)), poleff_pair(schind,ind0(chind))] = calc_pair_diff_pol(pha,phb,ea,eb);
            [Qpair(schind,ind0(chind)),Upair(schind,ind0(chind))] = calc_pair_diff_stokes(pha,phb,ea,eb);


    end
end



% Do the same thing for 2018 data

% Fix the phis so that they're at around 0 or 90
for  chind = 1:length(fd_2018.ch)
    ch_chi = atand(tand(p.theta(fd_2018.ch(chind))+p.chi(fd_2018.ch(chind))));
    if abs(ch_chi)<45
        fd_2018.phi(chind) = atand(tand(fd_2018.phi(chind)));
    else
        fd_2018.phi(chind) = atand(tand(fd_2018.phi(chind)-90))+90;
    end
end

% Make cuts based on where the median is WRT the FP orientation
cutind = true(size(fd_2018.ch));
for chind = 1:2640
    ci = find(fd_2018.ch==chind);
    if ~isempty(ci)
        A = atand(tand(p.theta(chind)+p.chi(chind)-fd_2018.phi(ci)'));
        B = atand(tand(p.theta(chind)+p.chi(chind)-nanmedian(fd_2018.phi(ci))));
        cutind(ci) = abs(A-B)<1;%0.45;
    end
end
fd_2018 = structcut(fd_2018,cutind);


%
% Order the a/b phis based on schedule
[phis_2018, xpols_2018, phis_err_2018, xs_2018, ys_2018] = deal(NaN(len,2640));
for schind = 1:len
    for chind = 1:2640
        ci = find(fd_2018.ch==chind & ismember(fd_2018.schnum,scheds{schind}));

        if ~isempty(ci)

            ch_chi = atand(tand(p.theta(chind)+p.chi(chind)));

            val = reshape(fd_2018.phi(ci),[],1);
            xp = reshape(fd_2018.xpol(ci),[],1);
            err = reshape(fd_2018.phi_err(ci),[],1);

            %err(isnan(val))=1e10;
            %xs_2018(schind,chind) = mean(fd_2018.x(ci));
            %ys_2018(schind,chind) = mean(fd_2018.y(ci));
            phis_2018(schind,chind) = mean(val);%wmean(val,1./err,1);
            xpols_2018(schind,chind) = mean(xp);%wmean(xp,1./err,1);
            phis_err_2018(schind,chind) = nanmin(err);

        end
    end
end

if 1
    % Calculate pair-diff angles
    % Loop over channels to account for MCE0;
    [phi_pair_2018, poleff_pair_2018, phi_pair_err_2018] = deal(NaN(len,2640));
    for schind = 1:len
        for chind = 1:length(ind0)

            pha = phis_2018(schind,ind0(chind));
            phb = phis_2018(schind,ind90(chind));
            ea = xpols_2018(schind,ind0(chind));
            eb = xpols_2018(schind,ind90(chind));

            phi_pair_err_2018(schind,ind0(chind)) = max(phis_err_2018(schind,[ind0(chind) ind90(chind)]));
            if p.mce(ind0(chind))~=0

                [phi_pair_2018(schind,ind0(chind)), poleff_pair_2018(schind,ind0(chind))] = calc_pair_diff_pol(pha,phb,ea,eb);
            else

                [phi_pair_2018(schind,ind0(chind)), poleff_pair_2018(schind,ind0(chind))] = calc_pair_diff_pol(phb,pha,eb,ea);
            end


        end

    end

else
    [phi_pair_2018, poleff_pair_2018] = deal(NaN(1,2640));
    for chind = 1:length(inda)
        ia = find(fd_2018.ch==inda(chind),1);
        ib = find(fd_2018.ch==indb(chind),1);
        if ~isempty(ia) & ~isempty(ib)
            pha = wmean(fd_2018.phi(ia),1./fd_2018.phi_err(ia),1);
            phb = wmean(fd_2018.phi(ib),1./fd_2018.phi_err(ib),1);
            ea = wmean(fd_2018.xpol(ia),1./fd_2018.phi_err(ia),1);
            eb = wmean(fd_2018.xpol(ib),1./fd_2018.phi_err(ib),1);


            if p.mce(inda(chind))~=0

                [phi_pair_2018(inda(chind)), poleff_pair_2018(inda(chind))] = calc_pair_diff_pol(pha,phb,ea,eb);

            else
                [poleff_pair_2018(inda(chind)), poleff_pair_2018(inda(chind))] = calc_pair_diff_pol(phb,pha,eb,ea);
            end

            if 0%~inrange(phi_pair_2018(inda(chind))+2.5,-2.5,2.5)
                phi_pair_2018(inda(chind)) = NaN;
                %poleff_pair_2018(inda(chind)) = NaN;
            end
        end
    end
end

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


%% Fig. 2.1.3 Averaged tile plots

clc
for valind = 1:2
    for pltind = 1:2
        vals = {nanmean(phi_pair,1) nanmean(poleff_pair,1)};
        valname = {'phi','xpol'};
        lims = {[-3.25 -1], [-1 1]*0.25; [0 1]*2e-4, [-1 1]*0.0001};
        labs = {'\phi_{pair} [Deg]','1-\epsilon_{pair}'};
        ttls = {'\phi_{pair} Tile Plot','\phi_{pair} Tile Plot - Median Subtracted';...
            '1-\epsilon_{pair} Tile Plot','1-\epsilon_{pair} Tile Plot - Median Subtracted'};
        pltname = {'','_medsub'};


        fig = figure(22+pltind);
        fig.Position(3:4) = [900 800];
        if pltind == 2
            for tileind = 1:20
                ind = p.tile==tileind;
                vals{valind}(ind) = vals{valind}(ind)-nanmedian(vals{valind}(ind));
            end
        end
        vals{valind}(ind90) = 0;
        plot_tiles(vals{valind},p,'fig',fig,'pair','sum','clim',lims{valind,pltind},'clab',labs{valind},'title',ttls{valind,pltind});
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

fd.phi_pair_medsub_notilt = fd.phi_pair_notilt;
fd.phi_pair_medsub = fd.phi_pair;
for chind = 1:length(p.gcp)
    idx = find(fd.ch==chind);
    if ~isempty(idx)
        fd.phi_pair_medsub(idx) = fd.phi_pair_medsub(idx)-nanmedian(fd.phi_pair_medsub(idx));
        fd.phi_pair_medsub_notilt(idx) = fd.phi_pair_medsub_notilt(idx)-nanmedian(fd.phi_pair_medsub_notilt(idx));
    end
end


idx = fd.phi_pair_medsub_notilt~=0 & fd.phi_pair_medsub~=0 & ~isnan(fd.phi_pair_medsub_notilt) & ~isnan(fd.phi_pair_medsub);
%idx = true(size(fd.ch));
edges = (-1:0.075:1)*0.6;
nexttile(1)
V = fd.phi_pair_medsub_notilt(idx);
N = histc(V,edges);
b = bar(edges,N,'histc');
b.FaceColor = cmlines(1,:);
M = nanmean(V);
S = nanstd(V);
L = length(find(~isnan(V)));
title({'No Tilt Meter Correction',...
    sprintf('M: %0.3f, STD: %0.3f, N: %0i',M,S,L)})
pbaspect([1 1 1])
xlim([-1 1]*0.6)
grid on
ylim([-0.1 2050])
xlabel('\phi_{pair} - median(\phi_{pair}) [Deg]')
ylabel('N')

nexttile(2)
V = fd.phi_pair_medsub(idx);
N = histc(V,edges);
b = bar(edges,N,'histc');
b.FaceColor = cmlines(1,:);
M = nanmean(V);
S = nanstd(V);
L = length(find(~isnan(V)));
title({'With Tilt Meter Correction',...
    sprintf('M: %0.3f, STD: %0.3f, N: %0i',M,S,L)})
pbaspect([1 1 1])
xlim([-1 1]*0.6)
grid on
ylim([-0.1 2050])
ax = gca;
ax.YTickLabel = {};
xlabel('\phi_{pair} - median(\phi_{pair}) [Deg]')

fname = 'tilt_corr_hist';
saveas(tiltcorrfig,fullfile(figdir,fname),'png')

%% Find the maximum STD for 2022 data subsets

Niter = 1000;

[M, S, N] = deal(NaN(Niter,1));
idx = NaN(Niter,5);
for iterind = 1:Niter
    ind1 = randperm(10,5);
    ind2 = find(~ismember(1:10,ind1));
    
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

for valind = 1%1:size(V,1)

for pltind = 1:size(V,2)
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
    saveas(fig,fullfile(figdir,fname))

end
end

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
%load('z:/dev/rps/rps_tilt_data_2022.mat')
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

%% Find the phis without tilt correction
[fd.phi_notilt, fd.xpol_notilt] = deal(NaN(size(fd.ch)));
for chind = 1:length(fd.ch)
    ind = find(fd_notilt.schnum==fd.schnum(chind) & fd_notilt.rowind==fd.rowind(chind) &...
        fd_notilt.ch==fd.ch(chind));
    if ~isempty(ind)
        fd.phi_notilt(chind) = nanmean(fd_notilt.phi(ind));
        fd.xpol_notilt(chind) = nanmean(fd_notilt.xpol(ind));
    end
end


% Create a bunch of extra stuff
% Also, create the phi and phi-diff per-obs 
clc
[fd.phi_pair, fd.poleff,fd.r,fd.theta,...
    fd.phi_pair_notilt, fd.poleff_notilt,...
    fd.xm,fd.ym,fd.thetam,fd.pair_idx] = deal(NaN(size(fd.ch)));
for chind = 1:length(fd.ch)
    
    % Phi/poleff pair
    isind0 = find(ismember(ind0,fd.ch(chind)));
    if ~isempty(isind0)
        has90 = find(fd.schnum == fd.schnum(chind) & fd.rowind==fd.rowind(chind) & fd.ch == ind90(isind0));
        if ~isempty(has90)
            fd.pair_idx(chind) = has90(1);
            [fd.phi_pair(chind) fd.poleff(chind)] = calc_pair_diff_pol(fd.phi(chind),...
                nanmean(fd.phi(has90)),fd.xpol(chind),nanmean(fd.xpol(has90)));
            [fd.phi_pair_notilt(chind) fd.poleff_notilt(chind)] = calc_pair_diff_pol(fd.phi_notilt(chind),...
                nanmean(fd.phi_notilt(has90)),fd.xpol_notilt(chind),nanmean(fd.xpol_notilt(has90)));
        end
    end

end

%%

[phis,phi_pair] = deal(NaN(length(dks),length(p.gcp)));
for obsind = 1:size(phi_pair,1)
    for chind = 1:size(phi_pair,2)
        ind = find(fd.obsnum==obsind & fd.ch==chind);
        if ~isempty(ind)
            if 1
            phis(obsind,chind) = nanmean(fd.phi(ind));
            phi_pair(obsind,chind) = nanmean(fd.phi_pair(ind));

            else
            phis(obsind,chind) = nanmean(fd.phi(ind)-fd.phi_s(ind)+90);
            phi_pair(obsind,chind) = nanmean(fd.phi_pair(ind)-fd.phi_s(ind)+90);
            end
        end
    end
end

%%

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


%% Do everything again, but for the type 9 scheds

for chind = 1:length(fd_type9.ch)
    % Tilt info
    if 1
        tiltcal = [0.284157 -0.0118];
        ind = lj_data.time >= fd_type9.t(chind) & lj_data.time<=fd_type9.t2(chind);
        fd_type9.tilt_out(chind) = nanmean(polyval(tiltcal,lj_data.AIN0(ind)));
        fd_type9.tilt_out_std(chind) = nanstd(polyval(tiltcal,lj_data.AIN0(ind)));
        fd_type9.tilt_temp(chind) = nanmean(lj_data.AIN2(ind)*100);
    end
end

% Create a bunch of extra stuff
% Also, create the phi and phi-diff per-obs 
clc
[fd_type9.obsnum,fd_type9.phi_pair, fd_type9.poleff,fd_type9.r,fd_type9.theta,...
    fd_type9.phi_pair_notilt, fd_type9.poleff_notilt,...
    fd_type9.xm,fd_type9.ym,fd_type9.thetam,fd_type9.pair_idx] = deal(NaN(size(fd_type9.ch)));
for chind = 1:length(fd_type9.ch)
    
    % Phi/poleff pair
    isind0 = find(ismember(ind0,fd_type9.ch(chind)));
    if ~isempty(isind0)
        has90 = find(fd_type9.schnum == fd_type9.schnum(chind) & fd_type9.rowind==fd_type9.rowind(chind) & fd_type9.ch == ind90(isind0));
        if ~isempty(has90)
            fd_type9.pair_idx(chind) = has90(1);
            [fd_type9.phi_pair(chind) fd_type9.poleff(chind)] = calc_pair_diff_pol(fd_type9.phi(chind),...
                nanmean(fd_type9.phi(has90)),fd_type9.xpol(chind),nanmean(fd_type9.xpol(has90)));
        end
    end

end

%

% Mirror coords
[fd_type9.theta, fd_type9.r] = cart2pol(fd_type9.x,fd_type9.y);
fd_type9.theta = wrapTo360(fd_type9.theta*180/pi);
fd_type9.thetam = wrapTo360(fd_type9.theta-fd_type9.dk_cen);
[fd_type9.xm, fd_type9.ym] = pol2cart(fd_type9.thetam*pi/180,fd_type9.r);

% Grab the moon-sun angles
% MJD, , ,raapp,decapp
fname = fullfile('z:/dev/sun_check_2022Aug12.txt');
f = fopen(fname);
hd_sun = textscan(f,'%f%s%s%f%f','delimiter',',','HeaderLines',60);
hd_sun{1} = hd_sun{1}-2400000.5;
fclose(f);
%fd_type9.az_cen_sun = interp1(hd_sun{1},hd_sun{4},fd_type9.t_cen);
fd_type9.az_cen_sun = wrapTo180(interp1(hd_sun{1},unwrap(hd_sun{4}*pi/180)*180/pi,fd_type9.t));
fd_type9.el_cen_sun = interp1(hd_sun{1},hd_sun{5},fd_type9.t);

%
clc
% MJD, , ,raapp,decapp
fname = fullfile('z:/dev/moon_check_2022Aug12.txt');
f = fopen(fname);
hd_moon = textscan(f,'%f%s%s%f%f','delimiter',',','HeaderLines',61);
hd_moon{1} = hd_moon{1}-2400000.5;
fclose(f);
%fd_type9.az_cen_moon = interp1(hd_moon{1},hd_moon{4},fd_type9.t_cen);
fd_type9.az_cen_moon = interp1(hd_moon{1},unwrap(hd_moon{4}*pi/180)*180/pi,fd_type9.t);
fd_type9.el_cen_moon = interp1(hd_moon{1},hd_moon{5},fd_type9.t);


%% phi vs. other stuff
% Need to run the above before we run this.

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
        'n1','n1','n1',...
        'n2','n2','n2'};

    yidx = {ismember(fd.ch,ind0), ismember(fd.ch,ind90), ismember(fd.ch,ind0),...
        ismember(fd.ch,ind0), ismember(fd.ch,ind90), ismember(fd.ch,ind0),...
        ismember(fd.ch,ind0), ismember(fd.ch,ind90), true(size(fd.ch)),...
        ismember(fd.ch,ind0), ismember(fd.ch,ind90), true(size(fd.ch))};
    
    yidx_type9 = {ismember(fd_type9.ch,ind0), ismember(fd_type9.ch,ind90), ismember(fd_type9.ch,ind0),...
        ismember(fd_type9.ch,ind0), ismember(fd_type9.ch,ind90), ismember(fd_type9.ch,ind0),...
        ismember(fd_type9.ch,ind0), ismember(fd_type9.ch,ind90), true(size(fd_type9.ch)),...
        ismember(fd_type9.ch,ind0), ismember(fd_type9.ch,ind90), true(size(fd_type9.ch))};

    ynames = {'phi_0','phi_90','phi_p',...
        'xpol_0','xpol_90','xpol_p',...
        'n1_0','n1_90','n1_p',...
        'n2_0','n2_90','n2_p',...
        };

    Ylabs = {'\phi_0 [Degrees]','\phi_{90} [Degrees]','\phi_{pair} [Degrees]',...
        '\epsilon_0','\epsilon_{90}','1-Pol Eff.',...
        'n1_0','n1_{90}','n1',...
        'n2_0','n2_{90}','n2'};
    
    %yttls = {'Pol 0','Pol 90','Pair-Diff'};
    if medind == 1
        ylims = {[-4.5 0], [-4.5 0]+90, [-4.5 0],...
            [-1 1]*0.02,[-1 1]*0.02, [-1 10]*1e-4,...
            [-1 1]*0.07,[-1 1]*0.07,[-1 1]*0.07,...
            [-1 1]*0.07,[-1 1]*0.07,[-1 1]*0.07,...
            };
    else
        ylims = {[-1 1]*1, [-1 1]*1, [-1 1]*0.5,...
            [-1 1]*0.01,[-1 1]*0.01, [-1 1]*6e-4,...
            [-1 1]*0.07,[-1 1]*0.07,[-1 1]*0.07...
            [-1 1]*0.07,[-1 1]*0.07,[-1 1]*0.07...
            };
    end

    X = {'t','obsnum','az_cen','el_cen','dk_cen','az_cen_sun','el_cen_sun',...
        'x','y','xm','ym','r','theta','thetam','tod',...
        'tilt_out','tilt_out_std','tilt_temp'};
    Xlabs = {'Time [dd-mmm]','Observation Number','Azimuth [Degrees]','Elevation [Degrees]',...
        'Deck [Degrees]','Sun Azimuth [Degrees]','Sun Elevation [Degrees]',...
        'Inst.-Fixed X [Degrees]','Inst.-Fixed Y [Degrees]',...
        'Mirror-Fixed X [Degrees]','Mirror-Fixed Y [Degrees]','Inst.-Fixed r [Degrees]',...
        'Inst.-Fixed \theta [Degrees]','Mirror-Fixed \theta [Degrees]','Time-of-Day [Days]',...
        'Tilt Out [Degrees]','Tilt Out STD [Degrees]','Tilt Temp [deg-C]'};
    

    for yind = 1:size(Y,2)

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


        for xind = [2 10:14]%1:length(X)
            clf; hold on;
            if xind == 1
                x = mjd2datenum(fd.(X{xind})(yidx{yind}));
            else
                x = fd.(X{xind})(yidx{yind});
            end
            y = fd.(Y{yind})(yidx{yind})-medvals(fd.ch(yidx{yind}));
            scatter(x,y,14,p.tile(fd.ch(yidx{yind}))','filled')

            % Type 9 plot
            if xind == 1
                x = mjd2datenum(fd_type9.(X{xind})(yidx_type9{yind}));
            else
                x = fd_type9.(X{xind})(yidx_type9{yind});
            end
            y = fd_type9.(Y{yind})(yidx_type9{yind})-medvals(fd_type9.ch(yidx_type9{yind}));
            plot(x,y,'x','Color',cmlines(2,:))
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
    z(obsind) = plot(V1(obsind,idx)',V2(obsind,idx)','.','MarkerSize',14);
end
xlim(lims)
ylim(lims)
grid on
xlabel({'\phi_{d,0}-Md(\phi_{d,0}) [Degrees]',''})
ylabel({'','\phi_{d,90}-Md(\phi_{d,90}) [Degrees]'})
leg = legend(z,legttls);
title(leg,'DK''s:')
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


%% Plot phis vs n1 & n2

[medvals] = deal(NaN(1,2640));
for chind = 1:length(medvals)
    ind = fd.ch==chind;
    medvals(chind) = median(fd.phi(ind));
end

fig = figure(1);
fig.Position(3:4) = [600 500];
clf; hold on;

medsub = fd.phi-medvals(fd.ch);
idx = ismember(fd.ch,ind90) & medsub~=0;
plot(fd.n1(idx),medsub(idx),'.')
plot([-1 1]*0.05,[-1 1],'k--')
grid on
xlabel('n2')
ylabel('\phi_{90}-Md(\phi_{90}) [Degrees]')
xlim([-1 1]*0.05)
ylim([-1 1]*1)
% leg = legend(z,legttls);
% title(leg,'DK''s:')
% title(corrttls{corrind})

%% Plot modecurve residuals!?

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

%%
idx = ismember(fd.ch,ind90);
cov([fd.phi_medsub(idx)',fd.xpol(idx)',fd.n1(idx)',fd.n2(idx)'])


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

%% See if there are differences in phi_s a/b
clf
ind = ~isnan(fd.pair_idx);
plot(fd.phi_medsub(find(ind)),fd.phi_medsub(fd.pair_idx(ind)),'.')



%% 

clf
ind = p.tile(ind0)==3;
plot(repmat(dks',1,length(ind0(ind)))',phis(:,ind0(ind))'-phis(:,ind90(ind))'+90,'.')












