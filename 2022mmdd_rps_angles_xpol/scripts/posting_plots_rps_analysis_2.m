function posting_plots_rps_analysis_2()

%%
addpath('z:/dev/rps/')
addpath('z:/pipeline/beammap')
addpath('z:/pipeline/util')
addpath('z:/dev/diff_polarization/')
addpath('z:/dev')
%%
clear all
close all
clc

% Load no tilt data
load('z:/dev/rps/rps_beam_fits_type5_notilt_cut.mat')
fd_notilt = fd;


% Load type 9 data
load('z:/dev/rps/type9_fit_dat_mirrorfitted_cut.mat')
load('z:/dev/rps/sch_type9.mat')
fd_type9 = fd;

% Load 2018 data
load('z:/dev/rps/rps_beam_fits_cut_2018.mat')
fd.xpol = fd.aparam(:,2);
fd.phi = fd.phi_d;
fd.phi_err = fd.aerr(:,1);
fd.schnum = fd.sch;
fd_2018 = fd;

load('z:/dev/rps/rps_beam_fits_type5_final_cut')
load('z:/dev/rps/fpu_data_obs.mat')
load('z:/pipeline/beammap/viridis_cm.mat')
load('z:/dev/sims/coaddopt_6600.mat')
load('z:/dev/rps/sch_type5.mat')
p_ind = coaddopt.ind;
figdir = fullfile('C:','Users','James','Documents','GitHub','postings','2022mmdd_rps_angles_xpol','figs','');

cmlines = colormap('lines');

[fd.theta, fd.r] = pol2cart(fd.x,fd.y);
fd.theta = fd.theta*180/pi;

fd.phi_corr = mirror_diff_pol_calc(fd.phi,fd.r,fd.theta,2200,2200,fd.dk_cen,fd.mirr_tilt,fd.mirr_roll);
fd.phi_corr = reshape(fd.phi_corr,size(fd.ch));

for chind = 1:length(fd.ch)
    % Obs Number
    for schind = 1:length(scheds)
        if ismember(fd.schnum(chind),scheds{schind})
            fd.obsnum(chind) = schind;
        end
    end
end

% If 1 use chflags from B18
if 1
inda = p_ind.rgl100a;
indb = p_ind.rgl100b;

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

len = length(scheds);

%% Grab the phis per obs
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

for pltind = 1:2%1:size(V,2)
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

%load('z:/dev/rps/rps_tilt_data_2022.mat')

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


%%
clc
[fd.t, fd.obsnum,fd.phi_pair, fd.poleff,fd.r,fd.theta,...
    fd.xm,fd.ym,fd.thetam,fd.phi_pair_corr,fd.phi_pair_notilt,...
    fd.poleff_corr,fd.poleff_notilt] = deal(NaN(size(fd.ch)));
for chind = 1:length(fd.ch)
    s = sch{fd.schnum(chind)};
    idx = s.index(fd.rowind(chind),1);
    fd.t(chind) = s.scans(idx).t1;
    
    % Obs Number
    for schind = 1:length(scheds)
        if ismember(fd.schnum(chind),scheds{schind})
            fd.obsnum(chind) = schind;
        end
    end

    % Phi/poleff pair
    isind0 = find(ismember(ind0,fd.ch(chind)));
    if ~isempty(isind0)
        has90 = find(fd.schnum == fd.schnum(chind) & fd.rowind==fd.rowind(chind) & fd.ch == ind90(isind0));
        if ~isempty(has90)
            [fd.phi_pair(chind) fd.poleff(chind)] = calc_pair_diff_pol(fd.phi(chind),...
                nanmean(fd.phi(has90)),fd.xpol(chind),nanmean(fd.xpol(has90)));
            [fd.phi_pair_corr(chind) fd.poleff_corr(chind)] = calc_pair_diff_pol(fd.phi_corr(chind),...
                nanmean(fd.phi_corr(has90)),fd.xpol(chind),nanmean(fd.xpol(has90)));
            [fd.phi_pair_notilt(chind) fd.poleff_notilt(chind)] = calc_pair_diff_pol(fd.phi_notilt(chind),...
                nanmean(fd.phi_notilt(has90)),fd.xpol_notilt(chind),nanmean(fd.xpol_notilt(has90)));
        end
    end

    % Tilt info
    if 0
        idx2 = max(s.index(fd.rowind(chind),:));
        fd.t2(chind) = s.scans(idx2).t2;
        ind = lj_data.time >= fd.t(chind) & lj_data.time<=fd.t2(chind);
        fd.tilt_out(chind) = nanmean(lj_data.AIN0(ind));
        fd.tilt_out_std(chind) = nanstd(lj_data.AIN0(ind));
        fd.tilt_temp(chind) = nanmean(lj_data.AIN2(ind));
    end
end

[phis,phi_pair,...
    phis_corr,phi_pair_corr,...
    phis_notilt,phi_pair_notilt] = deal(NaN(length(dks),length(p.gcp)));
for obsind = 1:size(phi_pair_corr,1)
    for chind = 1:size(phi_pair_corr,2)
        ind = find(fd.obsnum==obsind & fd.ch==chind);
        if ~isempty(ind)
            phis(obsind,chind) = nanmean(fd.phi(ind));
            phi_pair(obsind,chind) = nanmean(fd.phi_pair(ind));
            phis_corr(obsind,chind) = nanmean(fd.phi_corr(ind));
            phi_pair_corr(obsind,chind) = nanmean(fd.phi_pair_corr(ind));
            phis_notilt(obsind,chind) = nanmean(fd.phi_notilt(ind));
            phi_pair_notilt(obsind,chind) = nanmean(fd.phi_pair_notilt(ind));
        end
    end
end



% Mirror coords
[fd.theta, fd.r] = cart2pol(fd.x,fd.y);
fd.theta = wrapTo360(fd.theta*180/pi);
fd.thetam = wrapTo360(fd.theta-fd.dk_cen);
[fd.xm, fd.ym] = pol2cart(fd.thetam*pi/180,fd.r);

%% Grab the moon-sun angles
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


%% phi vs. other stuff

fd.az_cen = wrapTo360(fd.az_cen);
fd.tod = fd.t-floor(fd.t);
fd.xpol_corr = fd.xpol;

clc
medsubnames = {'','_medsub'};
for medind = 1:length(medsubnames)
    corrnames = {'','_corr','_notilt'};
    for corrind = 2%1:length(corrnames)
        Y = {'phi', 'phi','phi_pair';
            'xpol','xpol','poleff'};

        yidx = {ismember(fd.ch,ind0), ismember(fd.ch,ind90), ismember(fd.ch,ind0)};
        ynames = {'phi_a','phi_b','phi_p';...
            'xpol_a','xpol_b','xpol_p'};
        yttls = {'Pol A','Pol B','Pair-Diff'};
        if medind == 1
        ylims = {[-4.5 0], [-4.5 0]+90, [-4.5 0];...
            [-1 1]*0.02,[-1 1]*0.02, [-1 10]*1e-4};
        else
        ylims = {[-1 1]*1, [-1 1]*1, [-1 1]*0.5;...
            [-1 1]*0.01,[-1 1]*0.01, [-1 1]*6e-4};
        end

        X = {'t','obsnum','az_cen','el_cen','dk_cen','az_cen_sun','el_cen_sun',...
            'x','y','xm','ym','r','theta','thetam','tod'};
        Xlabs = {'Time [MJD]','Observation Number','Azimuth [Degrees]','Elevation [Degrees]',...
            'Deck [Degrees]','Inst. Fixed X [Degrees]','Inst. Fixed Y [Degrees]',...
            'Mirror-Fixed X [Degrees]','Mirror-Fixed Y [Degrees]','Inst. Fixed r [Degrees]',...
            'Inst. Fixed \theta [Degrees]','Time-of-Day [Days]'};
        fig = figure(51);
        fig.Position(3:4) = [900 400];
        clf


        for valind = 1:size(Y,1)
            for yind = 1:size(Y,2)
                Y{valind,yind} = [Y{valind,yind} corrnames{corrind}];
                
                if medind == 2
                [medvals] = deal(NaN(1,2640));
                for chind = 1:length(medvals)
                    ind = fd.ch==chind;
                    medvals(chind) = nanmedian(fd.(Y{valind,yind})(ind));
                end
                else
                    medvals = zeros(1,2640);
                end

                for xind = 1:length(X)
                    if 0
                        scatter(fd.(X{xind})(yidx{yind}),fd.(Y{valind,yind})(yidx{yind}),14,fd.obsnum(yidx{yind}),'filled')
                        ylim(ylims{valind,yind})
                    else
                        y = fd.(Y{valind,yind})(yidx{yind})-medvals(fd.ch(yidx{yind}));
                        scatter(fd.(X{xind})(yidx{yind}),y,14,fd.dk_cen(yidx{yind}))%,'filled')
                        ylim(ylims{valind,yind})
                    end
                    colormap(cm)
                    grid on
                    xlabel(Xlabs{xind})
                    ylabel(Y{valind,yind})
                    fname = sprintf('scatter_%s_vs_%s%s%s.png',ynames{valind,yind},X{xind},corrnames{corrind},medsubnames{medind});
                    saveas(fig,fullfile(figdir,fname))
                end
            end
        end
    end
end

%% Compare phis to diff-pol corrected

clc
% grab medians for subtraction
[medvals, medvals_corr] = deal(NaN(2640,1));
[fd.phi_medsub, fd.phi_corr_medsub] = deal(NaN(size(fd.ch)));

for chind = 1:length(medvals)
    ind = fd.ch==chind;
    medvals(chind) = nanmedian(fd.phi(ind));
    medvals_corr(chind) = nanmedian(fd.phi_corr(ind));
    fd.phi_medsub(ind) = fd.phi(ind)-medvals(chind);
    fd.phi_corr_medsub(ind) = fd.phi_corr(ind)-medvals_corr(chind);
end

%

fig = figure(1);
clf;
V = {fd.phi_medsub,fd.phi_corr_medsub};

for pltind = 1:2
subplot(1,2,pltind)
edges = (-1:0.075:1)*0.6;
N = histc(V{pltind},edges);
b = bar(edges,N,'histc');
b.FaceColor = cmlines(1,:);
grid on
xlabel('\phi-<\phi> [Degrees]')
title({'No Pol Correction',...
    sprintf('M: %0.4f | S: %0.4f | N: %i | EOM: %0.4f',...
    nanmean(V{pltind}),nanstd(V{pltind}),length(find(~isnan(V{pltind}))),...
    nanstd(V{pltind})./sqrt(length(find(~isnan(V{pltind})))))})
xlim([edges([1 end])])
end

%% Look at hists of phis vs correction
% This does a good job of showing that 

vals = {phis_notilt,phis,phis_corr};
valttls = {'No Correction','Tilt Only','Tilt + Diff-Pol'};

fig = figure(1);
fig.Position(3:4) = [1100,400];
clf;

for pltind = 1:length(vals)
    subplot(1,3,pltind)
    V = nanmean(vals{pltind}(:,ind0),2);
    edges = (-1:0.15:1)*0.15-2.36;
    N = histc(V,edges);
    b = bar(edges,N,'histc');
    b.FaceColor = cmlines(1,:);
    M = nanmean(V);
    S = nanstd(V);
    L = length(find(~isnan(V)));
    title({valttls{pltind},...
        sprintf('M: %0.4f | S: %0.4f | N: %i | EOM: %0.4f',M,S,L,S./sqrt(L))})
    grid on
end
sgtitle('Pol 0')

% Pol B
vals = {phis_notilt,phis,phis_corr};
valttls = {'No Correction','Tilt Only','Tilt + Diff-Pol'};

fig = figure(2);
fig.Position(3:4) = [1100,400];
clf;

for pltind = 1:length(vals)
    subplot(1,3,pltind)
    V = nanmean(vals{pltind}(:,ind90),2);
    edges = (-1:0.15:1)*0.15+87.4;
    N = histc(V,edges);
    b = bar(edges,N,'histc');
    b.FaceColor = cmlines(1,:);
    M = nanmean(V);
    S = nanstd(V);
    L = length(find(~isnan(V)));
    title({valttls{pltind},...
        sprintf('M: %0.4f | S: %0.4f | N: %i | EOM: %0.4f',M,S,L,S./sqrt(L))})
    grid on
end
sgtitle('Pol 90')


% Pair
vals = {phi_pair_notilt,phi_pair,phi_pair_corr};
valttls = {'No Correction','Tilt Only','Tilt + Diff-Pol'};

fig = figure(3);
fig.Position(3:4) = [1100,400];
clf;

for pltind = 1:length(vals)
    subplot(1,3,pltind)
    V = nanmean(vals{pltind}(:,ind0),2);
    edges = (-1:0.15:1)*0.15-2.5;
    N = histc(V,edges);
    b = bar(edges,N,'histc');
    b.FaceColor = cmlines(1,:);
    M = nanmean(V);
    S = nanstd(V);
    L = length(find(~isnan(V)));
    title({valttls{pltind},...
        sprintf('M: %0.4f | S: %0.4f | N: %i | EOM: %0.4f',M,S,L,S./sqrt(L))})
    grid on
end
sgtitle('Pol Pair')

%% Plot Mean phi vs DK






