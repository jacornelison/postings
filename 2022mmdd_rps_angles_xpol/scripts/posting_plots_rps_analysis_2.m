function posting_plots_rps_analysis_2()

%%
clear all
close all
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
p_ind = coaddopt.ind;
figdir = fullfile('C:','Users','James','Documents','GitHub','postings','2022mmdd_rps_angles_xpol','figs','');
%% Grab the chis

inda = p_ind.a;
indb = p_ind.b;

ind0 = [p_ind.rgl100a(ismember(p_ind.rgl100a,find(p.mce~=0))) ...
    p_ind.rgl100b(ismember(p_ind.rgl100b,find(p.mce==0)))];
ind90 = [p_ind.rgl100b(ismember(p_ind.rgl100b,find(p.mce~=0))) ...
    p_ind.rgl100a(ismember(p_ind.rgl100a,find(p.mce==0)))];

len = length(scheds);
clc

[phis, xpols, phis_err, xs, ys] = deal(NaN(len,2640));
for schind = 1:len
    for chind = 1:2640
        ci = find(fd.ch==chind & ismember(fd.schnum,scheds{schind}));

        if ~isempty(ci)

            %ch_chi = atand(tand(p.theta(chind)+p.chi(chind)));

            val = reshape(fd.phi(ci),[],1);
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
[phi_pair, poleff_pair, phi_pair_err] = deal(NaN(len,2640));
for schind = 1:len
    for chind = 1:length(inda)

        pha = phis(schind,inda(chind));
        phb = phis(schind,indb(chind));
        ea = xpols(schind,inda(chind));
        eb = xpols(schind,indb(chind));

        phi_pair_err(schind,inda(chind)) = max(phis_err(schind,[inda(chind) indb(chind)]));
        if p.mce(inda(chind))~=0

            [phi_pair(schind,inda(chind)), poleff_pair(schind,inda(chind))] = calc_pair_diff_pol(pha,phb,ea,eb);
        else

            [phi_pair(schind,inda(chind)), poleff_pair(schind,inda(chind))] = calc_pair_diff_pol(phb,pha,eb,ea);
        end

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
        for chind = 1:length(inda)

            pha = phis_2018(schind,inda(chind));
            phb = phis_2018(schind,indb(chind));
            ea = xpols_2018(schind,inda(chind));
            eb = xpols_2018(schind,indb(chind));

            phi_pair_err_2018(schind,inda(chind)) = max(phis_err_2018(schind,[inda(chind) indb(chind)]));
            if p.mce(inda(chind))~=0

                [phi_pair_2018(schind,inda(chind)), poleff_pair_2018(schind,inda(chind))] = calc_pair_diff_pol(pha,phb,ea,eb);
            else

                [phi_pair_2018(schind,inda(chind)), poleff_pair_2018(schind,inda(chind))] = calc_pair_diff_pol(phb,pha,eb,ea);
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


%% Fig 2.1 Angle vs. channel

phi_fiducial = atand(tand(p.chi_thetaref+p.chi));
yrnames = {'2022','2018'};
V = {phis, phis_2018};
V2 = {phi_pair, phi_pair_2018};
for yearind = 1:2


fig = figure(21);
fig.Position(3:4) = [900 500];
clf;

clc
cmlines = colormap('lines');
for dkind = 1:length(dks)
    X = {ind0, ind90, 1:2640};
    Y = {V{yearind}(dkind,ind0), V{yearind}(dkind,ind90), V2{yearind}(dkind,:)};
    ylab = {'\phi_0','\phi_{90}','\phi_{pair}'};
    ylims = {[-5 1], [85 91], [-5 1]};

    for pltind = 1:length(X)


        subplot(3,1,pltind)
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

yrnames = {'2022','2018'};
V = {xpols, xpols_2018};
V2 = {poleff_pair,poleff_pair_2018};
for yearind = 1:2


fig = figure(21);
fig.Position(3:4) = [900 500];
clf;


cmlines = colormap('lines');
for dkind = 1:length(dks)
    X = {ind0, ind90, 1:2640};
    Y = {V{yearind}(dkind,ind0), V{yearind}(dkind,ind90), V2{yearind}(dkind,:)};
    ylab = {'\epsilon_0','\epsilon_{90}','POLEFF_{pair}'};
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








