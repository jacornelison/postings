function rps_pointing_plots_Aug2022()


%load('z:/dev/rps/rps_beam_fits_type11_rerun_cut.mat')
load('z:/dev/rps/rps_beam_fits_type5_rerun_cut.mat')
%load('z:/dev/rps/rps_beam_fits_rerun_all_cut.mat')
figdir = 'c:/Users/James/Documents/GitHub/postings/2022mmdd_rps_pointing/figs/';
load('z:/dev/rps/fpu_data_obs.mat')
load('z:/dev/rps/pm.mat')
load('z:/dev/rps/source_fit_data.mat')
inrange = @(A,B,C) B <= A & A <= C;
outrange = @(A,B,C) A <= B | C <= A;
addpath('z:/dev')
addpath('z:/pipeline/beammap')
addpath('z:/pipeline/util')
addpath('z:/dev/rps')


prx = 2*sind(p.r/2).*cosd(p.theta)*180/pi;
pry = 2*sind(p.r/2).*sind(p.theta)*180/pi;

clc
% Moon-derived mirror params -- Use only these and nothing else!
mirror = struct();
mirror.height = 1.4592;
mirror.tilt= 44.8870;
mirror.roll = -0.070;
rpsopt.mirror = mirror;

rpsopt.source.distance = 195.5;
% Fit for the source params given our mirror info:
source = rps_fit_source(fd,rpsopt,p,'');
rpsopt.source = source;


% With new mirror and source parameters, update the pointing.
[fd.x,fd.y,phi] = beam_map_pointing_model(fd.az_cen,fd.el_cen,fd.dk_cen,model,'bicep3',mirror,source,[]);
fd.x = reshape(fd.x,size(fd.ch));
fd.y = reshape(fd.y,size(fd.ch));
fd.resx = reshape(prx(fd.ch),size(fd.ch))-fd.x;
fd.resy = reshape(pry(fd.ch),size(fd.ch))-fd.y;
[resth, resr] = cart2pol(fd.resx,fd.resy);
[fd.resx_rot, fd.resy_rot] = pol2cart(resth-fd.dk_cen*pi/180,resr);

% Order the a/b phis based on schedule
[phis, xpols, phis_err, xs, ys] = deal(NaN(len,2640));
for schind = 1:len
    for chind = 1:2640
        ci = find(fd.ch==chind & ismember(fd.schnum,scheds{schind}));

        if ~isempty(ci)

            ch_chi = atand(tand(p.theta(chind)+p.chi(chind)));

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
[phi_pair, xpol_pair, phi_pair_err] = deal(NaN(len,2640));
for schind = 1:len
    for chind = 1:length(inda)

        pha = phis(schind,inda(chind));
        phb = phis(schind,indb(chind));
        ea = xpols(schind,inda(chind));
        eb = xpols(schind,indb(chind));

        phi_pair_err(schind,inda(chind)) = max(phis_err(schind,[inda(chind) indb(chind)]));
        if p.mce(inda(chind))~=0

            [phi_pair(schind,inda(chind)), xpol_pair(schind,inda(chind))] = calc_pair_diff_pol(pha,phb,ea,eb);
        else

            [phi_pair(schind,inda(chind)), xpol_pair(schind,inda(chind))] = calc_pair_diff_pol(phb,pha,eb,ea);
        end

        if 0%~inrange(phi_pair(schind,inda(chind))+2.5,-2.5,2.5)
            phi_pair(schind,inda(chind)) = NaN;
            xpol_pair(schind,inda(chind)) = NaN;
            phi_pair_err(schind,inda(chind)) = 1e10;
        end
    end
    %ind = abs(phi_pair(schind,:)-nanmean(phi_pair(schind,:)))<1;
end



%%
% Organize the x/y pointing by obs
[xps, yps] = deal(NaN(len,2640));
for schind = 1:len
    for chind = 1:2640
        ci = find(fd.ch==chind & ismember(fd.schnum,scheds{schind}));

        if ~isempty(ci)

            %err(isnan(val))=1e10;
            xps(schind,chind) = mean(fd.x(ci));
            yps(schind,chind) = mean(fd.y(ci));

        end
    end
end

% Fit for a rotation and scaling

fpu = fit_fpu_angle_and_scaling_from_xy(fd,p);

%% Look at mean quiver plots with/without rot & scaling applied.

winscale = 1.5;
scaling = 10;
fig = figure(1);
fig.Position(3:4) = [500*winscale 450*winscale];
clf;
cm = colormap('turbo');
clridx = floor(linspace(1,size(cm,1),3*19));

% Things dealing with projection
xlims = {[-1 1]*15 [-1 1]*15 [-1 1]*0.5};
ylims = {[-1 1]*15 [-1 1]*15 [-0.5 0.8]};
projlabels = {' [Degrees]','_m [Degrees]','_m [Meters]'};
projnames = {'','_mirror','_mirror'};

% Things dealing with fits
fittype = {'overall','perdk'};%,'persch'};


casename = {'none','ang','scale','both'};
parms = [fpu.angle,fpu.scaling,fpu.xtrans,fpu.ytrans];
for caseind = 1:4
    switch caseind
        case 1
            params = [0,1,0,0];
            x = xps;
            y = yps;
            caselabel = sprintf('FPU Rot: %0.3f^o, Scaling: %0.3f',params(1),params(2));
        case 2
            params = parms.*[1,0,0,0];
            params(2) = 1;
            x = params(2).*(xps.*cosd(params(1))-yps.*sind(params(1)));
            y = params(2).*(xps.*sind(params(1))+yps.*cosd(params(1)));
            caselabel = sprintf('FPU Rot: %0.3f^o, Scaling: %0.3f',params(1),params(2));
        case 3
            params = parms.*[0,1,0,0];
            x = params(2).*(xps.*cosd(params(1))-yps.*sind(params(1)));
            y = params(2).*(xps.*sind(params(1))+yps.*cosd(params(1)));
            caselabel = sprintf('FPU Rot: %0.3f^o, Scaling: %0.3f',params(1),params(2));
        case 4
            params = parms;
            x = params(2).*(xps.*cosd(params(1))-yps.*sind(params(1)));
            y = params(2).*(xps.*sind(params(1))+yps.*cosd(params(1)));
            caselabel = sprintf('FPU Rot: %0.3f^o, Scaling: %0.3f',params(1),params(2));
    end

    projind = 1;
    corrind = 2;
    winscale = 1;
    scaling = 10;
    fig = figure(1);
    fig.Position(3:4) = [500*winscale 450*winscale];
    clf;
    ind = 1:length(scheds);
    quiver(prx,pry,(prx-nanmean(x(ind,:),1)')*scaling,(pry-nanmean(y(ind,:),1)')*scaling,0)
    grid on;
    xlim(xlims{projind})
    ylim(ylims{projind})
    xlabel(sprintf('X%s',projlabels{projind}))
    ylabel(sprintf('Y%s',projlabels{projind}))
    title({sprintf('Beam Center Residuals, x%i', scaling),...
        sprintf('All-DK average, Tilt: %1.2f  Roll: %1.3f',mirror.tilt,mirror.roll),...
        caselabel...
        })
    figname = fullfile(figdir,sprintf('quiver_mean_fit_%s.png',casename{caseind}));
    saveas(fig,figname)

end

%% Look at fits vs DK.

clc
[fpuparms, nchans, fdval, times] = deal([]);
unqsch = unique(fd.schnum);
% for schedind = 1:length(unqsch)
%     ind = ismember(fd.schnum,unqsch(schedind));
for schedind = 1:length(scheds)
    ind = ismember(fd.schnum,scheds{schedind});
    if ~isempty(find(ind))
        fpu = fit_fpu_angle_and_scaling_from_xy(structcut(fd,ind),p);
        fpuparms(end+1,:) = [fpu.angle,fpu.scaling,fpu.xtrans,fpu.ytrans];
        nchans(end+1) = length(find(ind));
        fdval(end+1) = nanmean(fd.dk_cen(ind));
        times(end+1) = nanmean(fd.t(ind));

    end
end

fig = figure(1);
fig.Position(3:4) = [900 300];
clf; hold on;
%scatter(times-floor(times),fpuparms(:,1),14,fdval,'filled')
scatter(-fdval,fpuparms(:,1),14,'filled');%times-floor(times),'filled')
grid on
%xlim([-100 200])
%ylim([-0.12 -0.04])
%C = colorbar();
%C.Label.String = ''
xlabel('DK {Deg]')
ylabel('FPU Rotation Angle [Deg]')
title('Per-Obs estimated Focal Plane angle Vs. DK')
saveas(fig,'figs\fpu_rot_vs_dk.png')

%% Fix source loc and fit for mirror.
% How does it look vs time compared to
clc
[mirrparms, nchans, times] = deal([]);
for schedind = 1:length(scheds)
    ind = ismember(fd.schnum,scheds{schedind});

    if ~isempty(find(ind))
        mirror = rps_fit_mirror(structcut(fd,ind),rpsopt,p,'');
        mirrparms(schedind,:) = [mirror.tilt,mirror.roll];
        nchans(schedind) = length(find(ind));
        times(end+1) = nanmean(fd.t(ind));

    end
end

fig = figure(1);
fig.Position(3:4) = [900 300];
clf; hold on;
scatter(times-floor(times),mirrparms(:,2),14,times,'filled')
grid on
%xlim([-100 200])
ylim([-0.12 -0.04])


%% Fix source loc and fit for mirror.
% How does it look vs time compared to
clc
[mirrparms, nchans, times, fdval] = deal([]);
unqsch = unique(fd.schnum);
for schedind = 1:length(unqsch)
    for rowind = 1:19
        ind = fd.schnum==unqsch(schedind) & fd.rowind==rowind;

        if ~isempty(find(ind))
            mirror = rps_fit_mirror(structcut(fd,ind),rpsopt,p,'');
            mirrparms(end+1,:) = [mirror.tilt,mirror.roll];
            nchans(end+1) = length(find(ind));
            times(end+1) = nanmean(fd.t(ind));
            fdval(end+1) = nanmean(fd.dk_cen(ind));
            %times(end+1) = sch{unqsch(schedind)}.scans(sch{unqsch(schedind)}.index(rowind,1)).t1;
            %fdval(end+1) = sch{unqsch(schedind)}.el_ofs(rowind,1);

        end
    end
end

%
fig = figure(1);
fig.Position(3:4) = [900 300];
clf; hold on;
scatter(times-0*floor(times),mirrparms(:,2),14,fdval,'filled')
grid on
%xlim([-100 200])
%ylim([-0.12 -0.04])

%% Look at the diff between CMB-derived phi and our phi

clc
[phia0, phib0] = deal(NaN(size(p.theta)));
phi_fiducial = atand(tand(p.theta+p.chi));
ia = ismember(1:2640,p_ind.a)'&p.mce==0;
ib = ismember(1:2640,p_ind.b)'&p.mce==0;
phib0(ia) = phi_fiducial(ib);
phia0(ia) = phi_fiducial(ia);

ia = ismember(1:2640,p_ind.rgla)'&p.mce~=0;
ib = ismember(1:2640,p_ind.rglb)'&p.mce~=0;
phib0(ia) = phi_fiducial(ia);
phia0(ia) = phi_fiducial(ib);

[phi_pair0 pol_eff0] = calc_pair_diff_pol(phia0,phib0,0,0);

dphi = nanmean(phi_pair,1)-phi_pair0';
ind = abs(dphi)<1.5;
%ind = true(size(dphi));
nanmean(dphi(ind))

ch = 1:2640;
fig = figure(1);
clf; hold on;
scatter(ch(ind),dphi(ind),14,'filled')


%% Fix source loc and fit for mirror. Moon and RPS
% How does it look vs time compared to
clear all
clc
load('z:/dev/rps/fpu_data_obs.mat')
load('z:/dev/rps/pm.mat')

[mirrparms, nchans, times, fdval] = deal([]);
for targind = 1:2
    switch targind
        case 1
            load('z:/dev/rps/moon_beam_fits_phase_corrected_cut.mat')
            unqsch = unique(fd.schind);
            for schedind = 1:length(unqsch)
                ind0 = fd.schind==unqsch(schedind);
                ind = ind0;
                for scanind = 1:19
                    ind = ind0 & fd.scanind==scanind;

                    if ~isempty(find(ind))
                        fd0 = moon_fit_mirror(structcut(fd,ind),'p',p,'p_ind',p_ind,'savedir','','pm',model);
                        mirrparms(end+1,:) = fd0.fitparam;
                        nchans(end+1) = length(find(ind));
                        times(end+1) = nanmean(fd.t_cen(ind));
                        %fdval(end+1) = -nanmean(fd.dk_cen(ind));
                        fdval(end+1) = targind;
                    end
                end
            end

        case 2
            load('z:/dev/rps/rps_beam_fits_rerun_all_cut.mat')
            %load('z:/dev/rps/rps_beam_fits_type5_rerun_cut.mat')
            load('z:/dev/rps/source_fit_data.mat')
            % Moon-derived mirror params -- Use only these and nothing else!
            mirror = struct();
            mirror.height = 1.4592;
            mirror.tilt= 44.8870;
            mirror.roll = -0.070;
            rpsopt.mirror = mirror;

            rpsopt.source.distance = 195.5;
            % Fit for the source params given our mirror info:
            source = rps_fit_source(fd,rpsopt,p,'');
            rpsopt.source = source;
            prx = 2*sind(p.r/2).*cosd(p.theta)*180/pi;
            pry = 2*sind(p.r/2).*sind(p.theta)*180/pi;


            % With new mirror and source parameters, update the pointing.
            [fd.x,fd.y,phi] = beam_map_pointing_model(fd.az_cen,fd.el_cen,fd.dk_cen,model,'bicep3',mirror,source,[]);
            fd.x = reshape(fd.x,size(fd.ch));
            fd.y = reshape(fd.y,size(fd.ch));
            fd.resx = reshape(prx(fd.ch),size(fd.ch))-fd.x;
            fd.resy = reshape(pry(fd.ch),size(fd.ch))-fd.y;
            [resth, resr] = cart2pol(fd.resx,fd.resy);
            [fd.resx_rot, fd.resy_rot] = pol2cart(resth-fd.dk_cen*pi/180,resr);


            unqsch = unique(fd.schnum);
            for schedind = 1:length(unqsch)
                ind0 = fd.schnum==unqsch(schedind);
                ind = ind0;
                for rowind = 1:19
                    ind = ind0 & fd.rowind==rowind;

                    if ~isempty(find(ind))
                        mirror = rps_fit_mirror(structcut(fd,ind),rpsopt,p,'');
                        mirrparms(end+1,:) = [mirror.tilt,mirror.roll];
                        nchans(end+1) = length(find(ind));
                        times(end+1) = nanmean(fd.t(ind));
                        %fdval(end+1) = -nanmean(fd.dk_cen(ind));
                        fdval(end+1) = targind;
                    end
                end
            end
    end
end
%% plot both Moon and RPS mirror fits
fig = figure(1);
fig.Position(3:4) = [900 300];
clf; hold on;
%scatter(times-1*floor(times),mirrparms(:,2),14,fdval,'filled')
scatter(datenum(mjd2datestr(times),'yyyy-mmm-dd:HH:MM:SS'),mirrparms(:,2),14,fdval,'filled')
datetick('x','dd-mmm','keeplimits')
grid on
colormap('jet')



%xlim([-100 200])
%ylim([-0.12 -0.04])









