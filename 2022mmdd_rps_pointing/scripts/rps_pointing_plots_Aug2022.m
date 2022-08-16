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

% Make a median subtraction array.
fd.phi_medsub = fd.phi;
for chind = 1:2640
    ind = fd.ch==chind;
    if ~isempty(find(ind))
        fd.phi_medsub(ind) = fd.phi_medsub(ind)-nanmedian(fd.phi(ind));
    end
end


%%
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

%% Grab the moon-sun angles
% MJD, , ,raapp,decapp
fname = fullfile('z:/dev/sun_check_2022Aug12.txt');
f = fopen(fname);
hd_sun = textscan(f,'%f%s%s%f%f','delimiter',',','HeaderLines',60);
hd_sun{1} = hd_sun{1}-2400000.5;
fclose(f);
%fd.az_cen_sun = interp1(hd_sun{1},hd_sun{4},fd.t_cen);
fd.az_cen_sun = interp1(hd_sun{1},unwrap(hd_sun{4}*pi/180)*180/pi,fd.t);
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

% Moon-derived mirror params -- Use only these and nothing else!
mirror = struct();
mirror.height = 1.4592;
mirror.tilt= 44.88;
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
        fdval(end+1) = -nanmean(fd.dk_cen(ind));
        times(end+1) = nanmean(fd.t(ind));

    end
end

fig = figure(1);
fig.Position(3:4) = [900 300];
clf; hold on;
scatter((times-0*floor(times))*1,fpuparms(:,1),14,fdval,'filled')
%scatter(-fdval,fpuparms(:,1),14,'filled');%times-floor(times),'filled')
grid on
%xlim([-100 200])
ylim([-0.12 0.4])
%C = colorbar();
%C.Label.String = ''
xlabel('DK [Deg]')
ylabel('FPU Rotation Angle [Deg]')
title('Per-Obs estimated Focal Plane angle Vs. DK')
%saveas(fig,'figs\fpu_rot_vs_dk.png')


%% Estimate fpu uncertainty by resampling with replacement
clc

for schedind = 1:length(scheds)
    ind = ismember(fd.schnum,scheds{schedind});
    if ~isempty(find(ind))
        fpu = fit_fpu_angle_and_scaling_from_xy(structcut(fd,ind),p);
        fpuparms(end+1,:) = [fpu.angle,fpu.scaling,fpu.xtrans,fpu.ytrans];
        nchans(end+1) = length(find(ind));
        fdval(end+1) = -nanmean(fd.dk_cen(ind));
        times(end+1) = nanmean(fd.t(ind));

    end
end

parms_wm = [];
len = length(scheds);
iter = len*10;
for iterind = 1:iter
    ind = randi([1, len],size(scheds));
    parms_wm(end+1,:) = wmean(fpuparms(ind,:),repmat(nchans(ind)'.^2,1,4),1);
end
clc
for valind = 1:4
    
fprintf('Mean: %0.3f, STD: %0.3f\n',mean(parms_wm(:,valind)),std(parms_wm(:,valind)))

end
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
scatter((times-floor(times))*24,mirrparms(:,2),14,'filled')
grid on
%xlim([-100 200])
ylim([-0.12 -0.04])
xlabel('Time-of-day [hrs]')
ylabel('Mirror Roll [Deg]')
title({'Per-Obs Roll Vs. Time-of-day','Type 5 Observations'})

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
            fdval(end+1) = -nanmean(fd.dk_cen(ind));
            %times(end+1) = sch{unqsch(schedind)}.scans(sch{unqsch(schedind)}.index(rowind,1)).t1;
            %fdval(end+1) = sch{unqsch(schedind)}.el_ofs(rowind,1);

        end
    end
end

%
fig = figure(1);
fig.Position(3:4) = [900 300];
clf; hold on;
%scatter((times-1*floor(times)),mirrparms(:,2),14,fdval,'filled')
scatter(mirrparms(:,1),mirrparms(:,2),18,(times-1*floor(times))*24,'filled')
%scatter3((times-0*floor(times))*24,mirrparms(:,1),mirrparms(:,2),14,fdval,'filled')
grid on
%xlim([-100 200])
%ylim([-0.12 -0.04])
%xlabel('Time-of-day [hrs]')
xlabel('tilt')
ylabel('roll')
colormap('parula')
c = colorbar();
c.Label.String = 'TOD';

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

%% Look at mirror stuff vs. stuff

%load('z:/dev/rps/rps_beam_fits_type5_rerun_cut.mat')
load('z:/dev/rps/rps_beam_fits_rerun_all_cut.mat')
load('z:/dev/rps/rps_tilt_data_2022.mat')

fd.tilt = polyval([0.2372 -0.0472],interp1(lj_data.time,lj_data.AIN0,fd.t));
fd.az_cen_sun = interp1(hd_sun{1},unwrap(hd_sun{4}*pi/180)*180/pi,fd.t);
fd.el_cen_sun = interp1(hd_sun{1},hd_sun{5},fd.t);
fd.az_cen_moon = interp1(hd_moon{1},unwrap(hd_moon{4}*pi/180)*180/pi,fd.t);
fd.el_cen_moon = interp1(hd_moon{1},hd_moon{5},fd.t);

mirror = struct();
mirror.height = 1.4592;
mirror.tilt= 44.88;
mirror.roll = -0.06;
rpsopt.mirror = mirror;

rpsopt.source.distance = 195.5;
% Fit for the source params given our mirror info:
%source = rps_fit_source(fd,rpsopt,p,'');
%rpsopt.source = source;

% Fit source with only data where sun is far away from source
ind = ~inrange(wrapTo180(fd.az_cen_sun),-100,100);
source = rps_fit_source(structcut(fd,ind),rpsopt,p,'');
rpsopt.source = source;

%
% With new mirror and source parameters, update the pointing.
[fd.x,fd.y,phi] = beam_map_pointing_model(fd.az_cen,fd.el_cen,fd.dk_cen,model,'bicep3',mirror,source,[]);
fd.x = reshape(fd.x,size(fd.ch));
fd.y = reshape(fd.y,size(fd.ch));
fd.resx = reshape(prx(fd.ch),size(fd.ch))-fd.x;
fd.resy = reshape(pry(fd.ch),size(fd.ch))-fd.y;
[resth, resr] = cart2pol(fd.resx,fd.resy);
[fd.resx_rot, fd.resy_rot] = pol2cart(resth-fd.dk_cen*pi/180,resr);
[fd.xm, fd.ym] = pol2cart((reshape(p.theta(fd.ch),size(fd.ch))-fd.dk_cen)*pi/180,reshape(p.r(fd.ch),size(fd.ch)));

%%
clc
unqsch = unique(fd.schnum);
[mirrparms, fpuparms, fpuparmsobs, nchans, sun_els, sun_azs, times, dksch,xs,ys,xms,yms,rs,thetas,tilts,thms, phias, phibs] = deal([]);
% for schedind = 1:length(scheds)
%     ind = ismember(fd.schnum,scheds{schedind});
for schedind = 1:length(unqsch)
    ind0 = fd.schnum==unqsch(schedind);
    ind = ind0;
    for rowind = 1:19
        ind = ind0 & fd.rowind==rowind;
        fd0 = structcut(fd,ind);
        if ~isempty(find(ind))
            mirror = rps_fit_mirror(fd0,rpsopt,p,'');
            mirrparms(end+1,:) = [mirror.tilt,mirror.roll];

            [fd0.x,fd0.y,phi] = beam_map_pointing_model(fd0.az_cen,fd0.el_cen,fd0.dk_cen,model,'bicep3',mirror,source,[]);
            fd0.x = reshape(fd0.x,size(fd0.ch));
            fd0.y = reshape(fd0.y,size(fd0.ch));

            fpu = fit_fpu_angle_and_scaling_from_xy(fd0,p);
            fpuparmsobs(end+1,:) = [fpu.angle,fpu.scaling,fpu.xtrans,fpu.ytrans];
            fpu = fit_fpu_angle_and_scaling_from_xy(structcut(fd,ind),p);
            fpuparms(end+1,:) = [fpu.angle,fpu.scaling,fpu.xtrans,fpu.ytrans];
            nchans(end+1) = length(find(ind));
            sun_els(end+1) = nanmean(fd.el_cen_sun(ind));
            sun_azs(end+1) = wrapTo180(nanmean(fd.az_cen_sun(ind)-source.azimuth));
            times(end+1) = nanmean(fd.t(ind));
            dksch(end+1) = -nanmean(fd.dk_cen(ind));
            xs(end+1) = nanmean(fd.x(ind));
            ys(end+1) = nanmean(fd.y(ind));
            xms(end+1) = nanmean(fd.xm(ind));
            yms(end+1) = nanmean(fd.ym(ind));
            rs(end+1) = nanmean(sqrt(fd.x(ind).^2+fd.y(ind).^2));
            thetas(end+1) = nanmean(atan2(fd.y(ind),fd.x(ind))*180/pi);
            thms(end+1) = nanmean(atan2(fd.y(ind),fd.x(ind))*180/pi-fd.dk_cen(ind));
            tilts(end+1) = nanmean(fd.tilt(ind));
            phias(end+1) = nanmean(fd.phi_medsub(ind & (ismember(fd.ch,p_ind.a) & (p.mce(fd.ch)~=0)') | (ismember(fd.ch,p_ind.b) & (p.mce(fd.ch)==0)')));
            phibs(end+1) = nanmean(fd.phi_medsub(ind & (ismember(fd.ch,p_ind.b) & (p.mce(fd.ch)~=0)') | (ismember(fd.ch,p_ind.a) & (p.mce(fd.ch)==0)')));
            
        end
    end
end

%%
clc
% X-axis
vals = {times-times(1), (times-floor(times))*24, dksch, sun_els, sun_azs, mirrparms(:,1),mirrparms(:,2),xs,ys,xms,yms,rs,wrapTo180(thetas),tilts,wrapTo180(thms),phias, phibs};
valnames = {'t','tod','dk','sun_els','sun_azs','tilt','roll','x','y','xm','ym','r','theta','tiltmeter','thm','phia','phib'};
vallabels = {'Time [Days]', 'Time-of-day [Hrs]','DK [Deg]','Sun Elevation [Deg]','Sun Azimuth [Deg]','Mirror Tilt [Deg]','Mirror Roll [Deg]','x [deg]' ,'y [deg]','x_m [deg]' ,'y_m [deg]',...
    'r [deg]','theta [deg]','tilt [deg]','theta_m [Deg]','\phi_A med-sub [Deg]','\phi_B med-sub [Deg]'};
vallims = {[-1 60], [0 24] [-100 200] [4 26] [-185 185] [-100 100],[44.85 44.91], [-0.12 -0.04]};

% Y-axis
if 0
    parms = {fpuparms(:,1), fpuparms(:,2)};
    parmnames = {'fpu_ang','fpu_scale'};
    parmlabels = {'FPU Angle [Deg]','FPU Scaling [Deg]'};
    parmlims = {[-0.1 0.4], [0.992 1.002]};

elseif 0
    parms = {mirrparms(:,1), mirrparms(:,2)};
    parmnames = {'tilt','roll'};
    parmlabels = {'Tilt [Deg]','Roll [Deg]'};
    parmlims = {[44.85 44.91], [-0.12 -0.02]};

else
    parms = {mirrparms(:,1), mirrparms(:,2), fpuparms(:,1), fpuparms(:,2), fpuparmsobs(:,1), fpuparmsobs(:,2)};
    parmnames = {'tilt','roll','fpu_ang','fpu_scale','fpu_ang_obs','fpu_scale_obs'};
    parmlabels = {'Tilt [Deg]','Roll [Deg]','FPU Angle [Deg]','FPU Scaling [Deg]','FPU Angle Obs [Deg]','FPU Scaling Obs [Deg]'};
    parmlims = {[44.76 44.92], [-0.12 0],[-0.75 0.4], [0.99 1.005],[-0.75 0.4], [0.99 1.005]};
end
% Pager click 3
%fitparms = {parms, {fpuparmsobs(:,1), fpuparmsobs(:,2)}};
fitnames = {'','_perobs'};


% Other stuff
clrs = [ones(1,11) ones(1,8)*2 ones(1,8)*3];



for fitind = 1%1:length(fitnames)
    %parms = fitparms{fitind};

for valind = 1:length(vals)
    for parmind = 1:length(parms)
        fig = figure(1);
        fig.Position(3:4) = [900 300];
        clf; hold on;
        
        if 1
            scatter(vals{valind},parms{parmind},14,floor(times),'filled')
            colormap('jet')
        else
            scatter(vals{valind},parms{parmind},14,1:length(scheds),'filled')
            colormap('jet')
        end
        grid on
        xlabel(vallabels{valind})
        ylabel(parmlabels{parmind})
        %xlim(vallims{valind})
        ylim(parmlims{parmind})

        figname = fullfile(figdir,sprintf('%s_vs_%s%s.png',parmnames{parmind},valnames{valind},fitnames{fitind}));
        saveas(fig,figname)
    end
end
end






