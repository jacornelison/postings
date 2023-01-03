%%
addpath('z:/dev')
addpath('z:/pipeline/beammap')
addpath('z:/pipeline/util')
addpath('z:/dev/rps')

%% Scatter plots mirrorfit per-rasterset / transformation fits per-rasterset
clc
clear all;
load('z:/dev/rps/fpu_data_obs.mat')
load('z:/dev/rps/pm.mat')
load('z:/dev/rps/source_fit_data.mat')
figdir = 'c:/Users/James/Documents/GitHub/postings/2022mmdd_rps_pointing/figs/';
inrange = @(A,B,C) B <= A & A <= C;
outrange = @(A,B,C) A <= B | C <= A;
[prx, pry] = pol2cart(p.theta*pi/180,p.r);

% Pager plots
rpsopt.model = model;

% Target Specific stuff
targnames = {'moon', 'rps','rps11'};

clc;
tic;

mirror = struct();
mirror.height = 1.4592;
mirror.tilt = 44.88;%44.88;
mirror.roll = -0.07;%-0.07;
load('z:/pipeline/beammap/viridis_cm.mat')
% Loop over targets
[fa, nchans] = deal([]);
[prx, pry] = pol2cart(p.theta*pi/180,p.r);
fpparms = {};
for targind = 1:3
    % Load Moon or RPS data
    switch targnames{targind}
        case 'moon'
            load('z:/dev/rps/moon_beam_fits_phase_corrected_cut_mirror_refit.mat')
            px = reshape(prx(fd.ch),size(fd.ch));
            py = reshape(pry(fd.ch),size(fd.ch));
            % Cut some wonky fits real quick:
            cutind = abs(px-fd.x)*20<1 & abs(py-fd.y)*20<1 & fd.schind ~= 29;% & ismember(fd.schind,1:29);
            fd = structcut(fd,cutind);

            fd.schnum = fd.schind;
            fd.rowind = fd.scanind;
            fd.t = fd.t_cen;

            source = struct();
            source.azimuth = reshape(fd.az_cen_src,[],1);
            source.distance = 3.8e8*cosd(reshape(fd.el_cen_src,[],1));
            source.height = 3.8e8*sind(reshape(fd.el_cen_src,[],1));
            source.elevation = reshape(fd.el_cen_src,[],1);


        case 'rps'
            load('z:/dev/rps/rps_beam_fits_type5_rerun_mirror_refit.mat')


            rpsopt.mirror = mirror;
            rpsopt.source.distance = 195.5;
            % Fit for the source params given our mirror info:
            source = rps_fit_source(fd,rpsopt,p,'');
            rpsopt.source = source;

        case 'rps11'
            load('z:/dev/rps/rps_beam_fits_type11_rerun_mirror_refit.mat')



            rpsopt.mirror = mirror;
            rpsopt.source.distance = 195.5;
            % Fit for the source params given our mirror info:
            source = rps_fit_source(fd,rpsopt,p,'');
            rpsopt.source = source;


    end

    % Resolution of data (Per-obs, or per-scan)
    % Fit the mirror stuff per-rasterset
    % Fit the fpu angle per-obs

    fprintf('%s\n',targnames{targind})
    [mirrparms, fpuparms, fpuparmsobs, sun_els, sun_azs, times, dksch,schedules,az,el] = deal([]);
    for obsind = 1:length(scheds)
        ind = ismember(fd.schnum,scheds{obsind});

        if ~isempty(find(ind))
            fd0 = structcut(fd,ind);

            s = scheds{obsind};
            dummyparms = [];
            schedloop = 1:length(s);
            for schedind = schedloop

                % Fit the mirror params per-rasterset
                for rastind = 1:19
                    idx = fd0.schnum== s(schedind) & fd0.rowind==rastind;

                    if ~isempty(find(idx))
                        fd_rast = structcut(fd0,idx);
                        mirrorperrast = struct();
                        mirrorperrast.height = 1.4592;
                        switch targnames{targind}
                            case 'moon'
                                fd_rast = moon_fit_mirror(fd_rast,'p',p,'p_ind',p_ind,'savedir','','pm',model);
                                mirrparms(end+1,:) = fd_rast.fitparam;
                                mirrorperrast.tilt = fd_rast.fitparam(1);
                                mirrorperrast.roll = fd_rast.fitparam(2);

                                source = struct();
                                source.azimuth = reshape(fd_rast.az_cen_src,[],1);
                                source.distance = 3.8e8*cosd(reshape(fd_rast.el_cen_src,[],1));
                                source.height = 3.8e8*sind(reshape(fd_rast.el_cen_src,[],1));
                                source.elevation = reshape(fd_rast.el_cen_src,[],1);

                            case {'rps','rps11'}
                                mirrorperrast = rps_fit_mirror(fd_rast,rpsopt,p,'');
                                mirrparms(end+1,:) = [mirrorperrast.tilt,mirrorperrast.roll];


                        end

                        [fd_rast.x,fd_rast.y,phi] = beam_map_pointing_model(fd_rast.az_cen,fd_rast.el_cen,fd_rast.dk_cen,model,'bicep3',mirrorperrast,source,[]);
                        fd_rast.x = reshape(fd_rast.x,size(fd_rast.ch));
                        fd_rast.y = reshape(fd_rast.y,size(fd_rast.ch));

                        fd0.x(idx) = fd_rast.x;
                        fd0.y(idx) = fd_rast.y;
                    end
                end
            end



            switch targnames{targind}
                case 'moon'
                    % cut wonky channels in the moon data real quick

                    sun_azs(end+1) = wrapTo180(nanmean(fd0.az_cen_moon-fd0.az_cen_sun));
                    sun_els(end+1) = nanmean(fd0.el_cen_moon-fd0.el_cen_sun);
                case {'rps','rps11'}
                    sun_azs(end+1) = wrapTo180(nanmean(source.azimuth-fd0.az_cen_sun));
                    sun_els(end+1) = nanmean(fd0.el_cen_sun-source.elevation);
            end
            az(end+1) = nanmedian(wrapTo180(fd0.az_cen));
            el(end+1) = nanmedian(fd0.el_cen);
            fpuobs = fit_fpu_angle_and_scaling_from_xy(fd0,p,'',[1,0,0,0]);
            fpuparmsobs(end+1,1) = [fpuobs.angle];
            fpuobs = fit_fpu_angle_and_scaling_from_xy(fd0,p,'',[0,1,1,1]);
            fpuparmsobs(end,2:4) = [fpuobs.scaling,fpuobs.xtrans,fpuobs.ytrans];
%             fpuobs = fit_fpu_angle_and_scaling_from_xy(fd0,p,'',[1,1,1,1]);
%             fpuparmsobs(end+1,:) = [fpuobs.angle,fpuobs.scaling,fpuobs.xtrans,fpuobs.ytrans];
            nchans(end+1) = length(find(ind));
            times(end+1) = nanmean(fd0.t);
            dksch(end+1) = -nanmean(fd0.dk_cen);
            schedules(end+1) = nanmean(fd0.schnum);
            fa(end+1) = fpuobs.angle;
            fd.x(ind) = fd0.x;
            fd.y(ind) = fd0.y;
        end
    end

    %
    % X-axis
    vals = {times-59548, (times-floor(times))*24, dksch, sun_els, sun_azs};
    valnames = {'t','tod','dk','sun_els','sun_azs'};
    vallabels = {'Time [Days]', 'Time-of-day [Hrs]','DK [Deg]','Sun-Target Elevation [Deg]',...
        'Sun-Target Azimuth [Deg]'};
    vallims = {[-1 60], [-0.5 24.5], [-100 200], [-20 20], [-185 185]};

    % Y-Axis
    %[thtrans, rtrans]  = cart2pol(fpuparmsobs(:,3),fpuparmsobs(:,4));
    %thtrans = wrapTo360(thtrans*180/pi);
    parms = {fpuparmsobs(:,1), fpuparmsobs(:,2),fpuparmsobs(:,3), fpuparmsobs(:,4)};%,rtrans,thtrans};
    parmnames = {'fpu_ang_obs','fpu_scale_obs','fpu_xtrans_obs','fpu_ytrans_obs'};%,'fpu_rtrans_obs','fpu_thtrans_obs'};
    parmlabels = {'FPU Angle Obs [Deg]','FPU Scaling Obs','FPU X Trans Obs [Deg]','FPU Y Trans Obs [Deg]'};%,'FPU R Trans Obs [Deg]','FPU Theta Trans Obs [Deg]'};
    parmlims = {[-1 1]*0.01, [0.9995 1.0002],[-1 1]*0.001, [-1 1]*0.001};%, [-1 1]*0.015,[0 360]};
    %parmlims = {[-1 1]*0.2, [0.993 1.002],[-1 1]*0.001, [-1 1]*0.001};%, [-1 1]*0.015,[0 360]};
    fpparms{targind} = [parms{1} parms{2}, parms{3}, parms{4}];
    C_real = cov(fpparms{targind});
    
    if 1
        for valind = 1:length(vals)
            for parmind = 1:length(parms)
                axlimnames = {'','_fixed'};
                for axind = 1:2
                    fig = figure(1);
                    fig.Position(3:4) = [900 300];
                    clf; hold on;

                    if 1
                        scatter(vals{valind},parms{parmind},14,schedules,'filled')
                        colormap(cm)
                    else
                        scatter(vals{valind},parms{parmind},14,1:length(scheds),'filled')
                        colormap(cm)
                    end

                    Soff = sqrt(C_real(parmind,parmind))*ones(size(vallims{valind}));
                    if parmind == 2
                        plot(vallims{valind},1+Soff,'--','Color',[1 1 1]*0.7)
                        plot(vallims{valind},1-Soff,'--','Color',[1 1 1]*0.7)
                    else
                        plot(vallims{valind},Soff,'--','Color',[1 1 1]*0.7)
                        plot(vallims{valind},-Soff,'--','Color',[1 1 1]*0.7)
                    end

                    grid on
                    xlabel(vallabels{valind})
                    ylabel(parmlabels{parmind})
                    xlim(vallims{valind})
                    if axind == 2
                        ylim(parmlims{parmind})
                    end
                    figname = fullfile(figdir,sprintf('%s_vs_%s_%s%s.png',parmnames{parmind},valnames{valind},targnames{targind},axlimnames{axind}));
                    saveas(fig,figname)
                end
            end
        end
    end
end
toc

%% Quiver plots per obs

clc
clear all;
close all;
load('z:/dev/rps/fpu_data_obs.mat')
load('z:/dev/rps/pm.mat')
load('z:/dev/rps/rpssch.mat')
load('z:/dev/rps/source_fit_data.mat')
figdir = 'c:/Users/James/Documents/GitHub/postings/2022mmdd_rps_pointing/figs/';

rpsopt.model = model;

prx = 2*sind(p.r/2).*cosd(p.theta)*180/pi;
pry = 2*sind(p.r/2).*sind(p.theta)*180/pi;
[prx, pry] = pol2cart(p.theta*pi/180,p.r);

winscale = 1.5;
scaling = 20;
fig = figure(1);
fig.Position(3:4) = [500*winscale 450*winscale];
clf;
%cm = colormap('turbo');
load('z:/pipeline/beammap/viridis_cm.mat')
clridx = floor(linspace(1,size(cm,1),3*19));

% Things dealing with projection
xlims = {[-1 1]*15 [-1 1]*15 [-1 1]*0.5};
ylims = {[-1 1]*15 [-1 1]*15 [-0.5 0.8]};
projlabels = {' [Degrees]','_m [Degrees]','_m [Meters]'};
projnames = {'','_mirror','_mirror'};

% Things dealing with fits
fittype = {'overall','perobs','perrast'};

% Things dealing with targets.
targnames = {'moon','rps','rps11'};

% Things to do with angle/scaling/translation correction
corrnames = {'none','ang','scale','xtrans','ytrans','all'};

% Other stuff
mirror = struct();
mirror.height = 1.4592;
mirror.tilt = 44.88;
mirror.roll = -0.07;

% Loop info:
targloop = 1:3;      % Target
projloop = 2;%1:2;      % Projection
fitloop = 3;%1:length(fittype);     % Overall/Per-Obs/Per-rast mirror fitting
%obsloop = 1:length(scheds);  % Observation to plot
%schedloop = 1;%1:3      % Schedule in the observation (if per-rast fit)
rastloop = 1:19;        % Rasterset in the schedule (if per-rast fit)
corrloop = 1:6;         % Angle/scaling/translation correction

scales = [];
schedules = [];
clc
for targind = targloop
    % Load Moon or RPS data
    switch targnames{targind}
        case 'moon'
            %load('z:/dev/rps/moon_beam_fits_phase_corrected_cut.mat')
            load('z:/dev/moon_analysis/moonsch.mat')
            load('z:/dev/rps/moon_beam_fits_phase_corrected_cut_mirror_refit.mat')
            px = reshape(prx(fd.ch),size(fd.ch));
            py = reshape(pry(fd.ch),size(fd.ch));
            % Cut some wonky channels real quick:
            cutind = abs(px-fd.x)*20<1 & abs(py-fd.y)*20<1;
            fd = structcut(fd,cutind);

            sch = moonsch;
            fd.schnum = fd.schind;
            fd.rowind = fd.scanind;
            fd.t = fd.t_cen;
            [fd.phi_medsub, fd.tilt] = deal(NaN(size(fd.ch)));

            source = struct();
            source.azimuth = reshape(fd.az_cen_src,[],1);
            source.distance = 3.8e8*cosd(reshape(fd.el_cen_src,[],1));
            source.height = 3.8e8*sind(reshape(fd.el_cen_src,[],1));
            source.elevation = reshape(fd.el_cen_src,[],1);


        case 'rps'

            load('z:/dev/rps/rps_beam_fits_type5_rerun_cut.mat')
            load('z:/dev/rps/rpssch.mat')
            rpsopt.mirror = mirror;
            rpsopt.source.distance = 195.5;
            % Fit for the source params given our mirror info:
            source = rps_fit_source(fd,rpsopt,p,'');
            %             source.elevation = source.elevation+1.02;
            %             source.height = source.distance.*tand(source.elevation);
            rpsopt.source = source;

        case 'rps11'

            load('z:/dev/rps/rps_beam_fits_type11_rerun_cut.mat')
            load('z:/dev/rps/sch_type11.mat')
            titles = {'0_1','0_2','90_1','90_2','-68','-23'};
            rpsopt.mirror = mirror;
            rpsopt.source.distance = 195.5;
            % Fit for the source params given our mirror info:
            source = rps_fit_source(fd,rpsopt,p,'');
            rpsopt.source = source;

    end

    for projind = projloop
        for fitind = fitloop


            [fd.x,fd.y,phi] = beam_map_pointing_model(fd.az_cen,fd.el_cen,fd.dk_cen,model,'bicep3',mirror,source,[]);
            fd.x = reshape(fd.x,size(fd.ch));
            fd.y = reshape(fd.y,size(fd.ch));

            [mirrparms, mirrerrs, nchans] = deal([]);

            %%%%%%%% Loop over observations %%%%%%%
            obsloop = 1:length(scheds);
            for obsind = obsloop
                ind = ismember(fd.schnum,scheds{obsind});

                fd0 = structcut(fd,ind);

                % Only keep channels that has its pair
                if 1
                    ind = true(size(fd0.ch));
                    for chind = 1:length(fd0.ch)
                        cia = find(p_ind.a==fd0.ch(chind));
                        cib = find(p_ind.b==fd0.ch(chind));
                        if ~isempty(cia)
                            idx = find(fd0.ch==p_ind.b(cia) & fd0.schnum==fd0.schnum(chind) & fd0.rowind==fd0.rowind(chind));

                        elseif ~isempty(cib)
                            idx = find(fd0.ch==p_ind.a(cib) & fd0.schnum==fd0.schnum(chind) & fd0.rowind==fd0.rowind(chind));
                        end

                        if isempty(idx)
                            ind(chind) = false;
                        end
                    end
                    fd0 = structcut(fd0,ind);
                end

                nchans(obsind) = length(find(ind));
                mirrorperobs = struct();
                mirrorperobs.height = 1.4592;
                switch targnames{targind}
                    case 'moon'
                        fd0 = moon_fit_mirror(fd0,'p',p,'p_ind',p_ind,'savedir','','pm',model);
                        mirrparms(end+1,:) = fd0.fitparam;
                        mirrorperobs.tilt = fd0.fitparam(1);
                        mirrorperobs.roll = fd0.fitparam(2);

                        source = struct();
                        source.azimuth = reshape(fd0.az_cen_src,[],1);
                        source.distance = 3.8e8*cosd(reshape(fd0.el_cen_src,[],1));
                        source.height = 3.8e8*sind(reshape(fd0.el_cen_src,[],1));
                        source.elevation = reshape(fd0.el_cen_src,[],1);


                    case {'rps','rps11'}
                        mirrorperobs = rps_fit_mirror(fd0,rpsopt,p,'');
                        mirrparms(end+1,:) = [mirrorperobs.tilt,mirrorperobs.roll];

                end

                if ismember(fitind,[1,2])
                    if fitind==2
                        mirror2 = mirrorperobs;
                    else
                        mirror2 = mirror;
                    end

                    [x, y, phi] = beam_map_pointing_model(fd0.az_cen,fd0.el_cen,fd0.dk_cen,model,'bicep3',mirror2,source,[]);
                    fd0.x = reshape(x,size(fd0.ch));
                    fd0.y = reshape(y,size(fd0.ch));


                elseif fitind == 3
                    % loop over rastersets in a schedule to fit the mirror
                    flds = fieldnames(fd0);
                    for fldind = 1:length(flds)
                        if ~isequal(size(fd0.ch),size(fd0.(flds{fldind})))
                            fd0 = rmfield(fd0,flds{fldind});
                        end
                    end


                    s = scheds{obsind};
                    dummyparms = [];
                    schedloop = 1:length(s);
                    for schedind = schedloop
                        for rastind = rastloop

                            idx = find(fd0.schnum == s(schedind) & fd0.rowind == rastind);
                            if ~isempty(idx)

                                fd_rast = structcut(fd0,idx);
                                mirrorperrast = struct();
                                mirrorperrast.height = 1.4592;
                                switch targnames{targind}
                                    case 'moon'
                                        fd_rast = moon_fit_mirror(fd_rast,'p',p,'p_ind',p_ind,'savedir','','pm',model);

                                        mirrorperrast.tilt = fd_rast.fitparam(1);
                                        mirrorperrast.roll = fd_rast.fitparam(2);
                                        %fd_rast = rmfield(fd_rast,{'guess','fitparam','fitgof','data','model','pm'});

                                        source = struct();
                                        source.azimuth = reshape(fd_rast.az_cen_src,[],1);
                                        source.distance = 3.8e8*cosd(reshape(fd_rast.el_cen_src,[],1));
                                        source.height = 3.8e8*sind(reshape(fd_rast.el_cen_src,[],1));
                                        source.elevation = reshape(fd_rast.el_cen_src,[],1);


                                    case {'rps','rps11'}
                                        mirrorperrast = rps_fit_mirror(fd_rast,rpsopt,p,'');

                                end
                                %mirrorperrast.tilt = mirrorperrast.tilt+0.012/2;
                                dummyparms(end+1,:) = [mirrorperrast.tilt,mirrorperrast.roll];
                                [x, y, phi] = beam_map_pointing_model(fd_rast.az_cen,fd_rast.el_cen,fd_rast.dk_cen,model,'bicep3',mirrorperrast,source,[]);
                                fd_rast.x = reshape(x,size(fd_rast.ch));
                                fd_rast.y = reshape(y,size(fd_rast.ch));

                                fd0.x(idx) = fd_rast.x;
                                fd0.y(idx) = fd_rast.y;


                                clear fd_rast
                            end
                        end
                    end
                    mirror2.tilt = nanmean(dummyparms(:,1));
                    mirror2.roll = nanmean(dummyparms(:,2));

                end



                fpu = fit_fpu_angle_and_scaling_from_xy(fd0,p);

                for corrind = corrloop
                    if ismember(fitind,[1,2,3])


                        switch corrnames{corrind}
                            case 'none'
                                x = fd0.x;
                                y = fd0.y;
                                corrtitle = sprintf('Angle=%0.3f^o Scale=%0.3f Xtrans=%0.3f^o Ytrans=%0.3f^o',0,1,0,0);
                            case 'ang'
                                x = fd0.x.*cosd(fpu.angle)-fd0.y.*sind(fpu.angle);
                                y = fd0.x.*sind(fpu.angle)+fd0.y.*cosd(fpu.angle);
                                corrtitle = sprintf('Angle=%0.3f^o Scale=%0.3f Xtrans=%0.3f^o Ytrans=%0.3f^o',fpu.angle,1,0,0);
                            case 'scale'
                                x = fpu.scaling.*fd0.x;
                                y = fpu.scaling.*fd0.y;
                                corrtitle = sprintf('Angle=%0.3f^o Scale=%0.3f Xtrans=%0.3f^o Ytrans=%0.3f^o',0,fpu.scaling,0,0);
                            case 'xtrans'
                                x = fd0.x+fpu.xtrans;
                                corrtitle = sprintf('Angle=%0.3f^o Scale=%0.3f Xtrans=%0.3f^o Ytrans=%0.3f^o',0,1,fpu.xtrans,0);
                            case 'ytrans'
                                y = fd0.y+fpu.ytrans;
                                corrtitle = sprintf('Angle=%0.3f^o Scale=%0.3f Xtrans=%0.3f^o Ytrans=%0.3f^o',0,1,0,fpu.ytrans);
                            case 'all'
                                x = fd0.x+fpu.xtrans;
                                y = fd0.y+fpu.ytrans;
                                x = fpu.scaling.*(x.*cosd(fpu.angle)-y.*sind(fpu.angle));
                                y = fpu.scaling.*(x.*sind(fpu.angle)+y.*cosd(fpu.angle));
                                corrtitle = sprintf('Angle=%0.3f^o Scale=%0.3f Xtrans=%0.3f^o Ytrans=%0.3f^o',fpu.angle,fpu.scaling,fpu.xtrans,fpu.ytrans);
                        end

                    elseif 0%fitind==3 % Fitting fpu params doesn't work for so few detectors.

                        s = scheds{obsind};
                        dummyparms = [];
                        schedloop = 1:length(s);
                        for schedind = schedloop
                            for rastind = rastloop

                                idx = find(fd0.schnum == s(schedind) & fd0.rowind == rastind);
                                if ~isempty(idx)
                                    fd_rast = structcut(fd0,idx);
                                    fpu = fit_fpu_angle_and_scaling_from_xy(fd_rast,p,'',[1 1 0 0]);
                                    switch corrnames{corrind}
                                        case 'none'
                                            x(idx) = fd_rast.x;
                                            y(idx) = fd_rast.y;
                                            corrtitle = sprintf('Angle=%0.3f^o Scale=%0.3f Xtrans=%0.3f^o Ytrans=%0.3f^o',0,1,0,0);
                                        case 'ang'
                                            x(idx) = fd_rast.x.*cosd(fpu.angle)-fd_rast.y.*sind(fpu.angle);
                                            y(idx) = fd_rast.x.*sind(fpu.angle)+fd_rast.y.*cosd(fpu.angle);
                                            corrtitle = sprintf('Angle=%0.3f^o Scale=%0.3f Xtrans=%0.3f^o Ytrans=%0.3f^o',fpu.angle,1,0,0);
                                        case 'scale'
                                            x(idx) = fpu.scaling.*fd_rast.x;
                                            y(idx) = fpu.scaling.*fd_rast.y;
                                            corrtitle = sprintf('Angle=%0.3f^o Scale=%0.3f Xtrans=%0.3f^o Ytrans=%0.3f^o',0,fpu.scaling,0,0);
                                        case 'xtrans'
                                            x(idx) = fd_rast.x+fpu.xtrans;
                                            corrtitle = sprintf('Angle=%0.3f^o Scale=%0.3f Xtrans=%0.3f^o Ytrans=%0.3f^o',0,1,fpu.xtrans,0);
                                        case 'ytrans'
                                            y(idx) = fd_rast.y+fpu.ytrans;
                                            corrtitle = sprintf('Angle=%0.3f^o Scale=%0.3f Xtrans=%0.3f^o Ytrans=%0.3f^o',0,1,0,fpu.ytrans);
                                        case 'all'
                                            x(idx) = fd_rast.x+fpu.xtrans;
                                            y(idx) = fd_rast.y+fpu.ytrans;
                                            x(idx) = fpu.scaling.*(x(idx).*cosd(fpu.angle)-y(idx).*sind(fpu.angle));
                                            y(idx) = fpu.scaling.*(x(idx).*sind(fpu.angle)+y(idx).*cosd(fpu.angle));
                                            corrtitle = sprintf('Angle=%0.3f^o Scale=%0.3f Xtrans=%0.3f^o Ytrans=%0.3f^o',fpu.angle,fpu.scaling,fpu.xtrans,fpu.ytrans);
                                    end



                                end
                            end
                        end


                    end
                    prxsch = prx(fd0.ch);
                    prysch = pry(fd0.ch);
                    x = reshape(x,size(prxsch));
                    y = reshape(y,size(prysch));

                    clf; hold on;
                    switch projind
                        case 1
                            x0 = prxsch;
                            y0 = prysch;
                            resx = prxsch-x;
                            resy = prysch-y;
                        case 999
                            xtrack = [1, -1]*10;
                            ytrack = [1, -1]*10;
                            [x_track_mirr, y_track_mirr] = get_mirror_coords(fd0.dk_cen,xtrack,ytrack,zeros(size(fd0.ch)),mount,mirror2);

                            mk = {'^','+'};
                            for j = 1:length(xtrack)
                                plot(x_track_mirr(j),y_track_mirr(j),'k','MarkerSize',14,'Marker',mk{j})
                            end

                            [x_mirr, y_mirr] = get_mirror_coords(fd0.dk_cen,prxsch',prysch',zeros(size(fd0.ch)),mount,mirror2);
                            [x_fit_mirr, y_fit_mirr] = get_mirror_coords(fd0.dk_cen,x',y',zeros(size(fd0.ch)),mount,mirror2);
                            x0 = x_mirr;
                            y0 = y_mirr;
                            resx = x_mirr-x_fit_mirr;
                            resy = y_mirr-y_fit_mirr;
                        case 2
                            r0 = p.r(fd0.ch);
                            th0 = (p.theta(fd0.ch) - fd0.dk_cen');
                            x0 = 2 * sind(r0 / 2) .* cosd(th0) * 180 / pi;
                            y0 = 2 * sind(r0 / 2) .* sind(th0) * 180 / pi;
                            resx = prxsch-x;
                            resy = prysch-y;
                            [resth, resr] = cart2pol(resx,resy);
                            resth = resth - fd0.dk_cen'*pi/180;
                            [resx, resy] = pol2cart(resth,resr);

                    end



                    if 1
                        s = scheds{obsind};
                        for si = 1:length(s)
                            for rowind = 1:19
                                ind = fd0.schnum==s(si) & fd0.rowind == rowind;
                                quiver(x0(ind),y0(ind),resx(ind)*scaling,resy(ind)*scaling,0,'color',cm(clridx((si-1)*19+rowind),:));
                            end
                        end
                    end


                    grid on;
                    xlim(xlims{projind})
                    ylim(ylims{projind})
                    xlabel(sprintf('X%s',projlabels{projind}))
                    ylabel(sprintf('Y%s',projlabels{projind}))
                    title({sprintf('Beam Center Residuals, x%i', scaling),...
                        sprintf('Tilt: %1.2f  Roll: %1.3f  Date: %s',mirror2.tilt,mirror2.roll,mjd2datestr(sch{scheds{obsind}(1)}.t1)),...
                        corrtitle})
                    figname = fullfile(figdir,sprintf('quiver_dk%s_fit_%s%s_%s_%s.png',titles{obsind},fittype{fitind},projnames{projind},corrnames{corrind},targnames{targind}));
                    saveas(fig,figname)
                end
            end
        end
    end
end


%% Look at distribution of beam center residuals

clc
%clear all;
load('z:/dev/rps/fpu_data_obs.mat')
load('z:/dev/rps/pm.mat')
load('z:/dev/rps/rpssch.mat')
load('z:/dev/rps/source_fit_data.mat')
load('z:/dev/rps/perdk_fpu_parms.mat')
figdir = 'c:/Users/James/Documents/GitHub/postings/2022mmdd_rps_pointing/figs/';
inrange = @(A,B,C) B <= A & A <= C;
outrange = @(A,B,C) A <= B | C <= A;
[prx, pry] = pol2cart(p.theta*pi/180,p.r);

targnames = {'moon','rps','rps11'};

for targind = 3%1:length(targnames)
    switch targnames{targind}
        case 'moon'
            load('z:/dev/rps/moon_beam_fits_phase_corrected_cut_mirror_refit.mat')

        case 'rps'
            load('z:/dev/rps/rps_beam_fits_type5_rerun_mirror_refit.mat')

        case 'rps11'
            load('z:/dev/rps/rps_beam_fits_type11_rerun_mirror_refit.mat')

    end


    px = reshape(prx(fd.ch),size(fd.ch));
    py = reshape(pry(fd.ch),size(fd.ch));

    fig = figure(1);
    fig.Position(3:4) = [900 300];
    clf; hold on;
    subplot(1,2,1)
    edges = [-1:0.05:1]*0.1;
    %edges = [0:0.075/2:1];
    N = histc((px-fd.x),edges);
    bar(edges,N,'histc')
    grid on
    xlim([min(edges) max(edges)])
    %ylim([0, 5500])
    xlabel('\DeltaX [Deg]')
    ylabel('N')
    title({'X residuals',sprintf('Mean: %0.3f STD: %0.3f', mean(px-fd.x),std(px-fd.x))})



    subplot(1,2,2)
    N = histc(py-fd.y,edges);
    bar(edges,N,'histc')
    grid on
    xlim([min(edges) max(edges)])
    %ylim([0, 5500])
    xlabel('\DeltaY [Deg]')
    title({'Y-residuals',sprintf('Mean: %0.3f STD: %0.3f', mean(py-fd.y),std(py-fd.y))})

    sgtitle('Histogram of Beam Center Residuals')

    saveas(fig,fullfile(figdir,sprintf('beam_cen_hists_%s.png',targnames{targind})))


% Sim the statistical uncertainty of the global parameters

clc
sig = std(py-fd.y);
load('z:/dev/rps/fpu_data_obs.mat')
[prx, pry] = pol2cart(p.theta*pi/180,p.r);
prx = rvec(prx);
pry = rvec(pry);

fpuparms = [];
N = 10000;

for iterind = 1:N

    dx = normrnd(0,sig,size(prx));
    dy = normrnd(0,sig,size(prx));

    fd = struct();
    fd.ch = 1:length(prx);
    fd.x = prx+dx;
    fd.y = pry+dy;

    cutind = ~isnan(fd.x) & ~isnan(fd.y);
    fd = structcut(fd,cutind);

    fpu = fit_fpu_angle_and_scaling_from_xy(fd,p);
    fpuparms = [fpuparms; [fpu.angle fpu.scaling fpu.xtrans fpu.ytrans]];

end

if 0
    figure(1)
    clf; hold on;
    quiver(prx,pry,dx*20,dy*20,0)
    title({'Realization of Beam Center Offsets','Scaling x20, \sigma=0.008^o'})
    grid on
    xlabel('X [Deg]')
    ylabel('Y [Deg]')
end


corrcoef(fpuparms)

end

%% Scatter plots mirrorfit per-rasterset / transformation fits per-rasterset
clc
clear all;
load('z:/dev/rps/fpu_data_obs.mat')
load('z:/dev/rps/pm.mat')
load('z:/dev/rps/source_fit_data.mat')
figdir = 'c:/Users/James/Documents/GitHub/postings/2022mmdd_rps_pointing/figs/';
inrange = @(A,B,C) B <= A & A <= C;
outrange = @(A,B,C) A <= B | C <= A;
[prx, pry] = pol2cart(p.theta*pi/180,p.r);

% Pager plots
rpsopt.model = model;

% Target Specific stuff
targnames = {'moon', 'rps','rps11'};

clc;
tic;

mirror = struct();
mirror.height = 1.4592;
mirror.tilt = 44.88;%44.88;
mirror.roll = -0.07;%-0.07;
load('z:/pipeline/beammap/viridis_cm.mat')
% Loop over targets
[fa, nchans] = deal([]);
[prx, pry] = pol2cart(p.theta*pi/180,p.r);
fpparms = {};
for targind = 1:3

    % Load Moon or RPS data
    switch targnames{targind}
        case 'moon'
            load('z:/dev/rps/moon_beam_fits_phase_corrected_cut_mirror_refit.mat')
            px = reshape(prx(fd.ch),size(fd.ch));
            py = reshape(pry(fd.ch),size(fd.ch));
            % Cut some wonky fits real quick:
            cutind = abs(px-fd.x)*20<1 & abs(py-fd.y)*20<1 & fd.schind ~= 29;% & ismember(fd.schind,1:29);
            fd = structcut(fd,cutind);

            fd.schnum = fd.schind;
            fd.rowind = fd.scanind;
            fd.t = fd.t_cen;

            source = struct();
            source.azimuth = reshape(fd.az_cen_src,[],1);
            source.distance = 3.8e8*cosd(reshape(fd.el_cen_src,[],1));
            source.height = 3.8e8*sind(reshape(fd.el_cen_src,[],1));
            source.elevation = reshape(fd.el_cen_src,[],1);


        case 'rps'
            load('z:/dev/rps/rps_beam_fits_type5_rerun_mirror_refit.mat')


            rpsopt.mirror = mirror;
            rpsopt.source.distance = 195.5;
            % Fit for the source params given our mirror info:
            source = rps_fit_source(fd,rpsopt,p,'');
            rpsopt.source = source;

        case 'rps11'
            load('z:/dev/rps/rps_beam_fits_type11_rerun_mirror_refit.mat')



            rpsopt.mirror = mirror;
            rpsopt.source.distance = 195.5;
            % Fit for the source params given our mirror info:
            source = rps_fit_source(fd,rpsopt,p,'');
            rpsopt.source = source;


    end

    % Resolution of data (Per-obs, or per-scan)
    % Fit the mirror stuff per-rasterset
    % Fit the fpu angle per-obs

    fprintf('%s\n',targnames{targind})
    [mirrparms, fpuparmsobs, sun_els, sun_azs, times, dksch,schedules,az,el] = deal([]);
    [mirrparms0,fpuparmsobs0, sun_els0, sun_azs0, times0, dksch0,schedules0,az0,el0] = deal([]);
    for obsind = 1:length(scheds)
        ind = ismember(fd.schnum,scheds{obsind});
        count = 0;
        if ~isempty(find(ind))
            fd0 = structcut(fd,ind);

            s = scheds{obsind};
            dummyparms = [];
            schedloop = 1:length(s);
            for schedind = schedloop

                % Fit the mirror params per-rasterset
                for rastind = 1:19
                    idx = fd0.schnum== s(schedind) & fd0.rowind==rastind;

                    if ~isempty(find(idx))
                        fd_rast = structcut(fd0,idx);
                        mirrorperrast = struct();
                        mirrorperrast.height = 1.4592;
                        switch targnames{targind}
                            case 'moon'
                                fd_rast = moon_fit_mirror(fd_rast,'p',p,'p_ind',p_ind,'savedir','','pm',model);
                                mirrparms(end+1,:) = fd_rast.fitparam;
                                mirrorperrast.tilt = fd_rast.fitparam(1);
                                mirrorperrast.roll = fd_rast.fitparam(2);

                                source = struct();
                                source.azimuth = reshape(fd_rast.az_cen_src,[],1);
                                source.distance = 3.8e8*cosd(reshape(fd_rast.el_cen_src,[],1));
                                source.height = 3.8e8*sind(reshape(fd_rast.el_cen_src,[],1));
                                source.elevation = reshape(fd_rast.el_cen_src,[],1);

                            case {'rps','rps11'}
                                mirrorperrast = rps_fit_mirror(fd_rast,rpsopt,p,'');
                                mirrparms(end+1,:) = [mirrorperrast.tilt,mirrorperrast.roll];


                        end

                        [fd_rast.x,fd_rast.y,phi] = beam_map_pointing_model(fd_rast.az_cen,fd_rast.el_cen,fd_rast.dk_cen,model,'bicep3',mirrorperrast,source,[]);
                        fd_rast.x = reshape(fd_rast.x,size(fd_rast.ch));
                        fd_rast.y = reshape(fd_rast.y,size(fd_rast.ch));

                        fd0.x(idx) = fd_rast.x;
                        fd0.y(idx) = fd_rast.y;



                        switch targnames{targind}
                            case 'moon'
                                % cut wonky channels in the moon data real quick

                                sun_azs(end+1) = wrapTo180(nanmean(fd_rast.az_cen_moon-fd_rast.az_cen_sun));
                                sun_els(end+1) = nanmean(fd_rast.el_cen_moon-fd_rast.el_cen_sun);
                            case {'rps','rps11'}
                                sun_azs(end+1) = wrapTo180(nanmean(source.azimuth-fd_rast.az_cen_sun));
                                sun_els(end+1) = nanmean(fd_rast.el_cen_sun-source.elevation);
                        end
                        az(end+1) = nanmedian(wrapTo180(fd_rast.az_cen));
                        el(end+1) = nanmedian(fd_rast.el_cen);
                        fpuobs = fit_fpu_angle_and_scaling_from_xy(fd_rast,p,'',[1,0,0,0]);
                        fpuparmsobs(end+1,1) = fpuobs.angle;
%                         fpuobs = fit_fpu_angle_and_scaling_from_xy(fd_rast,p,'',[0,1,0,0]);
%                         fpuparmsobs(end,2) = fpuobs.scaling;
%                         fpuobs = fit_fpu_angle_and_scaling_from_xy(fd_rast,p,'',[0,0,1,0]);
%                         fpuparmsobs(end,3) = fpuobs.xtrans;
%                         fpuobs = fit_fpu_angle_and_scaling_from_xy(fd_rast,p,'',[0,0,0,1]);
%                         fpuparmsobs(end,4) = fpuobs.ytrans;
                        fpuobs = fit_fpu_angle_and_scaling_from_xy(fd_rast,p,'',[0,1,1,1]);
                        fpuparmsobs(end,2:4) = [fpuobs.scaling,fpuobs.xtrans,fpuobs.ytrans];
%                         fpuobs = fit_fpu_angle_and_scaling_from_xy(fd_rast,p,'',[1,1,1,1]);
%                         fpuparmsobs(end+1,:) = [fpuobs.angle,fpuobs.scaling,fpuobs.xtrans,fpuobs.ytrans];
                        nchans(end+1) = length(find(ind));
                        times(end+1) = nanmean(fd_rast.t);
                        dksch(end+1) = -nanmean(fd_rast.dk_cen);
                        schedules(end+1) = nanmean(fd_rast.schnum);
                        %fd.x(ind) = fd_rast.x;
                        %fd.y(ind) = fd_rast.y;
                        count = count+1;
                    end
                end
            end
        end
    end

    %
    % X-axis
    vals = {times-59548, (times-floor(times))*24, dksch, sun_els, sun_azs};
    valnames = {'t','tod','dk','sun_els','sun_azs'};
    vallabels = {'Time [Days]', 'Time-of-day [Hrs]','DK [Deg]','Sun-Target Elevation [Deg]',...
        'Sun-Target Azimuth [Deg]'};
    vallims = {[-1 60], [-0.5 24.5], [-100 200], [-20 20], [-185 185]};

    % Y-Axis
    %[thtrans, rtrans]  = cart2pol(fpuparmsobs(:,3),fpuparmsobs(:,4));
    %thtrans = wrapTo360(thtrans*180/pi);
    parms = {fpuparmsobs(:,1), fpuparmsobs(:,2),fpuparmsobs(:,3), fpuparmsobs(:,4)};%,rtrans,thtrans};
    parmnames = {'fpu_ang_obs','fpu_scale_obs','fpu_xtrans_obs','fpu_ytrans_obs'};%,'fpu_rtrans_obs','fpu_thtrans_obs'};
    parmlabels = {'FPU Angle Obs [Deg]','FPU Scaling Obs','FPU X Trans Obs [Deg]','FPU Y Trans Obs [Deg]'};%,'FPU R Trans Obs [Deg]','FPU Theta Trans Obs [Deg]'};
    parmlims = {[-1 1]*0.2, [0.993 1.002],[-1 1]*0.001, [-1 1]*0.001};%, [-1 1]*0.015,[0 360]};
    %parmlims = {[-1 1]*0.01, [0.9995 1.0002],[-1 1]*0.001, [-1 1]*0.001};%, [-1 1]*0.015,[0 360]};
    fpparms{targind} = [parms{1} parms{2}, parms{3}, parms{4}];
    C_real = cov(fpparms{targind});
    
    if 1
        for valind = 1:length(vals)
            for parmind = 1:length(parms)
                axlimnames = {'','_fixed'};
                for axind = 1:2
                    fig = figure(1);
                    fig.Position(3:4) = [900 300];
                    clf; hold on;

                    if 1
                        scatter(vals{valind},parms{parmind},14,schedules,'filled')
                        colormap(cm)
                    else
                        scatter(vals{valind},parms{parmind},14,1:length(scheds),'filled')
                        colormap(cm)
                    end

                    Soff = sqrt(C_real(parmind,parmind))*ones(size(vallims{valind}));
                    if parmind == 2
                        plot(vallims{valind},1+Soff,'--','Color',[1 1 1]*0.7)
                        plot(vallims{valind},1-Soff,'--','Color',[1 1 1]*0.7)
                    else
                        plot(vallims{valind},Soff,'--','Color',[1 1 1]*0.7)
                        plot(vallims{valind},-Soff,'--','Color',[1 1 1]*0.7)
                    end

                    grid on
                    xlabel(vallabels{valind})
                    ylabel(parmlabels{parmind})
                    xlim(vallims{valind})
                    if axind == 2
                        ylim(parmlims{parmind})
                    end
                    figname = fullfile(figdir,sprintf('%s_vs_%s_%s%s_perrast.png',parmnames{parmind},valnames{valind},targnames{targind},axlimnames{axind}));
                    saveas(fig,figname)
                end
            end
        end
    end
end
toc

