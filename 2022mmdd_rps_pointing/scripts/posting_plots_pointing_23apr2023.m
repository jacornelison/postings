function posting_plots_pointing_23apr2023()
% Hopefully the last time we have to do this.
% We want to show that the beam center residuals have some small degree of
% flucatuation in rotation, but is mostly dominated by a translation. While
% this could be the result of any sort of systematic in the pointing model,
% this is profoundly solved by allowing the mirror roll to be a free
% parameter on ~1-hour timescales.
%
% We have two main cases -- mirror fixed and mirror free.
% After the beam centers are converted to inst.-fixed coordinates, fit an
% angle and/or translation to them. Since we're primarily concered with the
% translation induced by roll, we should be fitting the translation in
% mirror-fixed coordinates actually where x is always in the direction of
% mirror roll and y is always in the direction of mirror tilt.
%
% For plots, we want to see the distribution of best-fit parameters as well
% as the effect on the residuals. So we'll make histograms of the
% parameters and quiver plots of the residuals.

%% Initialization

addpath('z:/dev/')
addpath('z:/dev/rps/')
addpath(genpath('z:/pipeline/'))
figdir = 'C:\Users\James\Documents\GitHub\postings\2022mmdd_rps_pointing\figs';

%% Load relevant variables
load('z:/dev/rps/rps_phi_final_2022.mat') % RPS 2022 analysis data products
load('z:/dev/rps/rps_obs_info.mat') % which sch goes to which dk/obs
load('z:/dev/rps/pm.mat') % starpointing params
load('z:/dev/rps/fpu_data_obs.mat') % pointing from 2022

%% Set up the cases
mirr_fit = {'free','fixed'};
res_fit = {...
    [0 0 0 0],... No fits
    [1 0 0 0],... Fit Angle only
    [0 0 1 0],... Fit Trans only
    [1 0 1 0]... Fit Both simultaneously
    };

res_ttls = {...
    'None',...
    'Angle Only',...
    'Xtrans Only',...
    'Ang and Xtrans',...
    };

res_names = {...
    'none',...
    'ang',...
    'xtrans',...
    'ang_and_xtrans',...
    };

% Array of per-rasterset fit parameters
unqsch = unique(fd.schnum); unqsch = unqsch(~isnan(unqsch));
fit_array = NaN(...
    length(unqsch)*19,... Nrastersets = Nsch * 19 rastersets
    length(mirr_fit),... N mirror cases
    length(res_fit),... N residual fit cases
    3 ... mirr tilt, FP angle, FP xtrans, x-res
    );

% array of x-residuals
res_array = NaN(...
    length(fd.ch),... Nrastersets = Nsch * 19 rastersets
    length(mirr_fit),... N mirror cases
    length(res_fit),... N residual fit cases
    2 ... x and y
    );

% Init mirror/source
mirror = struct();
mirror.height = 1.4952;
mirror.tilt = 44.88;
mirror.roll = -0.070;

source = struct();
source.distance = 195.5;
source.azimuth = -177.5221;
source.elevation = 2.678;
source.height = source.distance*tand(source.elevation);

rpsopt = struct;
rpsopt.model = model;
rpsopt.source = source;
rpsopt.mirror = mirror;

%% Recompute the mirror tilt/roll per-rasterset

fd0 = fd;
[fd0.mirr_tilt, fd0.mirr_roll] = deal(NaN(size(fd0.ch)));
for schedind = 1:length(unqsch)
    for rowind = 1:19
        idx = logical(fd0.schnum==unqsch(schedind) & fd0.rowind==rowind);
        if ~isempty(find(idx))
            p0 = p;
            fd_rast = structcut(fd0,idx);

            mirror0 = rps_fit_mirror(fd_rast,rpsopt,p,'');
            fd0.mirr_tilt(idx) = mirror0.tilt;
            fd0.mirr_roll(idx) = mirror0.roll;
        end
    end
end
fd = fd0;


%% Compute the pointing et al in one big loop.
% Takes longer but is easier to read.

DOQUIV = false;
KEEPTILT = true;
for mirrind = 1:length(mirr_fit)
    mf = mirr_fit(mirrind);

    for resind = 1:length(res_fit)
        fd0 = fd;
        [fd0.x0, fd0.y0] = deal(NaN(size(fd0.ch)));

        rf = res_fit{resind};

        

        % apply the fixed params if we need to.
        if strcmp(mf,'fixed') & ~KEEPTILT
            [x,y,~] = beam_map_pointing_model(fd0.az_cen,fd0.el_cen,fd0.dk_cen,model,'bicep3',mirror,source,[]);
            fd0.x = reshape(x,size(fd0.ch));
            fd0.y = reshape(y,size(fd0.ch));
            fd0.mirr_tilt = ones(size(fd0.ch))*mirror.tilt;
            fd0.mirr_roll = ones(size(fd0.ch))*mirror.roll;
        end

        % Per-rasterset fits
        % phi_final is already mirror-fit per-rasterset, so there's no point in
        % doing that again when mf==fixed. But we still want to grab the
        % parameters and fit for the angles/translation
        count = 1;
        for schedind = 1:length(unqsch)
            for rowind = 1:19
                idx = logical(fd0.schnum==unqsch(schedind) & fd0.rowind==rowind);
                if ~isempty(find(idx))
                    p0 = p;
                    fd_rast = structcut(fd0,idx);
                    mirror0 = mirror;
                    if strcmp(mf,'fixed') 
                        if KEEPTILT
                        mirror0.tilt = nanmean(fd_rast.mirr_tilt);
                        end
                    else
                       mirror0.tilt = nanmean(fd_rast.mirr_tilt);
                       mirror0.roll = nanmean(fd_rast.mirr_roll);
                    end

                    [x,y,~] = beam_map_pointing_model(fd_rast.az_cen,fd_rast.el_cen,...
                        fd_rast.dk_cen,model,'bicep3',mirror0,source,[]);

                    fd_rast.x = reshape(x,size(fd_rast.ch));
                    fd_rast.y = reshape(y,size(fd_rast.ch));
                                      

                    % Convert to quasi-mirror-centered-coords:
                    % Quasi because we're just subtracting DK
                    [th, r] = cart2pol(fd_rast.x,fd_rast.y);
                    th = th-nanmean(fd_rast.dk_cen)*pi/180;
                    [fd_rast.x, fd_rast.y] = pol2cart(th,r);
                    [fd0.x(idx), fd0.y(idx)] = pol2cart(th,r);

                    p0.theta = p0.theta-nanmean(fd_rast.dk_cen);

                    % Fit xformation params to residuals
                    if resind == 1
                        fpu = struct;
                        fpu.angle = 0;
                        fpu.xtrans = 0;
                    else
                        fpu = fit_fpu_angle_and_scaling_from_xy(fd_rast,p0,'',rf);
                    end

                    % Record all of our fit params for the histograms:
                    fit_array(count,mirrind,resind,:) = [mirror0.roll, fpu.angle, fpu.xtrans];

                    % Now update the fd struct's cmb- and rps-derived
                    % pointings so we can plot the quivers in a bit.
                    param = [fpu.angle, 1, fpu.xtrans, 0];

                    [prx, pry] = pol2cart(p0.theta*pi/180,p0.r);
                    fd_rast.x0 = reshape(prx(fd_rast.ch),size(fd_rast.ch));
                    fd_rast.y0 = reshape(pry(fd_rast.ch),size(fd_rast.ch));

                    xyrot = fpu_ort_and_scaling_model(param,fd_rast);
                    xyrot = reshape(xyrot,[],2);

                    fd0.x0(idx) = reshape(xyrot(:,1),size(fd_rast.ch));
                    fd0.y0(idx) = reshape(xyrot(:,2),size(fd_rast.ch));
                    res_array(idx,mirrind,resind,1) = reshape(fd0.x0(idx)-fd0.x(idx),[],1);
                    res_array(idx,mirrind,resind,2) = reshape(fd0.y0(idx)-fd0.y(idx),[],1);
                end
                count = count+1;
            end
        end

        % Quivers happen per-obs so now let's loop over obs and plot.
        if DOQUIV
            scaling = 20;
            fig = figure(69420);
            fig.Position(3:4) = [500 450]*1.5;
            t = tiledlayout(1,1);
            t.Padding = 'compact';
            t.TileSpacing = 'compact';
            nexttile()

            for obsind = 1:length(dks)
                clf;
                idx = fd0.obsnum==obsind;

                fd_obs = structcut(fd0,idx);
                quiver(fd_obs.x0,fd_obs.y0,(fd_obs.x0-fd_obs.x)*scaling,(fd_obs.y0-fd_obs.y)*scaling,0)
                %plot(fd_obs.x0,fd_obs.y0,'.')
                grid on
                xlabel('X_m [Deg]')
                ylabel('Y_m [Deg]')
                legend(sprintf('CMB - RPS x%i',scaling))
                xlim([-1 1]*15)
                ylim([-1 1]*15)
                pbaspect([1 1 1])
                title({...
                    sprintf('Obs %i | DK %s',obsind, titles{obsind})...
                    sprintf('Mirror: %s | Xform Fit: %s',mirr_fit{mirrind},res_ttls{resind}),...
                    })
                
                fname = sprintf('quiver_finalcheck_mirr_%s_fit_%s_dk_%s',mirr_fit{mirrind},res_names{resind},titles{obsind});
                saveas(fig,fullfile(figdir,fname),'png')

            end
        end

    end
end

%% Now plot the histograms

fig = figure(69421);
fig.Position(3:4) = [1000 300];


lims = {[-1 1]*0.1-0.07 [-1 1]*0.4 [-1 1]*0.1 [-1 1]*0.1};
res = 30;

parmnames = {...
    'roll',...
    'ang',...
    'xtrans',...
    'sres',...
    };

parmttls = {...
    'Mirror Roll [Deg]',...
    'FPU Angle [Deg]',...
    'X_m Translation [Deg]',...
    'X_m Residuals [Deg]',...
    };

for mirrind = 1:length(mirr_fit)
    for resind = 1:length(res_fit)
        t = tiledlayout(1,4);
        t.Padding = 'compact';
        t.TileSpacing = 'compact';

        for parmind = 1:4
            nexttile()
            hold on
            lim = lims{parmind};
            if parmind == 4
            V = squeeze(res_array(:,mirrind,resind,1));
            else
            V = squeeze(fit_array(:,mirrind,resind,parmind));
            V = V(inrange(V,lim(1),lim(2)));
            end
            
            M = nanmean(V);
            S = nanstd(V);
            L = length(V(~isnan(V)));
            if 1
                edges = linspace(lims{parmind}(1),lims{parmind}(2),res);
                N = histc(V,edges);
                bar(edges,N,'histc')
            else
                hist(V,res)
            end
            xlim(lims{parmind})
            grid on
            title({...
                sprintf('M: %1.2E | S: %1.2E | N: %i',M,S,L),...
                })
            xlabel(parmttls{parmind})
            pbaspect([1 1 1])
            
            
        end
        sgtitle(sprintf('Mirror: %s | Xform Fit: %s',mirr_fit{mirrind},res_ttls{resind}))
        fname = sprintf('hist_finalcheck_mirr_%s_fit_%s',mirr_fit{mirrind},res_names{resind});
        saveas(fig,fullfile(figdir,fname),'png')
    end
end

%% avg number of dets per rasterset
Ndet = [];
for schedind = 1:length(unqsch)
    for rowind = 1:19
        idx = logical(fd0.schnum==unqsch(schedind) & fd0.rowind==rowind);
        if ~isempty(find(idx))
            Ndet(end+1) = length(find(idx));
        end
    end
end

%% Make histograms per resolution (WIP)


for resind = 1:4
    fd0 = fd;
    p0 = p;
    % make mirr-fixed obs angles
    r = reshape(p0.r(fd0.ch),size(fd.ch));
    th = reshape(p0.theta(fd0.ch),size(fd.ch))-fd0.dk_cen;
    [fd0.x0, fd0.y0] = pol2cart(th*pi/180,r);

    switch resind
        case 1 % Overall
            fd_rast = fd0;
            mirror0 = mirror;
            [x,y,phi_s] = beam_map_pointing_model(fd_rast.az_cen,fd_rast.el_cen,...
                        fd_rast.dk_cen,model,'bicep3',mirror0,source,[]);
            [th, r] = cart2pol(x,y);

        case 2 % Per Obs

        case 3 % Per Schedule

        case 4 % Per Rasterset

    end
end

%% Bootstrap the mirror fits to get an error

Niter = 100;

mirrparms = NaN(length(unqsch),19,Niter,2);
skiplist = false(10,19);
for iterind = 1:Niter
    
    fd0 = fd;
    [fd0.mirr_tilt, fd0.mirr_roll] = deal(NaN(size(fd0.ch)));
    for schedind = 1:length(unqsch)
        for rowind = 1:19
            idx = find(logical(fd0.schnum==unqsch(schedind) & fd0.rowind==rowind));
            if isempty(idx) %|| skiplist(schedind,rowind)
                continue
            end
            L = length(idx);
            chind = randi(L,L,1);
            p0 = p;
            fd_rast = structcut(fd0,idx(chind));

            mirror0 = rps_fit_mirror(fd_rast,rpsopt,p,'');
            %fd0.mirr_tilt(idx) = mirror0.tilt;
            %fd0.mirr_roll(idx) = mirror0.roll;
            %if mirror0.tilt==0 || mirror.tilt==0
                %skiplist(schedind,rowind) = true;
            %else
            mirrparms(schedind,rowind,iterind,:) = [mirror0.tilt mirror0.roll];
            %end
        end
    end
    
end


save('z:/dev/rps/mirrparms_bootstrap.mat','mirrparms')
%fd = fd0;




