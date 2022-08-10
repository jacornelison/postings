function moon_analysis_plots_2022may27()
% For making moon analysis posting plots.
% This was run in Matlab 2022a and is only sparsely commented.
% Contact James Cornelison if you need help.

datadir = 'data/';
%Load the stuff.
load(fullfile(datadir,'moon_beam_fits_phase_corrected.mat'))
%load(('z:/dev/moon_analysis/moon_beam_fits.mat'))
load(fullfile(datadir,'moonsch.mat'))
load(fullfile(datadir,'fpu_data_obs.mat'))
load(fullfile(datadir,'pm.mat'))

inrange = @(A,B,C) A<=C & A>=B;
figdir = 'figs\';

prx = 2 * sind(p.r / 2) .* cosd(p.theta) * 180 / pi;
pry = 2 * sind(p.r / 2) .* sind(p.theta) * 180 / pi;
%[prx, pry] = pol2cart(p.theta*pi/180,p.r);

addpath('z://pipeline/util')
addpath('z://pipeline/beammap')
addpath(genpath('z:dev/'))

% Best fit params for the pointing model.
mirror = struct();
mirror.height = 1.4592;
mirror.tilt = 44.887;
mirror.roll = -0.0697;

source = struct();
source.azimuth = reshape(fd.az_cen_src,[],1);
source.distance = 3.8e8*cosd(reshape(fd.el_cen_src,[],1));
source.height = 3.8e8*sind(reshape(fd.el_cen_src,[],1));

[x, y, phi] = beam_map_pointing_model(fd.az_cen,fd.el_cen,fd.dk_cen,model,'bicep3',mirror,source,[]);
fd.x = x'; fd.y = y';
fd.resx = reshape(prx(fd.ch)-x,size(fd.ch));
fd.resy = reshape(pry(fd.ch)-y,size(fd.ch));
[resth, resr] = cart2pol(fd.resx,fd.resy);
resth = resth - fd.dk_cen*pi/180;
[fd.resx_rot, fd.resy_rot] = pol2cart(resth,resr);
fd.resr = resr;


mount = struct();
mount.aperture_offr = 0;
mount.aperture_offz = 0;
mount.dk_offy = 0;
mount.dk_offx = 0;
mount.el_tilt = 0;
mount.el_offx = 0;
mount.el_offz = 0;
mount.az_tilt_ha = 0;
mount.az_tilt_lat = 0;
mount.az_offz = 0;
mount.az_tilt_org = 0; % 1.22l;
mount.el_tilt_org = 0;
mount.aperture_offz = 0.9970;
mount.az_offz = 1.8355;

% Corresponding schedules
scheds = {...
    [2 3 4],...
    [5 6 7],...
    [8 9 10],...
    [11 12 13],...
    [14 15 18],...
    [16 17 19],...
    [20],...
    [21],...
    [22 23 24],...
    [26 27 28],...
    [30 32 34],...
    [31 33 35],...
    [37 39 41],...
    [38 40 42],...
    [43 45 47],...
    [44 46 48],...
    [49 51],...
    [50 52],...
    [53 54 55],...
    [56 57 58],...
    [59 60 61],...
    [62 63 64],...
    [65 66 67],...
    [68 69 70],...
    [71 72 73],...
    [74 75 76]};

titles = {...
    '0_1',...
    '90_1',...
    '23_1',...
    '174_1',...
    '68_1',...
    '-81_1',...
    '68_2',...
    '68_3',...
    '45_1',...
    '135_1',...
    '23_2',...
    '23_3',...
    '174_2',...
    '174_3',...
    '-81_2',...
    '-81_3',...
    '68_4',...
    '68_5',...
    '0_2',...
    '90_2',...
    '45_2',...
    '135_2',...
    '-23',...
    '-68',...
    '112',...
    '157',...
    };

dks = [0,...
    90,...
    23,...
    174,...
    68,...
    -81,...
    68,...
    68,...
    45,...
    135,...
    23,...
    23,...
    174,...
    174,...
    -81,...
    -81,...
    68,...
    68,...
    0,...
    90,...
    45,...
    135,...
    -23,...
    -68,...
    112,...
    157];

times = [];
for schind = 1:length(scheds)
    si = moonsch{scheds{schind}(1)};
    times(schind) = datenum(mjd2datestr(si.t1),'yyyy-mmm-dd:HH:MM:SS');
end

%%%%%%%%
% Cuts!
%%%%%%%%
cutind = true(size(fd.ch));

% Timing: The as-fit moon position will be WAY off from the Horizons data
%cutind = cutind & inrange(fd.resaz,nanmedian(fd.resaz)-1,nanmedian(fd.resaz)+1);

% Center General Cut
threshold = 0.25;
cutind = cutind & (abs(fd.resx) <= threshold & abs(fd.resy)<=threshold);

% Known shitty schedules cut
%cutind = cutind & fd.schind ~= 55;
fd = structcut(fd,cutind);

cutind = false(size(fd.ch));
% Get rid of stuff not in the scheds list
for schind = 1:length(scheds)
    cutind = cutind | ismember(fd.schind,scheds{schind});
end
fd = structcut(fd,cutind);

% median cut
threshold = 0.5;
cutind = false(size(fd.ch));
for schind = 1:length(moonsch)
    sch = moonsch{schind};
    for scanind = 1:length(sch.scans)
        ind = fd.schind==schind ;%& fd.scanind == scanind;
        medx = nanmedian(fd.resx(ind));
        medy = nanmedian(fd.resy(ind));

        cutind(ind) = inrange(fd.resx(ind),medx-threshold,medx+threshold) &...
            inrange(fd.resy(ind),medy-threshold,medy+threshold);
    end
end
fd = structcut(fd,cutind);

% Cut Derotated residuals
thresh = 0.1;
cutind = false(size(fd.ch));
for schind = 1:length(moonsch)

    for scanind = 1:19
        ind = fd.schind==schind & fd.scanind == scanind;
        if ~isempty(find(ind))
            cutind(ind) =  abs(fd.resx_rot(ind)-nanmedian(fd.resx_rot(ind)))<thresh;
        end
    end
end
fd = structcut(fd,cutind);

% Finally, statistical cuts.
flds = {...
    'gof',...
    };

cutind = true(size(fd.ch));
for fldind = 1:length(flds)
    q = quantile(fd.(flds{fldind}),[0.003, 0.997]);
    cutind = cutind & inrange(fd.(flds{fldind}),q(1),q(2));
end


fd = structcut(fd,cutind);
thresh = 0.1;
cutind = true(size(fd.ch));
for schind = 1:length(moonsch)
    for scanind = 1:19
        ind = find(fd.schind==schind & fd.scanind==scanind);
        if ~isempty(ind)
            q1 = quantile(fd.resx_rot(ind),[thresh, 1-thresh]);
            q2 = quantile(fd.resy_rot(ind),[thresh, 1-thresh]);

            ind2 = inrange(fd.resx_rot(ind),q1(1),q1(2)) &...
                inrange(fd.resy_rot(ind),q2(1),q2(2));
            cutind(ind(~ind2)) = false;

        end
    end
end
%fd = structcut(fd,cutind);

[xps, yps, xs, ys, xrots, yrots] = deal(nan(length(scheds),2640));
for schind = 1:length(scheds)
for chind = 1:2640
    ind = find(fd.ch==chind & ismember(fd.schind,scheds{schind}));
    if ~isempty(ind)
        xps(schind,chind) = mean(fd.x(ind));
        yps(schind,chind) = mean(fd.y(ind));
        xs(schind,chind) = mean(fd.resx(ind));
        ys(schind,chind) = mean(fd.resy(ind));
        xrots(schind,chind) = mean(fd.resx_rot(ind));
        yrots(schind,chind) = mean(fd.resy_rot(ind));

    end
end
end

%%%%%%%%%%%
% REAL Posting plots
%%%%%%%%%%%
%% Quiver plots per "DK"

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

% Things dealing with corrections
corrname = {'','_phase_corrected'};
corrtitle = {'No Correction','Phase Corrected'};
corrparams = {[44.886 -0.0713], [44.882 -0.0751]};


clc
for projind = 1:2
    for fitind = 1:2
        for corrind = 1:2

            %load(sprintf('z:/dev/moon_analysis/perdk_mirror_parms%s.mat',corrname{corrind}))
            load(sprintf('z:/dev/moon_analysis/moon_beam_fits%s_cut.mat',corrname{corrind}))
            [mirrorparms, mirrorerrs, nchans] = deal([]);
            for schedind = 1:length(scheds)
                ind = ismember(fd.schind,scheds{schedind});
            
                if ~isempty(find(ind))
                    fd0 = moon_fit_mirror(structcut(fd,ind),'p',p,'p_ind',p_ind,'savedir','','pm',model);
                    mirrorparms(schedind,:) = fd0.fitparam;
                    mirrorerrs(schedind,:) = fd0.fiterr;
                    nchans(schedind) = length(find(ind));
            
                end
            end
            for schedind = 1:length(scheds)
                ind = ismember(fd.schind,scheds{schedind});

                fd0 = structcut(fd,ind);

                % Only keep channels that has its pair
                if 1
                    ind = true(size(fd0.ch));
                    for chind = 1:length(fd0.ch)
                        cia = find(p_ind.a==fd0.ch(chind));
                        cib = find(p_ind.b==fd0.ch(chind));
                        if ~isempty(cia)
                            idx = find(fd0.ch==p_ind.b(cia) & fd0.schind==fd0.schind(chind) & fd0.scanind==fd0.scanind(chind));

                        elseif ~isempty(cib)
                            idx = find(fd0.ch==p_ind.a(cib) & fd0.schind==fd0.schind(chind) & fd0.scanind==fd0.scanind(chind));
                        end

                        if isempty(idx)
                            ind(chind) = false;
                        end
                    end
                    fd0 = structcut(fd0,ind);
                end

                mirror = struct();
                mirror.height = 1.4592;
                switch fitind
                    case 1
                        mirror.tilt = wmean(mirrorparms(:,1),nchans.^2');
                        mirror.roll = wmean(mirrorparms(:,2),nchans.^2');

                    case 2
                        mirror.tilt = mirrorparms(schedind,1);
                        mirror.roll = mirrorparms(schedind,2);
                end

                source = struct();
                source.azimuth = reshape(fd0.az_cen_src,[],1);
                source.distance = 3.8e8*cosd(reshape(fd0.el_cen_src,[],1));
                source.height = 3.8e8*sind(reshape(fd0.el_cen_src,[],1));



                [x, y, phi] = beam_map_pointing_model(fd0.az_cen,fd0.el_cen,fd0.dk_cen,model,'bicep3',mirror,source,[]);

                prxsch = prx(fd0.ch);
                prysch = pry(fd0.ch);

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
                        [x_track_mirr, y_track_mirr] = get_mirror_coords(fd0.dk_cen,xtrack,ytrack,zeros(size(fd0.ch)),mount,mirror);

                        mk = {'^','+'};
                        for j = 1:length(xtrack)
                            plot(x_track_mirr(j),y_track_mirr(j),'k','MarkerSize',14,'Marker',mk{j})
                        end

                        [x_mirr, y_mirr] = get_mirror_coords(fd0.dk_cen,prxsch',prysch',zeros(size(fd0.ch)),mount,mirror);
                        [x_fit_mirr, y_fit_mirr] = get_mirror_coords(fd0.dk_cen,x',y',zeros(size(fd0.ch)),mount,mirror);
                        x0 = x_mirr;
                        y0 = y_mirr;
                        resx = x_mirr-x_fit_mirr;
                        resy = y_mirr-y_fit_mirr;
                    case 2
                        x0 = prxsch;
                        y0 = prysch;
                        resx = prxsch-x;
                        resy = prysch-y;
                        [resth, resr] = cart2pol(resx,resy);
                        resth = resth - fd0.dk_cen'*pi/180;
                        [resx, resy] = pol2cart(resth,resr);
                end


%                 clf;
%                 quiver(x0,y0,resx*scaling,resy*scaling,0)
%                 hold on;

            if 1
                s = scheds{schedind};
                for si = 1:length(s)
                    for rowind = 1:19
                        ind = fd0.schind==s(si) & fd0.scanind == rowind;
                        quiver(x0(ind),y0(ind),resx(ind)*scaling,resy(ind)*scaling,0,'color',cm(clridx((si-1)*19+rowind),:));
                    end
                end
            end
        

            grid on;
            xlim(xlims{projind})
            ylim(ylims{projind})
            xlabel(sprintf('X%s',projlabels{projind}))
            ylabel(sprintf('Y%s',projlabels{projind}))
                title({sprintf('Beam Center Residuals, x%i, %s', scaling,corrtitle{corrind}),...
                    sprintf('Tilt: %1.2f  Roll: %1.3f  Date: %s',mirror.tilt,mirror.roll,mjd2datestr(moonsch{scheds{schedind}(1)}.t1))...
                    })
                figname = fullfile(figdir,sprintf('quiver_dk%s_fit_%s%s%s.png',titles{schedind},fittype{fitind},corrname{corrind},projnames{projind}));
                saveas(fig,figname)
            end
        end
    end
end

%% Quiver of mean over all dks using mean tilt/roll

% Best fit params for the pointing model.
clc
mirror = struct();
mirror.height = 1.4592;
mirror.tilt = 44.887;
mirror.roll = -0.075;

source = struct();
source.azimuth = reshape(fd.az_cen_src,[],1);
source.distance = 3.8e8*cosd(reshape(fd.el_cen_src,[],1));
source.height = 3.8e8*sind(reshape(fd.el_cen_src,[],1));

[x, y, phi] = beam_map_pointing_model(fd.az_cen,fd.el_cen,fd.dk_cen,model,'bicep3',mirror,source,[]);
fd.resx = reshape(prx(fd.ch)-x,size(fd.ch));
fd.resy = reshape(pry(fd.ch)-y,size(fd.ch));
[resth, resr] = cart2pol(fd.resx,fd.resy);
resth = resth - fd.dk_cen*pi/180;
[fd.resx_rot, fd.resy_rot] = pol2cart(resth,resr);
fd.resr = resr;

moonopt = moon_fit_fpu_angle_and_scaling(fd,mirror,model,p,'');
moonopt.fpu.angle

%%

clc
scaling = 10;


%%
fdsch = struct();
% fdsch.ch = repmat(1:2640,size(xps,1),1);
% fdsch.ch = fdsch.ch(~isnan(xps));
% fdsch.x = xps(~isnan(xps));
% fdsch.y = yps(~isnan(yps));
fdsch.ch = 1:2640';
x0 = nanmean(xps,1)';
y0 = nanmean(yps,1)';
fdsch.ch = fdsch.ch(~isnan(x0)&~isnan(y0));
fdsch.x = x0(~isnan(x0)&~isnan(y0));
fdsch.y = y0(~isnan(y0)&~isnan(y0));

%fpu = fit_fpu_angle_and_scaling_from_xy(fdsch,p)
fpu = fit_fpu_angle_and_scaling_from_xy(fd,p)

%%
model = fpu_ort_and_scaling_model([0.0416 0.9977 0.0030 -0.0017],fdsch);
model = reshape(model,[],2);
x0 = prx(fdsch.ch);
y0 = pry(fdsch.ch);


fig = figure(1);
%quiver(x0,y0,(x0-fdsch.x)*scaling,(y0-fdsch.y)*scaling,0)
quiver(x0,y0,(x0-model(:,1))*scaling,(y0-model(:,2))*scaling,0)
%plot(fdsch.ch,fdsch.x,'.')

%% Bootstrap!
fdsch.ch = 1:2640';
x0 = nanmean(xps,1)';
y0 = nanmean(yps,1)';
fdsch.ch = fdsch.ch(~isnan(x0)&~isnan(y0));
fdsch.x = x0(~isnan(x0)&~isnan(y0));
fdsch.y = y0(~isnan(y0)&~isnan(y0));


iter = 2000;
parms = NaN(4,iter);
for iterind = 1:iter
    bootind = randi([1,length(fdsch.ch)],1,length(fdsch.ch));
    fdsch0.ch = fdsch.ch(bootind);
    fdsch0.x = fdsch.x(bootind);
    fdsch0.y = fdsch.y(bootind);

    fpu = fit_fpu_angle_and_scaling_from_xy(fdsch0,p,'',[1 1 0 0]);
    parms(:,iterind) = [fpu.angle,fpu.scaling,fpu.xtrans,fpu.ytrans];
end


%% Means with rotations applied


[fpuparms, nchans] = deal([]);
for schedind = 1:length(scheds)
    ind = ismember(fd.schind,scheds{schedind});

    if ~isempty(find(ind))
        fpu = fit_fpu_angle_and_scaling_from_xy(structcut(fd,ind),p);
        fpuparms(schedind,:) = [fpu.angle,fpu.scaling,fpu.xtrans,fpu.ytrans];
        nchans(schedind) = length(find(ind));

    end
end


casename = {'none','ang','scale','both'};
parms = [fpu.angle,fpu.scaling,fpu.xtrans,fpu.ytrans];
for caseind = 1:4
    switch caseind
        case 1
            params = [0,1,0,0];
            x = xps;
            y = yps;
            caselabel = sprintf('FPU Rot: %0.2f^o, Scaling: %0.3f',params(1),params(2));
        case 2
            params = parms.*[1,0,0,0];
            params(2) = 1;
            x = params(2).*(xps.*cosd(params(1))-yps.*sind(params(1)));
            y = params(2).*(xps.*sind(params(1))+yps.*cosd(params(1)));
            caselabel = sprintf('FPU Rot: %0.2f^o, Scaling: %0.3f',params(1),params(2));
        case 3
            params = parms.*[0,1,0,0];
            x = params(2).*(xps.*cosd(params(1))-yps.*sind(params(1)));
            y = params(2).*(xps.*sind(params(1))+yps.*cosd(params(1)));
            caselabel = sprintf('FPU Rot: %0.2f^o, Scaling: %0.3f',params(1),params(2));
        case 4
            params = parms;
            x = params(2).*(xps.*cosd(params(1))-yps.*sind(params(1)));
            y = params(2).*(xps.*sind(params(1))+yps.*cosd(params(1)));
            caselabel = sprintf('FPU Rot: %0.2f^o, Scaling: %0.3f',params(1),params(2));
    end

projind = 1;
corrind = 2;
winscale = 1;
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

% Things dealing with corrections
corrname = {'','_phase_corrected'};
corrtitle = {'No Correction','Phase Corrected'};
corrparams = {[44.886 -0.0713], [44.882 -0.0751]};
ind = 1:length(scheds);
quiver(prx,pry,(prx-nanmean(x(ind,:),1)')*scaling,(pry-nanmean(y(ind,:),1)')*scaling,0)
grid on;
xlim(xlims{projind})
ylim(ylims{projind})
xlabel(sprintf('X%s',projlabels{projind}))
ylabel(sprintf('Y%s',projlabels{projind}))
    title({sprintf('Beam Center Residuals, x%i, %s', scaling,corrtitle{corrind}),...
        sprintf('All-DK average, Tilt: %1.2f  Roll: %1.3f',mirror.tilt,mirror.roll),...
        caselabel...
        })
    figname = fullfile(figdir,sprintf('quiver_mean_fit_%s%s%s.png',casename{caseind},corrname{corrind},projnames{projind}));
    saveas(fig,figname)

end



%% calculate the stat and sys uncertainties


clc;
casesch = {[1:6, 8:length(scheds)], [1:6, 8:length(scheds)], [10:length(scheds)]};
timesch = {{[1:6, 8:9], 10:17, 18:25},{[1:6, 8:9], 10:17, 18:25},{10:17, 18:25}};
[mn, sstat, ssys] = deal(NaN(length(casesch),2));

txt{1,1} =  ['|         Uncorrected       |'...
    '          Corrected        |'...
    '      Corrected (Excl.)    |'];
[txt{2,1}, txt{3,1}, txt{4,1}] = deal('|');
txthdr2 = '  mean   |  stat  |   sys  |';
txtrslt = ' %02.4f | %02.4f | %02.4f |';
for caseind = 1:length(casesch)

    txt{2,1} = [txt{2,1} txthdr2];

    switch caseind
        case 1
            %load('z:/dev/moon_analysis/perdk_mirror_parms.mat')
            load('z:/dev/moon_analysis/moon_beam_fits_cut.mat')



        case {2,3}
            %load('z:/dev/moon_analysis/perdk_mirror_parms_phase_corrected.mat')
            load('z:/dev/moon_analysis/moon_beam_fits_phase_corrected_cut.mat')
    end

    [mirrparms, nchans] = deal([]);
    for schedind = 1:length(scheds)
        ind = ismember(fd.schind,scheds{schedind});

        if ~isempty(find(ind))
            fd0 = moon_fit_mirror(structcut(fd,ind),'p',p,'p_ind',p_ind,'savedir','','pm',model);
            mirrparms(schedind,:) = fd0.fitparam;
            mirrerrs(schedind,:) = fd0.fiterr;
            nchans(schedind) = length(find(ind));

        end
    end

    D = [];
    for rowind = 1:(len-1)
        for colind = (rowind+1)
            if dks(rowind)==dks(colind)
                D(end+1,:) = mirrorparms(rowind,:)-mirrorparms(colind,:);
            end
        end
    end
    m = nanmean(mirrparms(casesch{caseind},:),1);
    st = nanmean(mirrerrs(casesch{caseind},:),1);
    [m, wv, neff] = wmean(mirrparms(casesch{caseind},:),nchans(casesch{caseind})',1);
    st = std(abs(D),[],1);

    tsch = timesch{caseind};
    for parmind = 1:2
        [sy,sy1, sy2] = deal([]);
        for timeind = 1:length(tsch)
            sy1(end+1) = nanstd(mirrparms(tsch{timeind},parmind));
            sy2(end+1) = nanmean(mirrparms(tsch{timeind},parmind));
        end
        sy(parmind) = sqrt(sum(sy1.^2)+((max(sy2)-min(sy2))/2).^2);
        txt{2+parmind,1} = [txt{2+parmind,1} sprintf(txtrslt,m(parmind),st(parmind),sy(parmind))];
    end

    mn(caseind,:) = m;
    sstat(caseind,:) = st;
    ssys(caseind,:) = sy;


end

disp(txt)

%% See how each mirror/source param affects residuals

fig = figure(1);
fig.Position(3:4) = [750 675];
clf;

scaling = 10;
lims = [-1 1]*15;

res = 1.5;
az = (-12:res:12)+180;
el = (-12:res:12)+90;
[AZ, EL] = meshgrid(az,el);
AZ = reshape(AZ,[],1);
EL = reshape(EL,[],1);
DK = zeros(size(EL));

model0 = model;
flds = fieldnames(model);
for fldind = 1:length(flds)
    model0.(flds{fldind}) = 0;
end

paramnames = {'tilt','roll','az','el'};
paramlims = {[44.95 45.05],[-0.1 0.1],[-0.1 0.1],[-0.1 0.1]};
paramvars = {'mirror.tilt','mirror.roll','source.azimuth','source.elevation'};
spacing = 10;

% Things dealing with projection
xlims = {[-1 1]*15 [-1 1]*15 [-1 1]*0.5};
ylims = {[-1 1]*15 [-1 1]*15 [-0.5 0.8]};
projlabels = {' [Degrees]','_m [Degrees]','_m [Meters]'};
projnames = {'','_mirror','_mirror'};

% Mount stuff
mount = struct();
mount.aperture_offr = 0;
mount.aperture_offz = 0;
mount.dk_offy = 0;
mount.dk_offx = 0;
mount.el_tilt = 0;
mount.el_offx = 0;
mount.el_offz = 0;
mount.az_tilt_ha = 0;
mount.az_tilt_lat = 0;
mount.az_offz = 0;
mount.az_tilt_org = 0; % 1.22l;
mount.el_tilt_org = 0;
mount.aperture_offz = 0.9970;
mount.az_offz = 1.8355;

cm = colormap('lines');
for projind = 2%1:2
    for paramind = 1:4
        mirror = struct();
        mirror.height = 1.4592;
        mirror.tilt = 45;
        mirror.roll = 0;
        source = struct();
        source.distance = 195.5;
        source.azimuth = 0;
        source.elevation = 0;
        source.height = source.distance.*tand(source.elevation);

        [x0,y0,~] = beam_map_pointing_model(AZ,EL,DK,model0,'bicep3',mirror,source,[]);

        val = linspace(paramlims{paramind}(1),paramlims{paramind}(2),spacing);
        for valind = 1:spacing

            eval(sprintf('%s = val(valind);',paramvars{paramind}))

            source.height = source.distance.*tand(source.elevation);
            [x,y,~] = beam_map_pointing_model(AZ,EL,DK,model0,'bicep3',mirror,source,[]);

            switch projind
                case 1
                    resx = x0-x;
                    resy = y0-y;
                    x1 = x0;
                    x1 = x0;
                case 999

                    xtrack = [1, -1]*10;
                    ytrack = [1, -1]*10;
                    [x_track_mirr, y_track_mirr] = get_mirror_coords(mean(DK)*[1,1],xtrack,ytrack,[0,0],mount,mirror);
                    [x_mirr, y_mirr] = get_mirror_coords(DK,x0,y0,zeros(size(DK)),mount,mirror);
                    [x_fit_mirr, y_fit_mirr] = get_mirror_coords(DK,x,y,zeros(size(DK)),mount,mirror);
                    resx = x_mirr-x_fit_mirr;
                    resy = y_mirr-y_fit_mirr;
                    x1 = x_mirr;
                    y1 = y_mirr;
                    mk = {'^','+'};
                    clf; hold on;
                    for j = 1:length(xtrack)
                        plot(x_track_mirr(j),y_track_mirr(j),'k','MarkerSize',14,'Marker',mk{j})
                    end

                case 2
                        resx = x0-x;
                        resy = y0-y;
                        [resth, resr] = cart2pol(resx,resy);
                        resth = resth - DK*pi/180;
                        [resx, resy] = pol2cart(resth,resr);
                        x1 = x0;
                        y1 = y0;

            end

            quiver(x1,y1,resx*scaling,resy*scaling,0,'Color',cm(1,:))
            grid on
            xlim(xlims{projind})
            ylim(ylims{projind})
            xlabel(sprintf('X%s',projlabels{projind}))
            ylabel(sprintf('Y%s',projlabels{projind}))
            title({sprintf('Best Fit residuals, x%i',scaling),...
                sprintf('%s = %0.2f;',paramvars{paramind},val(valind)),...
                })

            figname = fullfile(figdir,sprintf('testquiver_%s_%i%s.png',paramnames{paramind},valind,projnames{projind}));
            saveas(fig,figname)
        end
    end
end

%% Mirror params overall
fd0 = moon_fit_mirror(fd,'p',p,'p_ind',p_ind,'savedir','','pm',model);
fd0.fitparam

%% fit for mirror per tod.
% How does it look vs time compared to 
clc 
[mirrparms, nchans] = deal([]);
for schedind = 1:length(scheds)
    ind = ismember(fd.schind,scheds{schedind});

    if ~isempty(find(ind))
        fd0 = moon_fit_mirror(structcut(fd,ind),'p',p,'p_ind',p_ind,'savedir','','pm',model);
        mirrparms(schedind,:) = fd0.fitparam;
        nchans(schedind) = length(find(ind));

    end
end


%
fig = figure(1);
fig.Position(3:4) = [900 300];
clf; hold on;
%plot(dks,mirrparms(:,2),'.')
scatter((times-floor(times))*24,mirrparms(:,2),14,times-times(1),'filled')
%plot(times-floor(times),mirrparms(:,2),'.')
grid on
%xlim([-100 200])
ylim([-0.12 -0.04])
xlabel('Time-of-day [hrs]')
ylabel('Mirror Roll [Deg]')
title('Estimated Roll Vs. Time-of-day')
C = colorbar();
C.Label.String = 'Time [Days]';
figname = fullfile(figdir,'roll_vs_tod.png');
saveas(fig,figname)

%% fit for mirror per dk.
% How does it look vs time compared to 
clc 
[mirrparms, nchans] = deal([]);
for schedind = 1:length(scheds)
    ind = ismember(fd.schind,scheds{schedind});

    if ~isempty(find(ind))
        fd0 = moon_fit_mirror(structcut(fd,ind),'p',p,'p_ind',p_ind,'savedir','','pm',model);
        mirrparms(schedind,:) = fd0.fitparam;
        nchans(schedind) = length(find(ind));

    end
end


%
fig = figure(1);
fig.Position(3:4) = [900 300];
clf; hold on;
%plot(dks,mirrparms(:,2),'.')
scatter(dks,mirrparms(:,2),14,times-times(1),'filled')
%plot(times-floor(times),mirrparms(:,2),'.')
grid on
%xlim([-100 200])
ylim([-0.12 -0.04])
xlabel('DK [Deg]')
ylabel('Mirror Roll [Deg]')
title('Estimated Roll Vs. DK')
C = colorbar();
C.Label.String = 'Time [Days]';
figname = fullfile(figdir,'roll_vs_dk.png');
saveas(fig,figname)



%% Look at fpu fits vs DK.
mirror = struct();
mirror.height = 1.4592;
mirror.tilt = 44.883;
mirror.roll = -0.0701;

source = struct();
source.azimuth = reshape(fd.az_cen_src,[],1);
source.distance = 3.8e8*cosd(reshape(fd.el_cen_src,[],1));
source.height = 3.8e8*sind(reshape(fd.el_cen_src,[],1));

[x, y, phi] = beam_map_pointing_model(fd.az_cen,fd.el_cen,fd.dk_cen,model,'bicep3',mirror,source,[]);
fd.x = x'; fd.y = y';
fd.resx = reshape(prx(fd.ch)-x,size(fd.ch));
fd.resy = reshape(pry(fd.ch)-y,size(fd.ch));
[resth, resr] = cart2pol(fd.resx,fd.resy);
resth = resth - fd.dk_cen*pi/180;
[fd.resx_rot, fd.resy_rot] = pol2cart(resth,resr);
fd.resr = resr;


[fpuparms, nchans] = deal([]);
for schedind = 1:length(scheds)
    ind = ismember(fd.schind,scheds{schedind});

    if ~isempty(find(ind))
        fpu = fit_fpu_angle_and_scaling_from_xy(structcut(fd,ind),p);
        fpuparms(schedind,:) = [fpu.angle,fpu.scaling,fpu.xtrans,fpu.ytrans];
        nchans(schedind) = length(find(ind));

    end
end

%
fig = figure(1);
fig.Position(3:4) = [900 300];
clf; hold on;
scatter((times-floor(times))*24,fpuparms(:,1),14,times-times(1),'filled')
%scatter(times,fpuparms(:,1),14,dks,'filled')
grid on
%xlim([-100 200])
xlabel('Time-of-day [hrs]')
ylabel('FPU Rotation Angle [Deg]')
ylim([-0.12 0.4])
title('Estimated FPU orientation Vs. Time-of-day')
C = colorbar();
C.Label.String = 'Time [Days]';
figname = fullfile(figdir,'fpu_ang_vs_tod.png');
saveas(fig,figname)


%% Look at fpu fits vs TOD.
mirror = struct();
mirror.height = 1.4592;
mirror.tilt = 44.883;
mirror.roll = -0.0701;

source = struct();
source.azimuth = reshape(fd.az_cen_src,[],1);
source.distance = 3.8e8*cosd(reshape(fd.el_cen_src,[],1));
source.height = 3.8e8*sind(reshape(fd.el_cen_src,[],1));

[x, y, phi] = beam_map_pointing_model(fd.az_cen,fd.el_cen,fd.dk_cen,model,'bicep3',mirror,source,[]);
fd.x = x'; fd.y = y';
fd.resx = reshape(prx(fd.ch)-x,size(fd.ch));
fd.resy = reshape(pry(fd.ch)-y,size(fd.ch));
[resth, resr] = cart2pol(fd.resx,fd.resy);
resth = resth - fd.dk_cen*pi/180;
[fd.resx_rot, fd.resy_rot] = pol2cart(resth,resr);
fd.resr = resr;


[fpuparms, nchans] = deal([]);
for schedind = 1:length(scheds)
    ind = ismember(fd.schind,scheds{schedind});

    if ~isempty(find(ind))
        fpu = fit_fpu_angle_and_scaling_from_xy(structcut(fd,ind),p);
        fpuparms(schedind,:) = [fpu.angle,fpu.scaling,fpu.xtrans,fpu.ytrans];
        nchans(schedind) = length(find(ind));

    end
end

%
fig = figure(1);
fig.Position(3:4) = [900 300];
clf; hold on;
scatter(dks,fpuparms(:,1),14,times-times(1),'filled')
%scatter(times,fpuparms(:,1),14,dks,'filled')
grid on
xlim([-100 200])
xlabel('DK [Deg]')
ylabel('FPU Rotation Angle [Deg]')
ylim([-0.12 0.4])
title('Estimated FPU orientation Vs. DK')
C = colorbar();
C.Label.String = 'Time [Days]';
colormap('jet')
figname = fullfile(figdir,'fpu_ang_vs_dk.png');
saveas(fig,figname)


%% Corner plots


fig = figure(50);
clf;
edges = (-1:0.05:1)*0.5;
t = tiledlayout(5,5);
len = length(scheds);
[M, S, L, F, n] = deal([]);
for rowind = 1:(len)
    for colind = (rowind+1):(len)
        if dks(rowind)==dks(colind)
        V1 = sqrt(xs(rowind,:).^2+ys(rowind,:).^2);
        V2 = sqrt(xs(colind,:).^2+ys(colind,:).^2);
        Vdiff = V1-V2;
        Vdiff = Vdiff(~isnan(Vdiff)&Vdiff~=0);
        N = histc(Vdiff,edges);
        M(end+1) = nanmean(Vdiff);
        S(end+1) = nanstd(Vdiff);
        L(end+1) = length(find(~isnan(Vdiff)));

        n(end+1) = nexttile;
        bar(edges,N,'histc')
        grid on;
        ttl = title({sprintf('DK%s - DK%s',titles{rowind},titles{colind}), ...
            sprintf('M:%0.3f S%0.3f N:%03i',M(end),S(end),L(end))});
        if abs(rowind-colind)==1
            ttl.Color = [0 0.75 0];
        end
        
        pbaspect([1 1 1])
        if length(S)>18
            xlabel('\Deltar [Deg]')
        end
        end
    end
end
t.TileSpacing = 'tight';
t.Padding = 'tight';
figname = fullfile(figdir,'diff_hists.png');
saveas(fig,figname)



%% Stat uncertainty 
% Get the stat uncert on the mirror params by finding the std of estimates
% at similar DKs that are next to each other.






