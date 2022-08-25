function rps_pointing_posting_plots_Jun2022()
%%
%clear all
%close all

load('z:/dev/rps/fpu_data_obs.mat')
load('z:/dev/rps/pm.mat')
load('z:/dev/rps/source_fit_data.mat')

addpath('z:/dev')
addpath('z:/pipeline/beammap')
addpath('z:/pipeline/util')


% Schedule groups
schtype = 5;
switch schtype

    case 5
        %load('z:/dev/rps/rps_beam_fits.mat')
        %load('z:/dev/rps/rps_beam_fits_mirror_persch.mat')
        load('z:/dev/rps/rps_beam_fits_type5_rerun.mat')
        load('z:/dev/rps/rpssch.mat')
        load('z:/dev/rps/fpuangle_fit_data.mat')
        load('z:/dev/rps/perdk_fpu_parms.mat')
        scheds = {...
            [1,5],...
            [3,4,7],...
            [8, 9, 10],...
            [11, 12, 13],...
            [14, 15, 16],...
            [17, 18, 19],...
            [21, 22, 23],...
            [24, 25, 26],...
            [27, 28, 29],...
            [30, 31, 32],...
            };

        dks = [0, 61, 23, 174, 68, -81, 90 45 135 68];

        titles = {'0' '61', '23','174','68_1','-81','90','45','135','68_2'};
        disp('type 5')
    case 6
        load('z:/dev/rps/rps_beam_fits_type6_rerun.mat')
        load('z:/dev/rps/sch_type6.mat')
        
        scheds = {...
            [1]...
            };

        dks = [135];

        titles = {'135'};
    
    case 10
        load('z:/dev/rps/rps_beam_fits_type10_rerun.mat')
        load('z:/dev/rps/sch_type10.mat')
        
        scheds = {...
            [1],[2],[3]...
            };

        dks = [0, 90, 45];

        titles = {'0', '90','45'};

    case 11
        load('z:/dev/rps/rps_beam_fits_type11_rerun.mat')
        load('z:/dev/rps/sch_type11.mat')
        %load('z:/dev/rps/fpuangle_fit_data_type11.mat')

        scheds = {...
            [1],[2],[3],[4],[5],[6]...
            };

        dks = [0, 0, 90, 90, -68, -23];

        titles = {'0', '0','90','90','-68','-23'};

        disp('Type 11')
    case 20 % No tilt correction data for DK=90
        load('z:/dev/rps/rps_beam_fits_notilt.mat')
        load('z:/dev/rps/rpssch.mat')
        load('z:/dev/rps/fpuangle_fit_data.mat')
        load('z:/dev/rps/perdk_fpu_parms.mat')
        scheds = {...
            [1,5],...
            [3,4,7],...
            [8, 9, 10],...
            [11, 12, 13],...
            [14, 15, 16],...
            [17, 18, 19],...
            [21, 22, 23],...
            [24, 25, 26],...
            [27, 28, 29],...
            [30, 31, 32],...
            };

        dks = [0,61, 23, 174, 68, -81, 90 45 135 68];

        titles = {'0', '61', '23','174','68_1','-81','90','45','135','68_2'};
        disp('Type5: No tilt')

end

len = length(scheds);
prx = 2*sind(p.r/2).*cosd(p.theta)*180/pi;
pry = 2*sind(p.r/2).*sind(p.theta)*180/pi;

figdir = 'figs/';

% Moon-derived mirror params -- Use only these and nothing else!
mirror = struct();
mirror.height = 1.4592;
mirror.tilt = 44.887;
mirror.roll = -0.0697;


fd.resx = reshape(prx(fd.ch),size(fd.ch))-fd.x;
fd.resy = reshape(pry(fd.ch),size(fd.ch))-fd.y;
[resth, resr] = cart2pol(fd.resx,fd.resy);
[fd.resx_rot, fd.resy_rot] = pol2cart(resth-fd.dk_cen*pi/180,resr);

inrange = @(A,B,C) B <= A & A <= C;
outrange = @(A,B,C) A <= B | C <= A;

fd.derot = fd.phi;
ind = inrange(atand(tand(p.chi(fd.ch)+p.chi_thetaref(fd.ch))),-45,45);
fd.phi_derot(ind) = atand(tand(fd.phi(ind)));
fd.phi_derot(~ind) = atand(tand(fd.phi(~ind)-90))+90;

fd.t = NaN(size(fd.ch));
for fdind = 1:length(fd.ch)
    s = sch{fd.schnum(fdind)};
    index = max(s.index(fd.rowind(fdind),:));
    fd.t(fdind) = s.scans(index).t1;
end


times = NaN(size(scheds));
for schind = 1:length(scheds)
    ind = ismember(fd.schnum,scheds{schind});
    t = nanmean(fd.t(ind));
    times(schind) = datenum(mjd2datestr(t),'yyyy-mmm-dd:HH:MM:SS');
end


% Cuts!
cutind = true(size(fd.ch));
% Obvious outlier cuts
switch schtype
    case 5
        cutind = ...
            inrange(fd.resx,-0.25,0.25) & ...
            inrange(fd.resy,-0.25,0.25) & ...
            inrange(fd.xerr,0,0.02) & ...
            inrange(fd.yerr,0,0.02) & ...
            inrange(fd.phi_err,0,1) & ...
            inrange(fd.agof,0,2e6) & ...
            (inrange(fd.phi_derot,-20,20) | ...
            inrange(fd.phi_derot,70,110)) &...
            ~ismember(fd.schnum,[2,20]) & ...
            fd.nrots>=8 &...
            true(size(fd.ch));
            
        % Special cut for something that went weird in sch 7
            cutind = cutind & ...
                ~(fd.schnum==7 & fd.rowind==17);
    
    case 6
        % Obvious outlier cuts
        cutind = ...
            inrange(fd.resx,-0.25,0.25) & ...
            inrange(fd.resy,-0.25,0.25) & ...
            inrange(fd.xerr,0,0.02) & ...
            inrange(fd.yerr,0,0.02) & ...
            inrange(fd.phi_err,0,1) & ...
            inrange(fd.agof,0,2e6) & ...
            (inrange(fd.phi_derot,-20,20) | ...
            inrange(fd.phi_derot,70,110)) &...
            fd.nrots>=8 &...
            true(size(fd.ch));
    case 10
        % Obvious outlier cuts
        cutind = ...
            inrange(fd.resx,-0.25,0.25) & ...
            inrange(fd.resy,-0.25,0.25) & ...
            inrange(fd.xerr,0,0.02) & ...
            inrange(fd.yerr,0,0.02) & ...
            inrange(fd.phi_err,0,1) & ...
            inrange(fd.agof,0,2e6) & ...
            (inrange(fd.phi_derot,-20,20) | ...
            inrange(fd.phi_derot,70,110)) &...
            fd.nrots>=8 &...
            true(size(fd.ch));
    
    case 11
        % Obvious outlier cuts
        cutind = ...
            inrange(fd.resx,-0.25,0.25) & ...
            inrange(fd.resy,-0.25,0.25) & ...
            inrange(fd.xerr,0,0.02) & ...
            inrange(fd.yerr,0,0.02) & ...
            inrange(fd.phi_err,0,1) & ...
            inrange(fd.agof,0,2e6) & ...
            (inrange(fd.phi_derot,-20,20) | ...
            inrange(fd.phi_derot,70,110)) &...
            fd.nrots>=8 &...
            true(size(fd.ch));
    case 20
        cutind = ...
            inrange(fd.resx,-0.25,0.25) & ...
            inrange(fd.resy,-0.25,0.25) & ...
            inrange(fd.xerr,0,0.02) & ...
            inrange(fd.yerr,0,0.02) & ...
            inrange(fd.phi_err,0,1) & ...
            inrange(fd.agof,0,2e6) & ...
            (inrange(fd.phi_derot,-20,20) | ...
            inrange(fd.phi_derot,70,110)) &...
            ~ismember(fd.schnum,[2,6,20]) & ...
            fd.nrots>=8 &...
            true(size(fd.ch));

        % Special cut for something that went weird in sch 7
            cutind = cutind & ...
                ~(fd.schnum==7 & fd.rowind==17);

end


fd = structcut(fd,cutind);


% Cut a percent of the data statistically
thresh = 2.5; % in percent
thresh = thresh/100;

cutind = true(size(fd.ch));
flds = {'xerr','yerr','phi_err','agof'};%,'resx','resy'};

for valind = 1:length(flds)
    val = fd.(flds{valind});
    q = quantile(val,[0*thresh 1-thresh]);
    cutind = cutind & inrange(val,q(1),q(2));

end

fd = structcut(fd,cutind);


% Which channels are along 0 and 90 (without looking at priors)
inda = inrange(atand(tand(p.chi(fd.ch)+p.chi_thetaref(fd.ch))),-45,45);
indb = ~inda;

fd.phi_medsub = atand(tand(fd.phi));
for chind = 1:length(fd.ch)
    ind = fd.ch==fd.ch(chind);
    if length(find(ind))>1
        if inda(chind)
            fd.phi_medsub(chind) = atand(tand(fd.phi(chind)))-nanmedian(atand(tand(fd.phi(ind))));
        else
            fd.phi_medsub(chind) = atand(tand(fd.phi(chind)-nanmedian(atand(tand(fd.phi(ind)+90)))-90));
        end
    else
        if indb(chind)
            fd.phi_medsub(chind) = atand(tand(fd.phi(chind)-90));
        end
    end
end

cutind = abs(fd.phi_medsub)<1;% & ...
%fd = structcut(fd,cutind);

% Which channels are along 0 and 90 (without looking at priors)
inda = inrange(atand(tand(p.chi(fd.ch)+p.chi_thetaref(fd.ch))),-45,45);
indb = ~inda;


unqsch = unique(fd.schnum);
[dksch, tsch] = deal([]);
for schind = 1:length(unqsch)
    ind = fd.schnum==unqsch(schind);
    tsch(schind) = datenum(mjd2datestr(nanmean(fd.t(ind))),'yyyy-mmm-dd:HH:MM:SS');

    for schind2 = 1:length(scheds)
        if ismember(unqsch(schind),scheds{schind2})
            dksch(schind) = dks(schind2);
        end
    end
end

% find a channel that corresponds to a particular group of schedules
% If there's more than one hit, take the mean.
fd_per_dk = struct();
flds = fieldnames(fd);

for valind = 1:length(flds)
    fd_per_dk.(flds{valind}) = NaN(length(p.gcp),length(scheds));
end

for chind = 1:length(p.gcp)
    for schind = 1:length(scheds)
        ind = ismember(fd.schnum,scheds{schind}) & fd.ch==chind;
        if ~isempty(find(ind))
            for valind = 1:length(flds)
                fd_per_dk.(flds{valind})(chind,schind) = nanmin(fd.(flds{valind})(ind));
            end
        end
    end
end


% Get pair-diff pols

[fd_per_dk.phidiff, fd_per_dk.xpoldiff] = deal(NaN(size(fd_per_dk.schnum)));
for schind = 1:size(fd_per_dk.schnum,2)

    %phi = nanmean(fdps.phi_derot,2);
    %xpol = nanmean(fdps.xpol,2);
    phi = fd_per_dk.phi_derot(:,schind);
    xpol = fd_per_dk.xpol(:,schind);

    Q = (cosd(2*phi(p_ind.a))-xpol(p_ind.a).*cosd(2*phi(p_ind.a)))-(cosd(2*phi(p_ind.b))-xpol(p_ind.b).*cosd(2*phi(p_ind.b)))./...
        (2+xpol(p_ind.a)+xpol(p_ind.b));

    U = (sind(2*phi(p_ind.a))-xpol(p_ind.a).*sind(2*phi(p_ind.a)))-(sind(2*phi(p_ind.b))-xpol(p_ind.b).*sind(2*phi(p_ind.b)))./...
        (2+xpol(p_ind.a)+xpol(p_ind.b));

    phidiff = atan2(U,Q)/2*180/pi;
    xpoldiff = 1-sqrt(Q.^2+U.^2);

    ind = ismember(p_ind.a,find(p.mce==0));
    phidiff(ind) = atand(tand(phidiff(ind)-90));

    fd_per_dk.phidiff(p_ind.a,schind) = phidiff;
    fd_per_dk.xpoldiff(p_ind.a,schind) = xpoldiff;

end

%
% Same shit, but per-schedule
unqsch = unique(fd.schnum);
fd_per_sch = struct();
for valind = 1:length(flds)
    fd_per_sch.(flds{valind}) = NaN(length(p.gcp),length(unqsch));
end

for chind = 1:length(p.gcp)
    for schedind = 1:length(unqsch)
        ind = find(fd.schnum==unqsch(schedind) & fd.ch==chind);
        if ~isempty(ind)
            for valind = 1:length(flds)
                fd_per_sch.(flds{valind})(chind,schedind) = nanmin(fd.(flds{valind})(ind));
            end
        end
    end
end

% Get pair-diff pols
[fd_per_sch.phidiff, fd_per_sch.xpoldiff] = deal(NaN(size(fd_per_sch.schnum)));
for schind = 1:size(fd_per_sch.schnum,2)

    %phi = nanmean(fdps.phi_derot,2);
    %xpol = nanmean(fdps.xpol,2);
    phi = fd_per_sch.phi_derot(:,schind);
    xpol = fd_per_sch.xpol(:,schind);

    Q = (cosd(2*phi(p_ind.a))-xpol(p_ind.a).*cosd(2*phi(p_ind.a)))-(cosd(2*phi(p_ind.b))-xpol(p_ind.b).*cosd(2*phi(p_ind.b)))./...
        (2+xpol(p_ind.a)+xpol(p_ind.b));

    U = (sind(2*phi(p_ind.a))-xpol(p_ind.a).*sind(2*phi(p_ind.a)))-(sind(2*phi(p_ind.b))-xpol(p_ind.b).*sind(2*phi(p_ind.b)))./...
        (2+xpol(p_ind.a)+xpol(p_ind.b));

    phidiff = atan2(U,Q)/2*180/pi;
    xpoldiff = 1-sqrt(Q.^2+U.^2);

    ind = ismember(p_ind.a,find(p.mce==0));
    phidiff(ind) = atand(tand(phidiff(ind)-90));

    fd_per_sch.phidiff(p_ind.a,schind) = phidiff;
    fd_per_sch.xpoldiff(p_ind.a,schind) = xpoldiff;

end


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

save(sprintf('z:/dev/rps/rps_beam_fits_type%01i_rerun_cut.mat',schtype),'fd','scheds','dks','titles');

%% Angle vs time

fig = figure(7);
fig.Position(3:4) = [1000 300];
clf; hold on;
plot(datenum(mjd2datestr(fd.t),'yyyy-mmm-dd:HH:MM:SS'),fd.phi_medsub,'.')
%plot(datenum(mjd2datestr(fd.t),'yyyy-mmm-dd:HH:MM:SS'),fd.resx_rot,'.')
%xlim([738515 738545])
%xticks([738516 738523 738530 738537 738544])
ylim([-1 1])
datetick('x','mmm-dd','keepticks')
grid on
ylabel({'Angle [Degrees]','median subtracted per-channel'})
xlabel('Date [mmm-dd]')
%edges = (-1:0.01:1)*20;
%N = histc(fd.phi_medsub,edges);
%bar(edges,N,'histc')

%% roll-axis vs time

fig = figure(8);
fig.Position = [2000 0 1000 300];
clf; hold on;
scatter(datenum(mjd2datestr(fd.t),'yyyy-mmm-dd:HH:MM:SS'),fd.resx_rot,5,-fd.dk_cen,'filled')
datetick('x','mmm-dd','keepticks')
grid on
ylabel({'Angle [Degrees]','median subtracted per-channel'})
xlabel('Date [mmm-dd]')
c = colorbar();
c.Title.String = 'DK [Deg]';

%% Source fits vs time (per dk)


figtitles = {'Source Az','Source El'};
fignames = {'az','el'};
for parmind = 1:2

    fig = figure(8+parmind);
    fig.Position(3:4) = [1000 300];
    clf; hold on;

    load('z:/dev/rps/persch_source_parms.mat')
    scatter(tsch,sourceparms(:,parmind),14,dksch,'filled')



    load('z:/dev/rps/perdk_source_parms.mat')

    scatter(times,sourceparms(:,parmind),30,dks,'filled','Marker','^')
    datetick('x','mmm-dd','keepticks')
    grid on
    ylabel(sprintf('%s [Degrees]',figtitles{parmind}))
    xlabel('Date [mmm-dd]')
    c = colorbar();
    c.Title.String = 'DK [Deg]';
    legend({'Fit per sch','Fit per DK'},'Location','southwest')
    %saveas(fig,fullfile(figdir,sprintf('source_fit_%s.png',fignames{parmind})))
end


%% Look at the mirror fits

fig = figure(1);
fig.Position = [2000 0 1000 300];
clf; hold on;


figtitles = {'Mirror Tilt','Mirror Roll'};
fignames = {'tilt','roll'};
for parmind = 1:2
    clf; hold on;
    load('z:/dev/rps/persch_mirror_parms.mat')
    scatter(tsch,mirrorparms(:,parmind)',14,dksch,'filled')



    load('z:/dev/rps/perdk_mirror_parms.mat')
    scatter(times,mirrorparms(:,parmind)',30,dks+10,'filled','Marker','^')
    datetick('x','mmm-dd','keeplimits')
    grid on
    ylabel(sprintf('%s [Degrees]',figtitles{parmind}))
    xlabel('Date [mmm-dd]')
    c = colorbar();
    c.Title.String = 'DK [Deg]';
    legend({'Fit per sch','Fit per DK'},'Location','southwest')
    saveas(fig,fullfile(figdir,sprintf('mirror_fit_%s.png',fignames{parmind})))

end

%% Look at the quivers!


fig = figure(1);
fig.Position(3:4) = [750 675];
clf;

scaling = 15;
lims = [-1 1]*15;
pername = {'perdk','persch'};
valname = {'source','mirror'};

load('z:/dev/rps/source_fit_data.mat')
for perind = 1:2
    for valind = 1:2
        for schedind = 1:length(scheds)
            if valind==2
                load('z:/dev/rps/perdk_mirror_parms.mat')

                source = rpsopt.source;
                mirror = rpsopt.mirror;
                ind = ismember(fd.schnum,scheds{schedind});
                fd0 = structcut(fd,ind);

                if perind == 1

                    mirror.tilt = mirrorparms(schedind,1);
                    mirror.roll = mirrorparms(schedind,2);
                    [fd0.x,fd0.y,phi] = beam_map_pointing_model(fd0.az_cen,fd0.el_cen,fd0.dk_cen,model,'bicep3',mirror,source,[]);
                else
                    load('z:/dev/rps/persch_mirror_parms.mat')
                    s = scheds{schedind};
                    for si = 1:length(s)
                        ind = unqsch==s(si);
                        mirror.tilt = mirrorparms(ind,1);
                        mirror.roll = mirrorparms(ind,2);
                        ind = fd0.schnum==s(si);
                        [fd0.x(ind),fd0.y(ind),phi] = beam_map_pointing_model(fd0.az_cen(ind),fd0.el_cen(ind),fd0.dk_cen(ind),model,'bicep3',mirror,source,[]);
                    end
                end




            else
                load('z:/dev/rps/perdk_source_parms.mat')
                mirror = rpsopt.mirror;
                source = rpsopt.source;
                ind = ismember(fd.schnum,scheds{schedind});
                fd0 = structcut(fd,ind);

                if perind == 1

                    source.azimuth = sourceparms(schedind,1);
                    source.elevation = sourceparms(schedind,2);
                    source.height = source.distance*tand(source.elevation);
                    [fd0.x,fd0.y,phi] = beam_map_pointing_model(fd0.az_cen,fd0.el_cen,fd0.dk_cen,model,'bicep3',mirror,source,[]);
                else
                    load('z:/dev/rps/persch_source_parms.mat')
                    s = scheds{schedind};
                    for si = 1:length(s)
                        ind = unqsch==s(si);
                        source.azimuth = sourceparms(ind,1);
                        source.elevation = sourceparms(ind,2);
                        source.height = source.distance*tand(source.elevation);
                        ind = fd0.schnum==s(si);
                        [fd0.x(ind),fd0.y(ind),phi] = beam_map_pointing_model(fd0.az_cen(ind),fd0.el_cen(ind),fd0.dk_cen(ind),model,'bicep3',mirror,source,[]);
                    end
                end

            end

            prx0 = reshape(prx(fd0.ch),size(fd0.ch));
            pry0 = reshape(pry(fd0.ch),size(fd0.ch));

            fd0.x = reshape(fd0.x,size(fd0.ch));
            fd0.y = reshape(fd0.y,size(fd0.ch));
            fd0.resx = reshape(prx(fd0.ch),size(fd0.ch))-fd0.x;
            fd0.resy = reshape(pry(fd0.ch),size(fd0.ch))-fd0.y;
            [resth, resr] = cart2pol(fd0.resx,fd0.resy);
            [fd0.resx_rot, fd0.resy_rot] = pol2cart(resth-fd0.dk_cen*pi/180,resr);


            quiver(prx0,pry0,fd0.resx_rot*scaling,fd0.resy_rot*scaling,0)
            grid on
            xlim(lims)
            ylim(lims)
            xlabel('x [Degrees]')
            ylabel('y [Degrees]')
            title(sprintf('Best Fit residuals, x%i',scaling))

            figname = fullfile(figdir,sprintf('quiver_%s_%s_dk_%s.png',valname{valind},pername{perind},titles{schedind}));
            saveas(fig,figname)
        end
    end
end


%% quiver of means

load('z:/dev/rps/rps_beam_fits_cut.mat')
fd1 = fd;
load('z:/dev/rps/rps_beam_fits_type11_cut.mat')
fd.schnum = fd.schnum+32;
fd = structcat(2,[fd1,fd]);
% Schedule groups
scheds = {...
    ...%[1,6],...
    [3,4,7],...
    ...[5],...
    [8, 9, 10],...
    [11, 12, 13],...
    [14, 15, 16],...
    [17, 18, 19],...
    [21, 22, 23],...
    [24, 25, 26],...
    [27, 28, 29],...
    [30, 31, 32],...
    33,...
    34,...
    35,...
    36,...
    37,...
    38,...
    };

dks = [61, 23, 174, 68, -81, 90 45 135 68 0 0 90 90 -68 -23];

titles = {'61', '23','174','68','-81','90_1','45','135','68_2', '0_1', '0_2','90_2','90_3','-68','-23'};

mirror = struct();
mirror.height = 1.4592;
mirror.tilt= 44.8870;
mirror.roll = -0.0750;

source = struct();
source.distance = 195.5000;
source.azimuth = -177.5247;
source.elevation = 2.6998;
source.height = 9.2189;

ort.angle = 0.0146;%0.0099;
ort.scaling = 0.9966;

fig = figure(1);
fig.Position(3:4) = [750 675];
clf;
lims = [-1 1]*15;
scaling = 10;


for fitind = 1:2
    clf;
    [fd.x,fd.y,~] = beam_map_pointing_model(fd.az_cen,fd.el_cen,fd.dk_cen,model,'bicep3',rpsopt.mirror,rpsopt.source,[]);
    fd.x = reshape(fd.x,size(fd.ch));%*0.996;
    fd.y = reshape(fd.y,size(fd.ch));%*0.996;

    if fitind ==2
        x = ort.scaling*(fd.x.*cosd(ort.angle)-fd.y.*sind(ort.angle));
        y = ort.scaling*(fd.x.*sind(ort.angle)+fd.y.*cosd(ort.angle));
        fd.x = x;
        fd.y = y;
    end

    fd.resx = reshape(prx(fd.ch),size(fd.ch))-fd.x;
    fd.resy = reshape(pry(fd.ch),size(fd.ch))-fd.y;
    [resth, resr] = cart2pol(fd.resx,fd.resy);
    [fd.resx_rot, fd.resy_rot] = pol2cart(resth-fd.dk_cen*pi/180,resr);

    % find a channel that corresponds to a particular group of schedules
    % If there's more than one hit, take the mean.
    fd_per_dk = struct();
    flds = fieldnames(fd);

    for valind = 1:length(flds)
        fd_per_dk.(flds{valind}) = NaN(length(p.gcp),length(scheds));
    end

    for chind = 1:length(p.gcp)
        for schind = 1:length(scheds)
            ind = ismember(fd.schnum,scheds{schind}) & fd.ch==chind;
            if ~isempty(find(ind))
                for valind = 1:length(flds)
                    fd_per_dk.(flds{valind})(chind,schind) = nanmin(fd.(flds{valind})(ind));
                end
            end
        end
    end

    if 1
        xm = wmean(fd_per_dk.resx_rot,fd_per_dk.xerr,2);
        ym = wmean(fd_per_dk.resy_rot,fd_per_dk.yerr,2);
    else
        x = wmean(fd_per_dk.x,fd_per_dk.xerr,2);
        y = wmean(fd_per_dk.y,fd_per_dk.yerr,2);
        xm = reshape(prx,size(x))-x;
        ym = reshape(pry,size(y))-y;
        ind = xm<0.2 & ym<0.2;
        xm = xm(ind);
        ym = ym(ind);
    end


    quiver(prx, pry,xm*scaling,ym*scaling,0)

    %
    fprintf('X M: %0.3f S: %0.3f\n',nanmean(xm),nanstd(xm))
    fprintf('Y M: %0.3f S: %0.3f\n',nanmean(ym),nanstd(ym))
    grid on
    xlim(lims)
    ylim(lims)
    xlabel('x [Degrees]')
    ylabel('y [Degrees]')
    title(sprintf('Best Fit residuals x%i, averaged over all DKs',scaling))

    figname = fullfile(figdir,sprintf('quiver_%s_mean.png',fitnames{fitind}));
    saveas(fig,figname)

end


%% Quivers
% with- and without orientation angle fits.
fig = figure(1);
fig.Position(3:4) = [750 675];
clf;

scaling = 10;

% Things dealing with fits
fitnames = {'none','angfit','scalefit','transfit','all'};
fittitles = {'Angle','Scale','X Trans','Y Trans'};

% Things dealing with projection
xlims = {[-1 1]*15 [-1 1]*0.5};
ylims = {[-1 1]*15 [-0.5 0.8]};
projlabels = {' [Degrees]','_m [Meters]'};
projnames = {'','_mirror'};

load('z:dev/rps/perdk_fpu_parms_type05.mat')
mirror = struct();
mirror.height = 1.4592;
mirror.tilt= 44.8870;
mirror.roll = -0.0697;

source = struct();
source.distance = 195.5000;
source.azimuth = -177.5287;
source.elevation = 2.6990;
source.height = 9.2161;

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

cm = colormap('turbo');
cm = [interp1([1,128,256],[0,0,1]*1,1:256)',...
    interp1([1,128,256],[0,1,0]*0.8,1:256)',...
    interp1([1,128,256],[1,0,0]*0.9,1:256)'];


clridx = floor(linspace(1,size(cm,1),3*19));


for projind = 1:2
    for fitind = 1:5
        for schedind = 1:length(scheds)
            clf; hold on;
            ind = ismember(fd.schnum,scheds{schedind});
            fd0 = structcut(fd,ind);

            [fd0.x,fd0.y,~] = beam_map_pointing_model(fd0.az_cen,fd0.el_cen,fd0.dk_cen,model,'bicep3',mirror,source,[]);

            prx0 = reshape(prx(fd0.ch),size(fd0.ch));
            pry0 = reshape(pry(fd0.ch),size(fd0.ch));

            fd0.x = reshape(fd0.x,size(fd0.ch));
            fd0.y = reshape(fd0.y,size(fd0.ch));
            titlestr = '';
            % Fit cases
            switch fitind
                case 2
                    fd0.x = (fd0.x.*cosd(fpuparms(schedind,1))-fd0.y.*sind(fpuparms(schedind,1)));
                    fd0.y = (fd0.x.*sind(fpuparms(schedind,1))+fd0.y.*cosd(fpuparms(schedind,1)));
                    titlestr = sprintf('%s=%0.3f ',fittitles{1},fpuparms(schedind,1));
                case 3
                    fd0.x = fd0.x*fpuparms(schedind,2);
                    fd0.y = fd0.y*fpuparms(schedind,2);
                    titlestr = sprintf('%s=%0.3f ',fittitles{2},fpuparms(schedind,2));
                case 4
                    fd0.x = fd0.x+fpuparms(schedind,3);
                    fd0.y = fd0.y+fpuparms(schedind,4);
                    titlestr = [sprintf('%s=%0.3f, ',fittitles{3},fpuparms(schedind,3)), ...
                        sprintf('%s=%0.3f ',fittitles{4},fpuparms(schedind,4))];
                case 5
                    fd0.x = fpuparms(schedind,2)*((fd0.x+fpuparms(schedind,3)).*cosd(fpuparms(schedind,1))-(fd0.y+fpuparms(schedind,4)).*sind(fpuparms(schedind,1)));
                    fd0.y = fpuparms(schedind,2)*((fd0.x+fpuparms(schedind,3)).*sind(fpuparms(schedind,1))+(fd0.y+fpuparms(schedind,4)).*cosd(fpuparms(schedind,1)));
                    titlestr = [...
                        sprintf('%s=%0.3f, ',fittitles{1},fpuparms(schedind,1)),...
                        sprintf('%s=%0.3f, ',fittitles{2},fpuparms(schedind,2)),...
                        sprintf('%s=%0.3f, ',fittitles{3},fpuparms(schedind,3)), ...
                        sprintf('%s=%0.3f ',fittitles{4},fpuparms(schedind,4))];
            end

            % Use either focal plane coords or mirror projection.
            switch projind
                case 1


                    fd0.resx = prx0-fd0.x;
                    fd0.resy = pry0-fd0.y;
                case 2

                    xtrack = [1, -1]*10;
                    ytrack = [1, -1]*10;
                    [x_track_mirr, y_track_mirr] = get_mirror_coords(-dks(schedind)*[1,1],xtrack,ytrack,[0,0],mount,mirror);

                    [x_mirr, y_mirr] = get_mirror_coords(fd0.dk_cen,prx0,pry0,zeros(size(fd0.ch)),mount,mirror);
                    [x_fit_mirr, y_fit_mirr] = get_mirror_coords(fd0.dk_cen,fd0.x,fd0.y,zeros(size(fd0.ch)),mount,mirror);
                    fd0.resx = x_mirr-x_fit_mirr;
                    fd0.resy = y_mirr-y_fit_mirr;
                    prx0 = x_mirr;
                    pry0 = y_mirr;
                    mk = {'^','+'};
                    for j = 1:length(xtrack)
                        plot(x_track_mirr(j),y_track_mirr(j),'k','MarkerSize',14,'Marker',mk{j})
                    end
            end

            
            resx = fd0.resx;
            resy = fd0.resy;
            %quiver(prx0,pry0,resx*scaling,resy*scaling,0)
            if 1
                s = scheds{schedind};
                for si = 1:length(s)
                    for rowind = 1:19
                        ind = fd0.schnum==s(si) & fd0.rowind == rowind;
                        quiver(prx0(ind),pry0(ind),resx(ind)*scaling,resy(ind)*scaling,0,'color',cm(clridx((si-1)*19+rowind),:));
                    end
                end
            end
        
            grid on
            xlim(xlims{projind})
            ylim(ylims{projind})
            xlabel(sprintf('X%s',projlabels{projind}))
            ylabel(sprintf('Y%s',projlabels{projind}))
            title({sprintf('Best Fit residuals, x%i',scaling),...
                titlestr...
                })

            figname = fullfile(figdir,sprintf('quiver_%s_dk_%s%s.png',fitnames{fitind},titles{schedind},projnames{projind}))
            saveas(fig,figname)


        end
    end
end
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
xlims = {[-1 1]*15 [-1 1]*0.5};
ylims = {[-1 1]*15 [-0.5 0.8]};
projlabels = {' [Degrees]','_m [Meters]'};
projnames = {'','_mirror'};

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
                case 2

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

%% Make a plot of the best-fit focal plane orientation

fig = figure(1);
fig.Position(3:4) = [750 675];
clf;

scatter(dks,fpuparms(:,1),14,'filled')
grid on
xlabel('DK Angle [Degrees]')
ylabel('Best-fit Orientation Angle [Degrees]')


