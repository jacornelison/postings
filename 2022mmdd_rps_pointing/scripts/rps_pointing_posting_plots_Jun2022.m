function rps_pointing_posting_plots_Jun2022()
%%
%clear all
%close all
load('z:/dev/rps/rps_beam_fits.mat')
%load('z:/dev/rps/rps_beam_fits_mirror_persch.mat')
load('z:/dev/rps/fpu_data_obs.mat')
load('z:/dev/rps/pm.mat')
load('z:/dev/rps/rpssch.mat')
load('z:/dev/rps/source_fit_data.mat')

addpath('z:/dev')
addpath('z:/pipeline/beammap')
addpath('z:/pipeline/util')

prx = 2*sind(p.r/2).*cosd(p.theta)*180/pi;
pry = 2*sind(p.r/2).*sind(p.theta)*180/pi;

figdir = 'figs/';
%mirror = rpsopt.mirror;
%source = rpsopt.source;
%source.azimuth = -177.5381;
%source.height = 9.2214;

%[fd.x,fd.y,phi] = beam_map_pointing_model(fd.az_cen,fd.el_cen,fd.dk_cen,model,'bicep3',mirror,source,[]);
%fd.x = reshape(fd.x,size(fd.ch));
%fd.y = reshape(fd.y,size(fd.ch));
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

% Schedule groups
scheds = {...
    ...%[1,6],...
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

dks = [61, 23, 174, 68, -81, 90 45 135 68];

titles = {'61', '23','174','68_1','-81','90','45','135','68_2'};


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
    ~ismember(fd.schnum,[1,2,5,6,20]) & ...
    true(size(fd.ch));

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
    saveas(fig,fullfile(figdir,sprintf('source_fit_%s.png',fignames{parmind})))
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

%% Calculate the effective FP orientation angle per DK
% FP orientation angle is the difference between CMB thetas and best-fit
% thetas

[fp_ort_per_dk, phiq_per_dk] = deal(NaN(length(scheds),1));
for schedind = 1:length(scheds)

    [th, r] = cart2pol(fd_per_dk.x(:,schedind),fd_per_dk.y(:,schedind));
    th = rad2deg(th);
    dth = wrapTo180(p.theta-th);
    fp_ort_per_dk(schedind) = nanmedian(dth);
    phiq_per_dk(schedind) = nanmedian(fd_per_dk.phidiff(:,schedind));

end

% Calculate the effective FP orientation angle per sch
unqsch = unique(fd.schnum);
dk_sch = NaN(size(unqsch));


[fp_ort_per_sch, phiq_per_sch] = deal(NaN(length(unqsch),1));
for schedind = 1:length(unqsch)
        [th, r] = cart2pol(fd_per_sch.x(:,schedind),fd_per_sch.y(:,schedind));
        th = rad2deg(th);
        dth = wrapTo180(p.theta-th);
        fp_ort_per_sch(schedind) = nanmedian(dth);
        phiq_per_sch(schedind) = nanmedian(fd_per_sch.phidiff(:,schedind));
        
    for schind = 1:length(scheds)
        if ismember(unqsch(schedind),scheds{schind})
            dk_sch(schedind) = dks(schind);
        end
    end
end

%% Plot the FP ort angle (mean theta's)
fig = figure(1);
fig.Position = [2000 0 400 400];
clf; hold on;


figtitles = {'Mean DTheta'};
fignames = {'theta'};
for parmind = 1
    clf; hold on;
    
    scatter(dksch,fp_ort_per_sch,20,'filled')%,'Marker','^')
    scatter(dks,fp_ort_per_dk,40,'filled','Marker','^')
    
    grid on
    ylabel(sprintf('%s [Degrees]',figtitles{parmind}))
    xlabel('DK [Degrees]')
    %c = colorbar();
    %c.Title.String = 'DK [Deg]';
    %legend({'Fit per sch','Fit per DK'},'Location','northeast')
    saveas(fig,fullfile(figdir,sprintf('%s_fit.png',fignames{parmind})))

end



