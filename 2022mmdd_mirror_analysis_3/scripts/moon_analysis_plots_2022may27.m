function moon_analysis_plots_apr2022()
% For making moon analysis posting plots.
% This was run in Matlab 2019a and is only sparsely commented.
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


% Best fit params for the pointing model.
mirror = struct();
mirror.height = 1.3;
mirror.tilt = 44.887;
mirror.roll = -0.063;

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


% Corresponding schedules
scheds = {...
    [2 3 4],...
    [5 6 7],...
    [8 9 10],...
    [11 12 13],...
    [14 15 18],...
    [16 17 19],...
    [20 21],...
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
    '45_1',...
    '135_1',...
    '23_2',...
    '23_3',...
    '174_2',...
    '174_3',...
    '-81_2',...
    '-81_3',...
    '68_3',...
    '68_4',...
    '0_2',...
    '90_2',...
    '45_2',...
    '135_2',...
    '-23',...
    '-68',...
    '112',...
    '157',...
    };

dks = [0,90,23,174,68,-81,68,45,135,23,23,174,174,-81,-81,68,68,0,90,45,135,-23,-68,112,157];

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
threshold = 3;
cutind = cutind & (abs(fd.resx) <= threshold & abs(fd.resy)<=threshold);


% Known shitty schedules cut
%cutind = cutind & fd.schind ~= 55;
fd = structcut(fd,cutind);

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

fd = structcut(fd,cutind);

%%%%%%%%%%%
% REAL Posting plots
%%%%%%%%%%%
%% First do a scatter hist

titles = {'','Derotated'};
vals = {fd.resx,fd.resy; fd.resx_rot, fd.resy_rot};
for valind = 1:length(titles)
    
    fig = figure(valind);
    fig.Position = [400 0*300 500 450]*1.2;
    clf;
    
    h = scatterhist(vals{valind,1},vals{valind,2},'NBins',[30,50],'Marker','.');
    grid on
    xlim([-1 1]*0.3)
    ylim([-1 1]*0.3)
    title({'Beam Center Residuals [CMB beams-fit]',titles{valind}})
    xlabel('X [Deg]')
    ylabel('Y [Deg]')
    fname = fullfile(figdir,['scatterhist_' titles{valind} '.png']);
    saveas(fig,fname)
end

clear vals



%% Plot and save residual plots
units = {...
    '',...
    '',...
    '',...
    '',...
    '[Deg]',...
    '[Deg]',...
    '[Deg]',...
    '[MJD]',...
    '[Deg]',...
    '[Deg]',...
    '[K]',...
    '[K]',...
    '[m/s]',...
    '[Deg]',...
    '[Deg-C]',...
    };

titles = {...
    'Schedule Number',...
    'Scan Number',...
    'Channel Number',...
    'Goodness-of-Fit (Unnormalized)',...
    'Boresight Az @ Beam Center',...
    'Boresight El @ Beam Center',...
    'Boresight DK @ Beam Center',...
    'Time @ Beam Center',...
    'Moon Az @ Beam Center',...
    'Moon El @ Beam Center',...
    'FPU Temp Mean',...
    'FPU Temp STD',...
    'Wind Speed',...
    'Wind Direction',...
    'Ambient Air Temp',...
    };

cals = [...
    1, 0;...
    1, 0;...
    1, 0;...
    1, 0;...
    1, -model.az_zero;...
    1, -model.el_zero;...
    1, 0;...
    1, 0;...
    1, 0;...
    1, 0;...
    1, 0;...
    1, 0;...
    1, 0;...
    1, -47.35;...
    1, 0;...
    ];


fig = figure(2);
fig.Position = [400 0*300 500 450]*1.2;
clf;
fnames = fieldnames(fd);

resname = {'x','y','r'};
rotname = {'','_derot'};
resval = {fd.resx,fd.resy,fd.resr;...
    fd.resx_rot,fd.resy_rot,fd.resr;};

for rotind = 1:2
    % First in XY
    for valind = 1:length(titles)
        clf;
        val = cals(valind,1)*fd.(fnames{valind})+cals(valind,2);
        q = quantile(val,[0.05 0.95]);
        [b ind] = sort(val);
        ind = ind(end:-1:1);
        pl = scatter(resval{rotind,1}(ind),resval{rotind,2}(ind),4,val(ind),'filled');
        axis equal
        c = colorbar();
        c.Title.String = units{valind};
        xlabel('X-Residual [Deg]')
        ylabel('Y-residual [Deg]')
        title({'Beam Center Residuals [CMB beams-fit]', titles{valind}})
        grid on
        xlim([-1 1]*0.3)
        ylim([-1 1]*0.3)
        colormap jet
        figname = fullfile(figdir,['scatter_' fnames{valind} '_xyres' rotname{rotind} '.png']);
        saveas(fig,figname)
    end
    
    % Next in R/Theta
    for resind = 1:3
        for valind = 1:length(titles)
            clf;
            val = cals(valind,1)*fd.(fnames{valind})+cals(valind,2);
            plot(val,resval{rotind,resind},'.')
            ylabel([resname{resind} '-Residual [Deg]'])
            xlabel([titles{valind} ' ' units{valind}])
            title({'Beam Center Residuals [CMB beams-fit]', titles{valind}})
            grid on
            %xlim([min(val) max(val)]*1.1)
            ylim([-1 1]*0.3)
            figname = fullfile(figdir,['scatter_' fnames{valind} '_' resname{resind} 'res' rotname{rotind} '.png']);
            saveas(fig,figname)
        end
    end
end

%% X Res versus Az

fig = figure(4);
fig.Position = [400 0*300 500 450]*1.2;
clf;

[b ind] = sort(fd.t_cen);
scatter(wrapTo360((fd.az_cen(ind)-model.az_zero)-fd.az_cen_src(ind))-180,fd.resx_rot(ind),4,fd.t_cen(ind)-nanmin(fd.t_cen),'filled')
grid on
xlabel({'Az offset WRT Moon [Deg]','(Boresight Az - Moon Az - 180)'})
ylabel('X-Residual [Deg]')
title({'Beam Center Residuals [CMB beams-fit]','Derotated'})
c = colorbar();
c.Title.String = 'Time [Days]';
ylim([-1 1]*0.15)
fname = fullfile(figdir,'scatter_az_cen_xres_derot_color_t_cen.png');
saveas(fig,fname)
colormap default

%% Y Res vs El
fig = figure(5);
fig.Position = [400 0*300 500 450]*1.2;
clf;

[b ind] = sort(fd.t_cen);
scatter(wrapTo180((fd.el_cen(ind)-model.el_zero)+fd.el_cen_src(ind))-90,fd.resy_rot(ind),4,fd.t_cen(ind)-nanmin(fd.t_cen),'filled')
grid on
xlabel({'El Offset WRT Moon','(Boresight El + Moon El - 90)'})
ylabel('Y-Residual [Deg]')
title({'Beam Center Residuals [CMB beams-fit]','Derotated'})
c = colorbar();
c.Title.String = 'Time [Days]';
ylim([-1 1]*0.3)
fname = fullfile(figdir,'scatter_el_cen_yres_derot_color_t_cen.png');
saveas(fig,fname)
colormap default


%% Look at fits as a function of time difference:

% 2-month time diffs
% Corresponding schedules
scha = [ 2  3  4  7  5  6 22 24 23 26 27 28];
schb = [53 54 55 56 57 58 59 60 61 62 63 64];

val = {fd.resx_rot,fd.resy_rot};

valname = {'xres','yres'};
titles = {'X-residual [Deg]','Y-Residual [Deg]'};
winscale = 1.2;
for valind = 1:2
    fig = figure(2+valind);
    fig.Position = [2000+500*(valind-1)*winscale 0*300 500*winscale 450*winscale];
    clf; hold on;
    pl = [];
    count = 1;
    for schind = 1:length(scha)
        if ismember(schind,[1 4 7 10])
            subplot(2,2,count)
            count = count+1;
            hold on;
        end
        [cha, chb] = deal(NaN(2640,1));
        for chind = 1:2640
            inda = fd.schind==scha(schind) & fd.ch==chind;
            indb = fd.schind==schb(schind) & fd.ch==chind;
            if ~isempty(find(inda)) && ~isempty(find(indb))
                cha(chind) = nanmin(val{valind}(inda));
                chb(chind) = nanmin(val{valind}(indb));
                dk = ceil(nanmedian(fd.dk_cen(inda)));
            end
        end
        size(find(~isnan(cha)));
        pl(schind) = scatter(cha,chb,10,'blue','filled','MarkerFaceAlpha',0.2);
        plot([-1 1],[-1 1],'k--')
        
        
        grid on
        xlim([-1 1]*0.15)
        ylim([-1 1]*0.15)
        xlabel('Schedule on 01 Dec 21')
        ylabel('Schedule on 25 Jan 22')
        title(sprintf('DK: %03i',-dk))
    end
    sgtitle({titles{valind},'Schedules spaced 2 Months apart'})
    figname = fullfile(figdir,['scatter_' valname{valind} '_long_tdiff.png']);
    saveas(fig,figname)
end
%legend(pl([1 4 7 10]),{'0', '90', '45', '135'})

%% 1-Week time diffs

% Corresponding schedules
scha = [30 32 34 37 39 41 43 45 47 49 51];
schb = [31 33 35 38 40 42 44 46 48 50 52];
schc = [0  10  9 11 13 12 16 17 19 14 15];


for midind = 1:2
    % Corresponding schedules
    scha = [30 32 34 37 39 41 43 45 47 49 51];
    schb = [31 33 35 38 40 42 44 46 48 50 52];
    schc = [0  10  9 11 13 12 16 17 19 14 15];
    
    
    if midind==1
        schb = scha;
    end
    scha = schc;
    
    
    val = {fd.resx_rot,fd.resy_rot};
    titles = {'X-residual [Deg]','Y-Residual [Deg]'};
    winscale = 1.2;
    for valind = 1:2
        fig = figure(2+valind);
        fig.Position = [2000+500*(valind-1)*winscale 0*300 500*winscale 450*winscale];
        clf; hold on;
        pl = [];
        count = 1;
        for schind = 1:length(scha)
            if ismember(schind,[1 4 7 10])
                subplot(2,2,count)
                count = count+1;
                hold on;
            end
            [cha, chb] = deal(NaN(2640,1));
            for chind = 1:2640
                inda = fd.schind==scha(schind) & fd.ch==chind;
                indb = fd.schind==schb(schind) & fd.ch==chind;
                if ~isempty(find(inda)) && ~isempty(find(indb))
                    cha(chind) = nanmin(val{valind}(inda));
                    chb(chind) = nanmin(val{valind}(indb));
                    dk = ceil(nanmedian(fd.dk_cen(inda)));
                end
            end
            size(find(~isnan(cha)));
            pl(schind) = scatter(cha,chb,10,'blue','filled','MarkerFaceAlpha',0.2);
            plot([-1 1],[-1 1],'k--')
            
            
            grid on
            xlim([-1 1]*0.15)
            ylim([-1 1]*0.15)
            xlabel('Schedule on 03 Dec 21')
            ylabel('Schedule on 10 Dec 21')
            title(sprintf('DK: %03i',-dk))
        end
        sgtitle({titles{valind},'Schedules spaced 1 Week apart'})
        figname = fullfile(figdir,['scatter_' valname{valind} '_mid' num2str(midind) '_tdiff.png']);
        saveas(fig,figname)
    end
    
    
end




%% Short time diffs

% Corresponding schedules
scha = [30 32 34 37 39 41 43 45 47 49 51];
schb = [31 33 35 38 40 42 44 46 48 50 52];


val = {fd.resx_rot,fd.resy_rot};
titles = {'X-residual [Deg]','Y-Residual [Deg]'};
winscale = 1.2;
for valind = 1:2
    fig = figure(valind);
    fig.Position = [2000+500*(valind-1)*winscale 0*300 500*winscale 450*winscale];
    clf; hold on;
    pl = [];
    count = 1;
    for schind = 1:length(scha)
        if ismember(schind,[1 4 7 10])
            subplot(2,2,count)
            count = count+1;
            hold on;
        end
        [cha, chb] = deal(NaN(2640,1));
        for chind = 1:2640
            inda = fd.schind==scha(schind) & fd.ch==chind;
            indb = fd.schind==schb(schind) & fd.ch==chind;
            if ~isempty(find(inda)) && ~isempty(find(indb))
                cha(chind) = nanmin(val{valind}(inda));
                chb(chind) = nanmin(val{valind}(indb));
                dk = ceil(nanmedian(fd.dk_cen(inda)));
            end
        end
        size(find(~isnan(cha)));
        pl(schind) = scatter(cha,chb,10,'blue','filled','MarkerFaceAlpha',0.2);
        plot([-1 1],[-1 1],'k--')
        
        
        grid on
        xlim([-1 1]*0.15)
        ylim([-1 1]*0.15)
        xlabel('Schedule on 10 Dec 21')
        ylabel('Schedule on 10 Dec 21+1hr')
        title(sprintf('DK: %03i',-dk))
    end
    sgtitle({titles{valind},'Schedules spaced back-to-back'})
    figname = fullfile(figdir,['scatter_' valname{valind} '_short_tdiff.png']);
    saveas(fig,figname)
end


%%
% So far, we've established from time spacing
% tilt = 44.887 pm 0.008;
% roll = -0.063 pm 0.02;
mirror = struct();
mirror.height = 1.3;
mirror.tilt = 44.887+0.001;
mirror.roll = -0.063;

source = struct();
source.azimuth = reshape(fd.az_cen_src,[],1);
source.distance = 3.8e8*cosd(reshape(fd.el_cen_src,[],1));
source.height = 3.8e8*sind(reshape(fd.el_cen_src,[],1));

[x1, y1, phi1] = beam_map_pointing_model(fd.az_cen,fd.el_cen,fd.dk_cen,model,'bicep3',mirror,source,[]);

mirror.tilt = 44.887-0.001;
[x2, y2, phi2] = beam_map_pointing_model(fd.az_cen,fd.el_cen,fd.dk_cen,model,'bicep3',mirror,source,[]);

nanmean(phi2-phi1)/2

%% Look diffs from DK to DK


[res_mat, res_mat_rot] = deal(NaN(2640,length(scheds),2));

%res_mat dims: channel, schedule set, residual axis
for chind = 1:2640
    for schind = 1:length(scheds)
        ind = ismember(fd.schind,scheds{schind}) & fd.ch==chind;
        if ~isempty(find(ind))
            res_mat(chind,schind,1) = nanmin(fd.resx(ind));
            res_mat(chind,schind,2) = nanmin(fd.resy(ind));
            res_mat_rot(chind,schind,1) = nanmin(fd.resx_rot(ind));
            res_mat_rot(chind,schind,2) = nanmin(fd.resy_rot(ind));
        end
    end
end

%% Per-DK roll vs time
% Without and with moon correction

corrname = {'nocorr','withcorr'};
corrtitle = {'No Phase Correction','With Phase Correction'};
valname = {'tilt','roll'};
lims = {[44.84 44.92], [-0.13 -0.02]};
for valind = 1:2
    for corrind = 1:2
        if corrind == 1
            load('z:/dev/moon_analysis/perdk_mirror_parms_rtheta.mat')
        else
            load('z:/dev/moon_analysis/perdk_mirror_parms_phase_corrected_rtheta.mat')
            
        end
        
        val = mirrorparms(:,valind)';
        
        times = NaN(1,length(val));
        for schedind = 1:length(scheds)
            t = moonsch{scheds{schedind}(1)}.t1;
            times(schedind) = datenum(mjd2datestr(t),'yyyy-mmm-dd:HH:MM:SS');
            
        end
        
        
        fig = figure(7+corrind);
        fig.Position = [2000 200 1000 280];
        clf; hold on;
        scatter(times,val,14,dks,'filled')
        grid on
        datetick('x','mmm-dd')
        ylabel([valname{valind} ' Fit [Degrees]'])
        xlabel('Date [mmm-dd]')
        c = colorbar();
        c.Title.String = 'DK [Deg]';
        ylim(lims{valind})
        title({['Best-Fit Mirror ' valname{valind}  ' per-Obs'],corrtitle{corrind}})
        fname = fullfile(figdir,[valname{valind} '_vs_time_perdk_' corrname{corrind} '.png']);
        saveas(fig,fname)
        
    end
end

%% DK vs. DK scatter plots.
plttitles = {'X-residual','Y-Residual'};
resname = {'xres','yres'};
winscale = 1.2;
for valind = 1:2
    fig = figure(valind);
    fig.Position = [2000+500*(valind-1)*winscale 0*300 500*winscale 450*winscale];
    
    
    for schinda = 1:length(scheds)
        
        for schindb = 1:length(scheds)
            clf; hold on;
            scatter(res_mat_rot(:,schinda,valind),res_mat_rot(:,schindb,valind),10,'blue','filled','MarkerFaceAlpha',1);
            %scatter(res_mat(:,schinda,pltind),res_mat(:,schindb,pltind),10,p.r,'filled','MarkerFaceAlpha',1);
            plot([-1 1],[-1 1],'k--')
            
            grid on
            xlim([-1 1]*0.15)
            ylim([-1 1]*0.15)
            xlabel(['dk' titles{schinda} ' [Deg]'])
            ylabel(['dk' titles{schindb} ' [Deg]'])
            title({plttitles{valind},['dk' titles{schinda} 'vs dk' titles{schindb}]})
            figname = fullfile(figdir,['scatter_dk' titles{schinda} 'vs_dk' titles{schindb} '_' resname{valind} '.png']);
            saveas(fig,figname)
        end
    end
end


%% Per-dk variance:

unqdks = unique(dks);


[d, stdsx, stdsy] = deal([]);
for dkind = 1:length(unqdks)
    ind = find(dks==unqdks(dkind));
    
    if length(ind)>1
        for inda = 1:(length(ind)-1)
            for indb = 2:length(ind)
                if inda~=indb
                    res_mat_ax = res_mat_rot(:,ind(inda),1);
                    res_mat_bx = res_mat_rot(:,ind(indb),1);
                    stdsx(end+1) = nanstd(res_mat_ax-res_mat_bx)/sqrt(2);
                    res_mat_ay = res_mat_rot(:,ind(inda),2);
                    res_mat_by = res_mat_rot(:,ind(indb),2);
                    stdsy(end+1) = nanstd(res_mat_ay-res_mat_by)/sqrt(2);
                    d(end+1) = dks(ind(inda));
                end
            end
        end
    end
end


%% Quiver plots per "DK"

winscale = 1.5;
scaling = 15;
fig = figure(1);
fig.Position(3:4) = [500*winscale 450*winscale];
clf;

fittype = {'','_rtheta'};

corrname = {'nocorr','withcorr'};
corrtitle = {'No Correction','Phase Corrected'};

for fitind = 2%1:2
    for corrind = 2%1:2
        if corrind == 1
            load(sprintf('z:/dev/moon_analysis/perdk_mirror_parms%s.mat',fittype{fitind}))
            load('z:/dev/moon_analysis/moon_beam_fits_cut.mat')
            
            
        else
            load(sprintf('z:/dev/moon_analysis/perdk_mirror_parms_phase_corrected%s.mat',fittype{fitind}))
            load('z:/dev/moon_analysis/moon_beam_fits_phase_corrected_cut.mat')
            
        end
        
        [fd.resx_rot_perdk,...
            fd.resy_rot_perdk,...
            fd.resr_perdk] = deal(NaN(size(fd.ch)));
        for schind = 1:length(scheds)
            ind = ismember(fd.schind,scheds{schind});
            
            fdsch = structcut(fd,ind);
            mirror = struct();
            mirror.height = 1.3;
            mirror.tilt = mirrorparms(schind,1);
            mirror.roll = mirrorparms(schind,2);
            
            source = struct();
            source.azimuth = reshape(fdsch.az_cen_src,[],1);
            source.distance = 3.8e8*cosd(reshape(fdsch.el_cen_src,[],1));
            source.height = 3.8e8*sind(reshape(fdsch.el_cen_src,[],1));
            
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
            
            
            
            [x, y, phi] = beam_map_pointing_model(fdsch.az_cen,fdsch.el_cen,fdsch.dk_cen,model,'bicep3',mirror,source,[]);
            
            prxsch = prx(fdsch.ch);
            prysch = pry(fdsch.ch);
            resx = reshape(prxsch-x,size(fdsch.ch));
            resy = reshape(prysch-y,size(fdsch.ch));
            [resth, resr] = cart2pol(resx,resy);
            resth = resth - fdsch.dk_cen*pi/180;
            [resx_rot, resy_rot] = pol2cart(resth,resr);
            
            fd.resx_rot_perdk(ind) = resx_rot;
            fd.resy_rot_perdk(ind) = resy_rot;
            fd.resr_perdk(ind) = resr;
            
            %[prxsch, prysch] = pol2cart((p.theta(fdsch.ch)-fdsch.dk_cen')*pi/180,p.r(fdsch.ch));
            xtrack = [1, -1]*10;
            ytrack = [1, -1]*10;
            [x_track_mirr, y_track_mirr] = get_mirror_coords(fdsch.dk_cen,xtrack,ytrack,zeros(size(fdsch.ch)),mount,mirror);
            
            [x_mirr, y_mirr] = get_mirror_coords(fdsch.dk_cen,prxsch',prysch',zeros(size(fdsch.ch)),mount,mirror);
            [x_fit_mirr, y_fit_mirr] = get_mirror_coords(fdsch.dk_cen,x',y',zeros(size(fdsch.ch)),mount,mirror);
            resx_mirr = x_mirr-x_fit_mirr;
            resy_mirr = y_mirr-y_fit_mirr;
            
            clf;
            %quiver(prxsch',prysch',resx_rot*scaling,resy_rot*scaling,0)
            %quiver(prxsch',prysch',resx*scaling,resy*scaling,0)
            quiver(x_mirr,y_mirr, resx_mirr*scaling,resy_mirr*scaling,0)
            hold on;
            
            mk = {'^','+'};
            for j = 1:length(xtrack)
                plot(x_track_mirr(j),y_track_mirr(j),'k','MarkerSize',14,'Marker',mk{j})
            end
            %plot(x_mirr,y_mirr,'.')
            %hold on; plot(x_fit_mirr,y_fit_mirr,'.')
            grid on;
            %xlim([-1 1]*15)
            %ylim([-1 1]*15)
            xlim([-1 1]*0.5)
            ylim([-1 1]*0.8)
            xlabel('X [Deg]')
            ylabel('Y [Deg]')
            title({sprintf('Beam Center Residuals, x%i, %s', scaling,corrtitle{corrind}),...
                sprintf('Tilt: %1.2f  Roll: %1.3f  Date: %s',mirrorparms(schind,1),mirrorparms(schind,2),mjd2datestr(moonsch{scheds{schind}(1)}.t1))...
                })
            figname = fullfile(figdir,sprintf('quiver_dk%s_fitperdk_%s%s.png',titles{schind},corrname{corrind},fittype{fitind}));
            saveas(fig,figname)
        end
    end
end


%% Find the difference in focal plane orientation between uncorrected and corrected:


winscale = 1.5;
scaling = 15;
fig = figure(1);
fig.Position(3:4) = [500*winscale 450*winscale];
clf;

fittype = {'','_rtheta'};
corrname = {'nocorr','withcorr'};
corrtitle = {'No Correction','Phase Corrected'};
axisname = {'_x','_y'};
axistitle = {'X','Y'};

theta_means = NaN(2,length(scheds));
for fitind = 2%1:2
    for schind = 1:length(scheds)
        
        [theta_ch, times_ch] = deal(NaN(2,2640));
        
        for valind = 1:2
            if valind==1
                load(sprintf('z:/dev/moon_analysis/perdk_mirror_parms%s.mat',fittype{fitind}))
                load('z:/dev/moon_analysis/moon_beam_fits_cut.mat')
                
                
            else
                load(sprintf('z:/dev/moon_analysis/perdk_mirror_parms_phase_corrected%s.mat',fittype{fitind}))
                load('z:/dev/moon_analysis/moon_beam_fits_phase_corrected_cut.mat')
                
            end
            
            ind = ismember(fd.schind,scheds{schind});
            
            fdsch = structcut(fd,ind);
            mirror = struct();
            mirror.height = 1.3;
            mirror.tilt = mirrorparms(schind,1);
            mirror.roll = mirrorparms(schind,2);
            
            source = struct();
            source.azimuth = reshape(fdsch.az_cen_src,[],1);
            source.distance = 3.8e8*cosd(reshape(fdsch.el_cen_src,[],1));
            source.height = 3.8e8*sind(reshape(fdsch.el_cen_src,[],1));
            [x, y, chi] = beam_map_pointing_model(fdsch.az_cen,fdsch.el_cen,fdsch.dk_cen,model,'bicep3',mirror,source,[]);
            %[r, th, chi] = beam_map_pointing_model(fdsch.az_cen,fdsch.el_cen,fdsch.dk_cen,model,'bicep3',mirror,source,[],'rtheta');
            [th, r] = cart2pol(x,y);
            th = rad2deg(th);
            
            th_diff = wrapTo180(p.theta(fdsch.ch)-th);
            ind = p.r(fdsch.ch)>6 & abs(th_diff)<0.4;
            theta_means(valind,schind) = nanmedian(th_diff(ind));
            
            
            if 1
                prsch = p.r(fdsch.ch)';
                pthsch = p.theta(fdsch.ch)'-fdsch.dk_cen;
                [xsch ysch] = pol2cart(pthsch*pi/180,prsch);
                xval = {xsch, ysch};
                for axind = 1:2
                    
                    %scatter(p.r(fdsch.ch),th_diff,14,(fdsch.t_cen-nanmin(fdsch.t_cen))*24,'filled')
                    scatter(xval{axind},th_diff,14,(fdsch.t_cen-nanmin(fdsch.t_cen))*24,'filled')
                    grid on
                    c = colorbar();
                    %caxis([0,3])
                    %caxis([-10,10])
                    c.Title.String = 'Time [Hours]';
                    ylim([-1 1]*1)
                    xlim([-1 1]*15)
                    xlabel(sprintf('%s Derotated CMB-Observed [Degrees]',axistitle{axind}))
                    ylabel('\theta_{CMB} - \theta_{fit} [Degrees]')
                    title({sprintf('Beam Center Residuals, Derotated, %s', corrtitle{valind}),...
                        sprintf('Tilt: %1.2f  Roll: %1.3f  Date: %s',mirrorparms(schind,1),mirrorparms(schind,2),mjd2datestr(moonsch{scheds{schind}(1)}.t1))...
                        })
                    colormap parula
                    figname = fullfile(figdir,sprintf('thetares_vs_r_dk%s_fitperdk_%s%s%s.png',titles{schind},corrname{valind},fittype{fitind},axisname{axind}));
                    saveas(fig,figname)
                    
                end
            end
        end
        
    end
end

%% FP ort vs time.

winscale = 1;
scaling = 15;
fig = figure(3);
fig.Position(3:4) = [1000*winscale 280*winscale];
clf;

fittype = {'','_rtheta'};
corrname = {'nocorr','withcorr'};
corrtitle = {'No Correction','Phase Corrected'};

theta_means = NaN(2,length(scheds));
for fitind = 2%1:2
    
    for valind = 1:2
        if valind==1
            load(sprintf('z:/dev/moon_analysis/perdk_mirror_parms%s.mat',fittype{fitind}))
            load('z:/dev/moon_analysis/moon_beam_fits_cut.mat')
            
            
        else
            load(sprintf('z:/dev/moon_analysis/perdk_mirror_parms_phase_corrected%s.mat',fittype{fitind}))
            load('z:/dev/moon_analysis/moon_beam_fits_phase_corrected_cut.mat')
            
        end
        
        for schind = 1:length(scheds)
            
            
            [theta_ch, times_ch] = deal(NaN(2,2640));
            
            ind = ismember(fd.schind,scheds{schind});
            
            fdsch = structcut(fd,ind);
            mirror = struct();
            mirror.height = 1.3;
            mirror.tilt = mirrorparms(schind,1);
            mirror.roll = mirrorparms(schind,2);
            
            source = struct();
            source.azimuth = reshape(fdsch.az_cen_src,[],1);
            source.distance = 3.8e8*cosd(reshape(fdsch.el_cen_src,[],1));
            source.height = 3.8e8*sind(reshape(fdsch.el_cen_src,[],1));
            [x, y, chi] = beam_map_pointing_model(fdsch.az_cen,fdsch.el_cen,fdsch.dk_cen,model,'bicep3',mirror,source,[]);
            %[r, th, chi] = beam_map_pointing_model(fdsch.az_cen,fdsch.el_cen,fdsch.dk_cen,model,'bicep3',mirror,source,[],'rtheta');
            [th, r] = cart2pol(x,y);
            th = rad2deg(th);
            for chind = 1:size(theta_ch,2)
                ci = find(fdsch.ch==chind);
                if ~isempty(ci)
                    theta_ch(valind,chind) = nanmean(wrapTo180(th(ci)));
                    times_ch(valind,chind) = nanmean(fdsch.t_cen(ci));
                end
            end
            
            th_diff = wrapTo180(p.theta'-theta_ch(valind,:));
            ind = p.r'>6 & abs(th_diff)<0.4;
            theta_means(valind,schind) = nanmedian(th_diff(ind));
            
            
            
        end
        if 1
            scatter(times,theta_means(valind,:),14,dks,'filled')
            grid on
            c = colorbar();
            c.Title.String = 'DK [Deg]';
            ylim([-1 1]*0.2)
            xlabel('Time')
            ylabel('Mean \theta_{CMB}-\theta_{fit} [Degrees]')
            title({sprintf('Mean FP Orientation Angle , %s', corrtitle{valind}),...
                })
            datetick('x','mmm-dd','KeepTicks')
            figname = fullfile(figdir,sprintf('meantheta_vs_time_perdk_%s%s.png',corrname{valind},fittype{fitind}));
            saveas(fig,figname)
        end
    end
end



%% Test quivers

fig = figure(1);
clf

[x0,y0] = pol2cart(p.theta*pi/180,p.r);

scaling = 1;
rscale = 1.05;
throt = x0*5/15;
[x, y] = pol2cart((p.theta+throt)*pi/180,p.r.*rscale);

xdiff = x0-x;
ydiff = y0-y;

quiver(x,y,xdiff*scaling,ydiff*scaling,0)
grid on;

%% calculate the stat and sys uncertainties



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
            load('z:/dev/moon_analysis/perdk_mirror_parms.mat')
        case {2,3}
            load('z:/dev/moon_analysis/perdk_mirror_parms_phase_corrected.mat')
    end
    
    
    m = nanmean(mirrorparms(casesch{caseind},:),1);
    st = nanmean(mirrorerrs(casesch{caseind},:),1);
    
    
    tsch = timesch{caseind};
    for parmind = 1:2
        [sy,sy1, sy2] = deal([]);
        for timeind = 1:length(tsch)
            sy1(end+1) = nanstd(mirrorparms(tsch{timeind},parmind));
            sy2(end+1) = nanmean(mirrorparms(tsch{timeind},parmind));
        end
        sy(parmind) = sqrt(sum(sy1.^2)+((max(sy2)-min(sy2))/2).^2);
        txt{2+parmind,1} = [txt{2+parmind,1} sprintf(txtrslt,m(parmind),st(parmind),sy(parmind))];
    end
    
    mn(caseind,:) = m;
    sstat(caseind,:) = st;
    ssys(caseind,:) = sy;
    
    
end

disp(txt)





