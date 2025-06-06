function moon_analysis_plots_apr2022()
% For making moon analysis posting plots.
% This was run in Matlab 2019a and if only sparsely commented.
% Contact James Cornelison if you need help.

%Load the stuff.
load('z://dev/moon_analysis/moon_beam_fits.mat')
load('z://dev/moon_analysis/moonsch.mat')
load('z://dev/moon_analysis/fpu_data_obs.mat')
load('z://dev/moon_analysis/pm.mat')
% JPL Horizons data
fname = 'z:/dev/moon_analysis/horizons_result.txt';
f = fopen(fname);
horz_data = textscan(f,'%f%s%s%f%f%f%f%f%f','delimiter',',','HeaderLines',61);
fclose(f);


inrange = @(A,B,C) A<=C & A>=B;
figdir = 'C:\Users\James\Documents\GitHub\postings\20220428_mirror_analysis\figs\';

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


% Interp the horizons data to make it compatible with our stuff.
fd.az_cen_horz = interp1(horz_data{1}-2400000.5,horz_data{end-1},fd.t_cen);
fd.el_cen_horz = interp1(horz_data{1}-2400000.5,horz_data{end},fd.t_cen);
fd.resaz = wrapTo180(fd.az_cen_src-fd.az_cen_horz);
fd.resel = fd.el_cen_src-fd.el_cen_horz;

% Best fit Mirror params for Horizons moon.
mirror = struct();
mirror.height = 1.3;
mirror.tilt = 44.887;
mirror.roll = -0.229;

source = struct();
source.azimuth = reshape(fd.az_cen_horz,[],1);
source.distance = 3.8e8*cosd(reshape(fd.el_cen_horz,[],1));
source.height = 3.8e8*sind(reshape(fd.el_cen_horz,[],1));

[x_horz, y_horz, phi_horz] = beam_map_pointing_model(fd.az_cen,fd.el_cen,fd.dk_cen,model,'bicep3',mirror,source,[]);
fd.resx_horz = reshape(prx(fd.ch)-x_horz,size(fd.ch));
fd.resy_horz = reshape(pry(fd.ch)-y_horz,size(fd.ch));
[resth, resr] = cart2pol(fd.resx_horz,fd.resy_horz);
resth = resth - fd.dk_cen*pi/180;
[fd.resx_rot_horz, fd.resy_rot_horz] = pol2cart(resth,resr);
fd.resr_horz = resr;

%%%%%%%%
% Cuts!
%%%%%%%%
cutind = true(size(fd.ch));

% Timing: The as-fit moon position will be WAY off from the Horizons data
cutind = cutind & inrange(fd.resaz,nanmedian(fd.resaz)-1,nanmedian(fd.resaz)+1);

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
for pltind = 1:length(titles)
    
    fig = figure(pltind);
    fig.Position = [400 0*300 500 450]*1.2;
    clf;
    
    h = scatterhist(vals{pltind,1},vals{pltind,2},'NBins',[30,50],'Marker','.');
    grid on
    xlim([-1 1]*0.3)
    ylim([-1 1]*0.3)
    title({'Beam Center Residuals [CMB beams-fit]',titles{pltind}})
    xlabel('X [Deg]')
    ylabel('Y [Deg]')
    fname = fullfile(figdir,['scatterhist_' titles{pltind} '.png']);
    saveas(fig,fname)
end

clear vals


%% Show how derotation works

source = struct();
source.azimuth = reshape(fd.az_cen_src,[],1);
source.distance = 3.8e8*cosd(reshape(fd.el_cen_src,[],1));
source.height = 3.8e8*sind(reshape(fd.el_cen_src,[],1));

tilts = linspace(44.8,45,5);
rolls = linspace(-0.2,0.1,5);
titles = {'','Derotated'};

for tiltind = 1:length(tilts)
    for rollind = 1:length(rolls)
        mirror = struct();
        mirror.height = 1.3;
        mirror.tilt = tilts(tiltind);
        mirror.roll = rolls(rollind);
        
        [x, y, phi] = beam_map_pointing_model(fd.az_cen,fd.el_cen,fd.dk_cen,model,'bicep3',mirror,source,[]);
        resx = reshape(prx(fd.ch)-x,size(fd.ch));
        resy = reshape(pry(fd.ch)-y,size(fd.ch));
        [resth, resr] = cart2pol(resx,resy);
        resth = resth - fd.dk_cen*pi/180;
        [resx_rot, resy_rot] = pol2cart(resth,resr);
        
        vals = {resx,resy; resx_rot, resy_rot};
        
        for pltind = 2%1:2
            
            fig = figure(pltind);
            fig.Position = [400 0*300 500 450]*1.2;
            clf;
            
            scatter(vals{pltind,1},vals{pltind,2},10,fd.dk_cen,'filled')
            grid on
            xlim([-1 1]*0.3)
            ylim([-1 1]*0.3)
            title({'Beam Center Residuals [CMB beams-fit]',sprintf('Tilt: %2.1f  Roll: %1.2f',tilts(tiltind),rolls(rollind))})
            xlabel('X [Deg]')
            ylabel('Y [Deg]')
            c = colorbar();
            c.Title.String = 'DK [Deg]';
            mirrname = sprintf('T_%i_R_%i_',tiltind,rollind);
            colormap default
            fname = fullfile(figdir,['derot_demo_' mirrname titles{pltind} '.png']);
            saveas(fig,fname)
            
        end
    end
end

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
    for pltind = 1:length(titles)
        clf;
        val = cals(pltind,1)*fd.(fnames{pltind})+cals(pltind,2);
        q = quantile(val,[0.05 0.95]);
        [b ind] = sort(val);
        ind = ind(end:-1:1);
        pl = scatter(resval{rotind,1}(ind),resval{rotind,2}(ind),4,val(ind),'filled');
        axis equal
        c = colorbar();
        c.Title.String = units{pltind};
        xlabel('X-Residual [Deg]')
        ylabel('Y-residual [Deg]')
        title({'Beam Center Residuals [CMB beams-fit]', titles{pltind}})
        grid on
        xlim([-1 1]*0.3)
        ylim([-1 1]*0.3)
        colormap jet
        figname = fullfile(figdir,['scatter_' fnames{pltind} '_xyres' rotname{rotind} '.png']);
        saveas(fig,figname)
    end
    
    % Next in R/Theta
    for resind = 1:3
        for pltind = 1:length(titles)
            clf;
            val = cals(pltind,1)*fd.(fnames{pltind})+cals(pltind,2);
            plot(val,resval{rotind,resind},'.')
            ylabel([resname{resind} '-Residual [Deg]'])
            xlabel([titles{pltind} ' ' units{pltind}])
            title({'Beam Center Residuals [CMB beams-fit]', titles{pltind}})
            grid on
            %xlim([min(val) max(val)]*1.1)
            ylim([-1 1]*0.3)
            figname = fullfile(figdir,['scatter_' fnames{pltind} '_' resname{resind} 'res' rotname{rotind} '.png']);
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
for pltind = 1:2
    fig = figure(2+pltind);
    fig.Position = [2000+500*(pltind-1)*winscale 0*300 500*winscale 450*winscale];
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
                cha(chind) = nanmin(val{pltind}(inda));
                chb(chind) = nanmin(val{pltind}(indb));
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
    sgtitle({titles{pltind},'Schedules spaced 2 Months apart'})
    figname = fullfile(figdir,['scatter_' valname{pltind} '_long_tdiff.png']);
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
    for pltind = 1:2
        fig = figure(2+pltind);
        fig.Position = [2000+500*(pltind-1)*winscale 0*300 500*winscale 450*winscale];
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
                    cha(chind) = nanmin(val{pltind}(inda));
                    chb(chind) = nanmin(val{pltind}(indb));
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
        sgtitle({titles{pltind},'Schedules spaced 1 Week apart'})
        figname = fullfile(figdir,['scatter_' valname{pltind} '_mid' num2str(midind) '_tdiff.png']);
        saveas(fig,figname)
    end
    
    
end




%% Short time diffs

% Corresponding schedules
scha = [30 32 34 37 39 41 43 45 47 49 51];
schb = [31 33 35 38 40 42 44 46 48 50 52];


val = {fd.resx_rot,fd.resy_rot};
%val = {fd.resx_rot_horz,fd.resy_rot_horz};
titles = {'X-residual [Deg]','Y-Residual [Deg]'};
winscale = 1.2;
for pltind = 1:2
    fig = figure(pltind);
    fig.Position = [2000+500*(pltind-1)*winscale 0*300 500*winscale 450*winscale];
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
                cha(chind) = nanmin(val{pltind}(inda));
                chb(chind) = nanmin(val{pltind}(indb));
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
    sgtitle({titles{pltind},'Schedules spaced back-to-back'})
    figname = fullfile(figdir,['scatter_' valname{pltind} '_short_tdiff.png']);
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
%

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

%% DK vs. DK scatter plots.
plttitles = {'X-residual','Y-Residual'};
resname = {'xres','yres'};
winscale = 1.2;
for pltind = 1:2
    fig = figure(pltind);
    fig.Position = [2000+500*(pltind-1)*winscale 0*300 500*winscale 450*winscale];
    
    
    for schinda = 1:length(scheds)
        
        for schindb = 1:length(scheds)
            clf; hold on;
            scatter(res_mat_rot(:,schinda,pltind),res_mat_rot(:,schindb,pltind),10,'blue','filled','MarkerFaceAlpha',1);
            %scatter(res_mat(:,schinda,pltind),res_mat(:,schindb,pltind),10,p.r,'filled','MarkerFaceAlpha',1);
            plot([-1 1],[-1 1],'k--')
            
            grid on
            xlim([-1 1]*0.15)
            ylim([-1 1]*0.15)
            xlabel(['dk' titles{schinda} ' [Deg]'])
            ylabel(['dk' titles{schindb} ' [Deg]'])
            title({plttitles{pltind},['dk' titles{schinda} 'vs dk' titles{schindb}]})
            figname = fullfile(figdir,['scatter_dk' titles{schinda} 'vs_dk' titles{schindb} '_' resname{pltind} '.png']);
            saveas(fig,figname)
        end
    end
end

%% Quiver plots - median subtracted.

fig = figure(pltind);
fig.Position = [2000+500*(pltind-1)*winscale 0*300 500*winscale 450*winscale];
clf; hold on;
scaling = 10;
for schind = 1:length(scheds)
    clf;
    rx = (res_mat_rot(:,schind,1)-0*nanmedian(res_mat_rot(:,schind,1)))*scaling;
    ry = (res_mat_rot(:,schind,2)-0*nanmedian(res_mat_rot(:,schind,2)))*scaling;
    ind = abs(rx) < 0.1*scaling & abs(ry) < 0.1*scaling;
    quiver(prx(ind),pry(ind),rx(ind),ry(ind),0)
    xlim([-1 1]*15)
    ylim([-1 1]*15)
    xlabel('X [Deg]')
    ylabel('Y [Deg]')
    title(['DK ' titles{schind}])
    grid on
    figname = fullfile(figdir,['quiver_dk' titles{schind} '_mdesub.png']);
    saveas(fig,figname)
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

%% Look at fits as a function of time difference:

% 2-month time diffs
% Corresponding schedules
scha = [ 2  3  4  7  5  6 22 24 23 26 27 28];
schb = [53 54 55 56 57 58 59 60 61 62 63 64];
schb = [53 54 55 56 57 58 59 60 61 62 63 64];

%Back-to-back
% Corresponding schedules
scha = [30 32 34 37 39 41 43 45 47 49 51];
schb = [31 33 35 38 40 42 44 46 48 50 52];


% Intrument-fixed without mirror?
% source = struct();
% source.azimuth = reshape(fd.az_cen_src,[],1);
% source.distance = 3.8e8*cosd(reshape(fd.el_cen_src,[],1));
% source.height = 3.8e8*sind(reshape(fd.el_cen_src,[],1));
% source.azimuth = reshape(fd.az_cen_horz,[],1);
% source.distance = 3.8e8*cosd(reshape(fd.el_cen_horz,[],1));
% source.height = 3.8e8*sind(reshape(fd.el_cen_horz,[],1));
%
% [xish,yish,phiish] = beam_map_pointing_model(fd.az_cen,fd.el_cen,fd.dk_cen,model,'bicep3',[],source,[]);

val = {fd.resx_rot,fd.resy_rot};
val = {fd.resx_rot_horz,fd.resy_rot_horz};
%val = {xish, yish};
valname = {'xres','yres'};
titles = {'X-residual [Deg]','Y-Residual [Deg]'};
winscale = 1.2;
for pltind = 1:2
    fig = figure(6+pltind);
    fig.Position = [2000+500*(pltind-1)*winscale 0*300 500*winscale 450*winscale];
    clf; hold on;
    pl = [];
    count = 1;
    for schind = 1:length(scha)
        if ismember(schind,[1 4 7 10])
            subplot(2,2,count)
            count = count+1;
            hold on;
            [cha, chb] = deal(NaN(2640,1));
        end
        
        for chind = 1:2640
            inda = fd.schind==scha(schind) & fd.ch==chind;
            indb = fd.schind==schb(schind) & fd.ch==chind;
            if ~isempty(find(inda)) && ~isempty(find(indb))
                cha(chind) = nanmin(val{pltind}(inda));
                chb(chind) = nanmin(val{pltind}(indb));
                dk = ceil(nanmedian(fd.dk_cen(inda)));
            end
        end
        size(find(~isnan(cha)));
        %pl(schind) = scatter(cha,chb,10,'blue','filled','MarkerFaceAlpha',0.2);
        if ismember(schind,[3,6,9,11])
            edge = 0.1;
            edges = linspace(-1,1,100)*edge;
            N = histc(cha-chb,edges);
            bar(edges,N,'histc');
            %ylim([0,500])
            xlim([-1 1]*edge)
        end
        
        
        grid on
        %        xlim([-1 1]*0.15)
        %        ylim([-1 1]*0.15)
        xlabel('Schedule on 01 Dec 21')
        ylabel('Schedule on 25 Jan 22')
        title(sprintf('DK: %03i',-dk))
    end
    sgtitle({titles{pltind},'Schedules spaced 2 Months apart'})
    %figname = fullfile(figdir,['scatter_' valname{pltind} '_long_tdiff.png']);
    %saveas(fig,figname)
end
%legend(pl([1 4 7 10]),{'0', '90', '45', '135'})

%% Comparing a single timestream to horizons data

load('z:/dev/moon_analysis/timestream_moonsch_4_scan_1.mat')
% JPL Horizons data
fname = 'z:/dev/moon_analysis/horizons_moonsch_4_scan_1.txt';
f = fopen(fname);
horz_data_sch = textscan(f,'%f%s%s%f%f%f%f%f%f','delimiter',',','HeaderLines',61);
fclose(f);

fname = 'z:/dev/moon_analysis/horizons_moonsch_4_scan_1_geo.txt';
f = fopen(fname);
horz_data_sch_geo = textscan(f,'%f%s%s%f%f%f%f%s%s','delimiter',',','HeaderLines',61);
fclose(f);


%%
tarray = double(d.array.frame.utc);
tpmac = double(d.antenna0.time.utcslow);
ttrack = double(d.antenna0.tracker.utc);
daz = double(d.antenna0.tracker.horiz_topo(:,1))./3.6e6;
del = double(d.antenna0.tracker.horiz_topo(:,2))./3.6e6;
ts = (tarray-tarray(1))*24*3600;

winscale = 1;
fig = figure(1);
fig.Position = [2000+500*winscale 0*300 900*winscale 450*winscale];
clf; hold on;

subplot(2,3,1)
plot(ts,daz)
grid on
%xlabel('Time [sec]')
ylabel({'pmac.tracker.horiz\_topo[1]'})
title('GCP Moon Az [Deg]')

subplot(2,3,2)
plot(ts,del)
grid on
%xlabel('Time [sec]')
ylabel({'pmac.tracker.horiz\_topo[2]'})
title('GCP Moon El [Deg]')

subplot(2,3,3)
hold on;
plot(ts,(tarray-tpmac)*24*3600)
grid on
xlabel('Time [sec]')
ylabel({'frame.utc-time.utcslow [sec]'})
title('Array Time Minus PMAC Time')

subplot(2,3,4)
plot(ts(2:end),diff(daz))
grid on
xlabel('Time [sec]')
ylabel({'Diff Az [Deg]'})
ylim([-5 1]*1e-3)
%title('GCP Moon Az [Deg]')

subplot(2,3,5)
plot(ts(2:end),diff(del))
ylim([-14 2]*1e-4)
yticks([-14:2:2]*1e-4)
xticks([0:30:210])
grid on
xlabel('Time [sec]')
ylabel({'Diff El [Deg]'})
%title('GCP Moon El [Deg]')

fname = fullfile(figdir,'timestream_azeltimediff.png');
saveas(fig,fname)

%% Offset over the course of the scan
% -44.65E -88.99N 2.843
%2459549.532418982,*,m,   202.48816,   -4.87726,    202.76492,    -4.98750,  165.812519,     4.978869,
%2459549.532418982, , ,   202.48734,   -5.86587,    202.76488,    -5.97603,        n.a.,         n.a.,
a = 14.5539;
b = -0.0157;
refract = @(el) (a.*cosd(el).*sind(el).^3 + b.*sind(el).*cosd(el).^3) ./...
    (sind(el).^4 + (a.*sind(el).^2+3.0*b.*cosd(el).^2));


az_horz_arr = interp1(horz_data_sch{1}-2400000.5,horz_data_sch{end-1},tarray);
el_horz_arr = interp1(horz_data_sch{1}-2400000.5,horz_data_sch{end},tarray);
ra_horz_arr = interp1(horz_data_sch{1}-2400000.5,horz_data_sch{6},tarray);
dec_horz_arr = interp1(horz_data_sch{1}-2400000.5,horz_data_sch{7},tarray);
az_horz_trk = interp1(horz_data_sch{1}-2400000.5,horz_data_sch{end-1},ttrack);
el_horz_trk = interp1(horz_data_sch{1}-2400000.5,horz_data_sch{end},ttrack);
ra_horz_trk = interp1(horz_data_sch{1}-2400000.5,horz_data_sch{6},ttrack);
dec_horz_trk = interp1(horz_data_sch{1}-2400000.5,horz_data_sch{7},ttrack);
az_horz_pmc = interp1(horz_data_sch{1}-2400000.5,horz_data_sch{end-1},tpmac);
el_horz_pmc = interp1(horz_data_sch{1}-2400000.5,horz_data_sch{end},tpmac);
%az_horz_arr_no_ref = interp1(horz_data_sch_geo{1}-2400000.5,horz_data_sch_geo{end-1},tarray);
%el_horz_arr_no_ref = interp1(horz_data_sch_geo{1}-2400000.5,horz_data_sch_geo{end},tarray);
ra_horz_arr_no_ref = interp1(horz_data_sch_geo{1}-2400000.5,horz_data_sch_geo{6},tarray);
dec_horz_arr_no_ref = interp1(horz_data_sch_geo{1}-2400000.5,horz_data_sch_geo{7},tarray);
%az_horz_pmc_no_ref = interp1(horz_data_sch_geo{1}-2400000.5,horz_data_sch_geo{end-1},tpmac);
%el_horz_pmc_no_ref = interp1(horz_data_sch_geo{1}-2400000.5,horz_data_sch_geo{end},tpmac);

tarraystr = mjd2datestr(tarray);
tarraystr = datestr(datenum(tarraystr,'yyyy-mmm-dd:HH:MM:SS'),'yyyy/mm/dd HH:MM:SS');
[az_horz_arr_no_ref2 el_horz_arr_no_ref2] = RADec2AzEl_no_c(ra_horz_arr_no_ref, dec_horz_arr_no_ref,tarraystr,...
    -89.9911,-44.65,2.843);
% [az_horz_arr_no_ref2 el_horz_arr_no_ref2] = RADec2AzEl_no_c(dra, ddec,tarraystr,...
%     -89.9911,-44.65,2.843);

% [az_horz_arr_no_ref2 el_horz_arr_no_ref2] = radec2azel(dra, ddec,tarray,...
% -89.9911,-44.65);

dra = double(d.antenna0.tracker.equat_geoc(:,1))./3.6e6;
ddec = double(d.antenna0.tracker.equat_geoc(:,2))./3.6e6;

fig = figure(2);
%fig.Position = [2000+500*0*winscale -1*200 900*winscale 900*winscale];
clf; hold on;

subplot(2,2,1)
plot(ts,daz-az_horz_arr)
title('With Refraction')

subplot(2,2,2)
plot(ts,daz-az_horz_arr_no_ref2)
title('Without Refraction')

subplot(2,2,3)
plot(ts,del-el_horz_arr)

subplot(2,2,4)
plot(ts,del-el_horz_arr_no_ref2)



fig = figure(3);
fig.Position = [2000+500*1.7*winscale -1*200 900*winscale 900*winscale];
clf; hold on;

subplot(2,2,1)
plot(ts,dra-ra_horz_arr)
title('With Refraction')

subplot(2,2,2)
plot(ts,dra-ra_horz_arr_no_ref)
title('Without Refraction')

subplot(2,2,3)
plot(ts,ddec-dec_horz_arr)

subplot(2,2,4)
plot(ts,ddec-dec_horz_arr_no_ref)




%% Offset over the course of the whole season

tdatenum = datenum(mjd2datestr(fd.t_cen),'yyyy-mmm-dd:HH:MM:SS');
ind = true(size(fd.ch));
%ind = fd.schind<53;

fig = figure(3);
fig.Position = [2000+500*winscale 0*300 900*winscale 650*winscale];
clf; hold on;

subplot(3,1,1)
plot(tdatenum(ind),fd.resaz(ind),'.')
datetick('x','mmm-dd')
grid on
xlabel('Time [Mmm-dd]')
ylabel('Az Res  (GCP-Horzns) [Deg]')

subplot(3,1,2)
plot(tdatenum(ind),fd.resel(ind),'.')
datetick('x','mmm-dd')
grid on
xlabel('Time [Mmm-dd]')
ylabel('El Res  (GCP-Horzns) [Deg]')

subplot(3,1,3)
plot(fd.el_cen_src,fd.resel(ind),'.')
grid on
xlabel('GCP Moon Elevation [Deg]')
ylabel('El Res  (GCP-Horzns) [Deg]')

sgtitle('GCP Moon Position Compared to JPL Horizons Data')
fname = fullfile(figdir,'gcp_vs_horizons_all.png');
saveas(fig,fname)


%% Wind dir and Azimuth pointing

val = fd.resr_perdk;

fig = figure(5);
clf;
fig.Position = [2000 0*300 600*1.75 450*1.75];
[b ind] = sort(val);
%ind = ind(end:-1:1);
scaling = 10;
subplot(2,2,1)
angdiff = ((fd.winddir-47.35))*pi/180;
polarscatter(angdiff(ind),fd.windspeed(ind),scaling,val(ind),'filled')
ax = gca;
ax.RAxis.Label.String = 'Wind Sp. [m/s]';
ax.ThetaDir = 'clockwise';
ax.ThetaZeroLocation = 'top';
caxis([0,0.2])
c = colorbar();
c.Title.String = 'r-Res [Deg]';
%colormap jet
title('Wind Direction WRT Az=0')

subplot(2,2,2)
angdiff = ((fd.az_cen-model.az_zero))*pi/180;
polarscatter(angdiff(ind),fd.windspeed(ind),scaling,val(ind),'filled')
ax = gca;
ax.RAxis.Label.String = 'Wind Sp. [m/s]';
ax.ThetaDir = 'clockwise';
ax.ThetaZeroLocation = 'top';
caxis([0,0.2])
c = colorbar();
c.Title.String = 'r-Res [Deg]';
%colormap jet
title('Boresight Az')

subplot(2,2,3)
angdiff = ((fd.winddir-47.35)-(fd.az_cen-model.az_zero))*pi/180;
p1 = polarscatter(angdiff(ind),fd.windspeed(ind),scaling,val(ind),'filled');
ax = gca;
%ax.RAxis.Label.String = 'Wind Sp. [m/s]';
ax.ThetaDir = 'clockwise';
ax.ThetaZeroLocation = 'top';
caxis([0,0.2])
c = colorbar();
c.Title.String = 'r-Res [Deg]';
%colormap jet
title('Wind Dir - Boresight Az')

subplot(2,2,4)
angdiff = wrapTo180((fd.winddir-47.35)-(fd.az_cen-model.az_zero));
[b ind] = sort(fd.windspeed);
ind = ind(end:-1:1);
scatter(angdiff(ind),val(ind),10,fd.windspeed(ind),'filled')
c = colorbar();
c.Title.String = 'Wnd Sp. [Deg]';
%colormap jet
grid on
xlabel('Wind Dir - Boresight Az [Deg]')
title('Wind Dir - Boresight Az')
ylabel('r_{cmb} - r_{fit} [Deg]')
sgtitle('Beam Center Residuals [CMB beams-fit]')
ylim([0 1]*0.2)
saveas(fig,fullfile(figdir,'polarscatter_winddir_vs_az.png'))

%% Quiver plots per "DK"
% Copy this into Cannon Matlab.
% scheds = {...
%     [2 3 4],...
%     [5 6 7],...
%     [8 9 10],...
%     [11 12 13],...
%     [14 15 18],...
%     [16 17 19],...
%     [20 21],...
%     [22 23 24],...
%     [26 27 28],...
%     [30 32 34],...
%     [31 33 35],...
%     [37 39 41],...
%     [38 40 42],...
%     [43 45 47],...
%     [44 46 48],...
%     [49 51],...
%     [50 52],...
%     [53 54 55],...
%     [56 57 58],...
%     [59 60 61],...
%     [62 63 64],...
%     [65 66 67],...
%     [68 69 70],...
%     [71 72 73],...
%     [74 75 76]};
% 
% parms = NaN(length(scheds),2);
% 
% for schind = 1:length(scheds)
% ind = ismember(fd.schind,scheds{schind});
% fdsch = structcut(fd,ind);
% 
% fdsch = moon_fit_mirror(fdsch);
% 
% parms(schind,:) = fdsch.fitparam;
% 
% end
% 


load('z:/dev/moon_analysis/perdk_mirror_parms.mat')

winscale = 1.5;
scaling = 15;
fig = figure(1);
fig.Position = [2000+500*0*winscale 0*300 500*winscale 450*winscale];
clf; hold on;

clear fdsch;
[fd.resx_rot_perdk,...
fd.resy_rot_perdk,...
fd.resr_perdk] = deal(NaN(size(fd.ch)));
for schind = 1:length(scheds)
    ind = ismember(fd.schind,scheds{schind});
    
    fdsch = structcut(fd,ind);
    mirror = struct();
    mirror.height = 1.3;
    mirror.tilt = parms(schind,1);
    mirror.roll = parms(schind,2);
    
    source = struct();
    source.azimuth = reshape(fdsch.az_cen_src,[],1);
    source.distance = 3.8e8*cosd(reshape(fdsch.el_cen_src,[],1));
    source.height = 3.8e8*sind(reshape(fdsch.el_cen_src,[],1));
    
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
    
    clf;
    quiver(prxsch',prysch',resx_rot*scaling,resy_rot*scaling,0)
    grid on;
    xlim([-1 1]*15)
    ylim([-1 1]*15)
    xlabel('X [Deg]')
    ylabel('Y [Deg]')
    title({sprintf('Beam Center Residuals, Derotated, x%i',scaling),sprintf('Tilt: %1.2f  Roll: %1.2f',parms(schind,1),parms(schind,2))})
    figname = fullfile(figdir,['quiver_dk' titles{schind} '_fitperdk.png']);
    saveas(fig,figname)
end

%%
ind = true(size(fd.ch));
ind = inrange(((fd.el_cen(ind)-model.el_zero)+fd.el_cen_src(ind))-90,-1,1) &...
      fd.resy_rot>0.1;
ind = find(ind);
[b i] = min(fd.resy_rot(ind))

fd.ch(ind(i))
moonsch{fd.schind(ind(i))}.name
fd.scanind(ind(i))

%%
clf


hist(fd.resy_rot_perdk(fd.schind==35&fd.scanind==10),20)


