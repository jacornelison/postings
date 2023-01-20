function grid_cal_posting_plots()
%% init variables
addpath(genpath('z:/pipeline/'))
addpath('z:/dev/rps')
figdir = fullfile('C:','Users','James','Documents','GitHub','postings','2022mmdd_rps_grid_cal','figs','');
datadir = fullfile('C:','Users','James','Documents','GitHub','postings','2022mmdd_rps_grid_cal','data','');

%% Tilt cals
%% Tilt Cals
clc
close all
figure(5)
hold on
level_cal = 85/3600.; % degrees per tick

rps_tilt_cals_all = {};
isaac_tilt_cals_all = {};

% First cal, indoors
rps_tilt_ticks = [0 1 2 1 0 0];
rps_tilt_meter = [0.113 0.210 0.302 0.204 0.111 0.110]; % in Volts
rps_tilt_cal = polyfit(rps_tilt_meter,level_cal*rps_tilt_ticks,1);
fprintf('RPS 03 Dec 2021: %f deg/V %+0.2fdeg\n', rps_tilt_cal(1),rps_tilt_cal(2))
plot(rps_tilt_meter,level_cal*rps_tilt_ticks,'.')


rps_tilt_ticks = [0 2 1 0 -1 -2 -1];
rps_tilt_meter = [0.1360 0.2610 0.1840 0.1160 0.0345 -0.0445 0.0390];
rps_tilt_cal = polyfit(rps_tilt_meter,level_cal*rps_tilt_ticks,1);
fprintf('RPS ?? Dec 2021: %f deg/V %+0.2fdeg\n', rps_tilt_cal(1),rps_tilt_cal(2))
plot(rps_tilt_meter,level_cal*rps_tilt_ticks,'.')

isaac_tilt_ticks = [0 2 0 1 -1 -2 0];
isaac_tilt_meter = [0.553 0.681 0.558 0.615 0.479 0.407 0.569]+0*(0.644-0.553); % in Volts
isaac_tilt_cal = polyfit(isaac_tilt_meter,level_cal*isaac_tilt_ticks,1);
fprintf('ISAAC 9 Dec 2021: %f deg/V %+0.2fdeg\n', isaac_tilt_cal(1),isaac_tilt_cal(2))
plot(isaac_tilt_meter,level_cal*isaac_tilt_ticks,'.')
isaac_tilt_ticks = [0 -1 -2 -2 -2 -1 0 1 2 1 0];
isaac_tilt_meter = [0.5585 0.4890 0.3885 0.3925 0.3975 0.4725 0.5485 0.6180 0.7180 0.6380 0.5485];
isaac_tilt_cal = polyfit(isaac_tilt_meter,level_cal*isaac_tilt_ticks,1);
fprintf('ISAAC 16 Dec 2021: %f deg/V %+0.2fdeg\n', isaac_tilt_cal(1),isaac_tilt_cal(2))
plot(isaac_tilt_meter,level_cal*isaac_tilt_ticks,'.')

rps_tilt_cals_all{end+1} = rps_tilt_cal;
isaac_tilt_cals_all{end+1} = isaac_tilt_cal;

% RPS - outdoor
inc_tilt_meter = [0.247,0.405,0.550,0.184];
inc_readout = [-0.04 -0.01 +0.03 -0.06]+0.055;
rps_tilt_cal = polyfit(inc_tilt_meter,inc_readout,1);
fprintf('RPS 24 Dec 2021: %f deg/V %+0.2fdeg\n', rps_tilt_cal(1),rps_tilt_cal(2))
plot(inc_tilt_meter,inc_readout,'.')

% isaac
inc_tilt_meter = [0.602,0.276,0.250,0.632,0.315,0.515];
inc_readout = -1*[-0.03 +0.06 +0.08 -0.03 +0.04 -0.03]+0.055;
isaac_tilt_cal = polyfit(inc_tilt_meter,inc_readout,1);
fprintf('ISAAC 24 Dec 2021: %f deg/V %+0.2fdeg\n', isaac_tilt_cal(1),isaac_tilt_cal(2))
plot(inc_tilt_meter,inc_readout,'.')
xlabel('Tilt meter [V]')
ylabel('Cal''d Tilt [^o]')
grid on
%legend()

rps_tilt_cals_all{end+1} = rps_tilt_cal;
isaac_tilt_cals_all{end+1} = isaac_tilt_cal;

clf
% Jun 25 calibrations
isaac_tilt_ticks = [0 -1 -2 -1 0 2 2];
isaac_tilt_meter = [0.4718 0.4054 0.314 0.389 0.458 0.545 0.5714];
isaac_inc_readout = [-0.05 -0.04 -0.09 -0.06 -0.03 0.01 0.01];
isaac_tilt_cal = polyfit(isaac_tilt_meter,level_cal*isaac_tilt_ticks,1);
fprintf('ISAAC 25 Jun 2022: %f deg/V %+0.2fdeg\n', isaac_tilt_cal(1),isaac_tilt_cal(2))

rps_tilt_ticks = [0 1 2 1 0 -1 -2 -2 -1 0 0 -2 -2 -1 -1 1 1 2 2];
rps_tilt_meter = [0.0278 0.119 0.2049 0.1208 0.051 -0.039 -0.118 -0.153 -0.035 0.068 0.026 -0.0987 -0.132 -0.051 -0.019 0.109 0.141 0.1875 0.216];
rps_inc_readout = [-0.02 -0.01 +0.02 -0.02 -0.05 -0.06 -0.09 -0.09 -0.07 -0.02 -0.02 -0.1 -0.1 -0.06 -0.07];
rps_tilt_cal = polyfit(rps_tilt_meter,level_cal*rps_tilt_ticks,1);
fprintf('RPS 25 Jun 2022: %f deg/V %+0.2fdeg\n', rps_tilt_cal(1),rps_tilt_cal(2))


plot(rps_tilt_meter,level_cal*rps_tilt_ticks,'.')

rps_tilt_cals_all{end+1} = rps_tilt_cal;
isaac_tilt_cals_all{end+1} = isaac_tilt_cal;

%% Bootstrap the cal data to get an uncertainty:

len = length(rps_tilt_ticks);
Niter = len-1;
parms = NaN(Niter,2);

for iterind = 1:Niter
    newidx = randi(len,1,len);
    parms(iterind,:) = polyfit(rps_tilt_meter(newidx),level_cal*rps_tilt_ticks(newidx),1);
end

parm_mean = mean(parms,1);
parm_std = std(parms,[],1);

%% Plot the Tilt Calibration
clc
fig = figure(31);
fig.Position(3:4) = [500 400];
clf; hold on;

lims = [-1 1]*0.25;
plot(rps_tilt_meter,level_cal*rps_tilt_ticks,'.')
plot(lims,polyval(rps_tilt_cals_all{end} ,lims),'k')
xlim(lims)

grid on
xlabel('Tilt Meter Output [V]')
ylabel('Starrett Level Angle [Degrees]')
text(-0.23,0.04,{'Tilt Calibration:',...
    sprintf('Slope [Deg/V]: %0.2f +/- %0.2f',rps_tilt_cals_all{end}(1),parm_std(1)),...
    sprintf('Offset [Deg]: %0.3f +/- %0.3f',rps_tilt_cals_all{end}(2),parm_std(2))})
ylim([-0.06 0.06])
title('Tilt Meter Calibration')
fname = 'tilt_cal.png';
saveas(fig,fullfile(figdir,fname))

%% Grab two days and plot them

% Load the 200MB tilt data file and downselect just two days for the posting.
% load('z:/dev/rps/rps_tilt_data_2022.mat')
% days = [59583,59587]; % 4th and 8th of Jan 2022
% ind = ismember(floor(lj_data.time),days);
% ljd = structcut(lj_data,ind);
% save(fullfile(datadir,'tilt_data_two_day.mat'),'ljd')
load(fullfile(datadir,'tilt_data_two_day.mat'))

days = [59583,59587];

clr = colormap('lines');
Tcal = 24;

% Account for the shift in scale factor and tilt zero due to temperature drift.
Tshift = (ljd.AIN0*100-Tcal)*0.0002; % 0.0002 deg/C
Tscale = (ljd.AIN0*100-Tcal)*0.0002; % 0.02%/C

fs = 2; % Hz sample rate
fc = 0.01; % Hz cutoff;
[b,a] = butter(6,fc/(fs/2));

ljd.AIN0 = filter(b,a,ljd.AIN0);
ljd.AIN2 = filter(b,a,ljd.AIN2);

% Account for the uncertainty of the calibration
tilt0 = polyval(parm_mean,ljd.AIN0)+Tshift;
tilt = tilt0;
tilt = [tilt polyval(parm_mean+parm_std.*[1,1],ljd.AIN0)+Tshift]; % slope+std off+std
tilt = [tilt polyval(parm_mean+parm_std.*[1,-1],ljd.AIN0)+Tshift]; % slope+std off-std
tilt = [tilt polyval(parm_mean+parm_std.*[-1,1],ljd.AIN0)+Tshift]; % slope-std off+std
tilt = [tilt polyval(parm_mean+parm_std.*[-1,-1],ljd.AIN0)+Tshift]; % slope-std off-std
maxtilt = max(tilt,[],2);
mintilt = min(tilt,[],2);
clear tilt;

Tday = (ljd.time-floor(ljd.time))*24;

clc
for dayind = 1:length(days)

    fig = figure(32);
    fig.Position(3:4) = [750 300];
    clf;
    subplot(2,1,1)
    hold on;
    ind = find(floor(ljd.time)==days(dayind));
    ind = ind(1000:end-1000);
    scatter(Tday(ind),maxtilt(ind),4,[1 1 1]*0.65,'filled')
    scatter(Tday(ind),mintilt(ind),4,[1,1,1]*0.65,'filled')
    scatter(Tday(ind),tilt0(ind),3,clr(1,:),'filled')

    grid on
    ylim([-1 1]*0.15);
    xlim([0,24])
    ylabel('Tilt [Deg]')

    subplot(2,1,2)
    hold on;
    scatter(Tday(ind),ljd.AIN2(ind)*100,14,'filled')
    grid on
    ylim([-25 0])
    xlim([0,24])
    ylabel('Tilt Temp. [^oC]')
    xlabel('Time of Day (Hours)')

    sgtitle(sprintf('Tilt output for: %s',mjd2datestr(days(dayind))))

    fname = sprintf('tilt_monitor_plot_day_%i.png',days(dayind));
    saveas(fig,fullfile(figdir,fname))
end


