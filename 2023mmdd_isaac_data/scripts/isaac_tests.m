function isaac_tests()
%%
set(groot,'defaulttextinterpreter','latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(groot,'defaultAxesFontSize',12)

%% Tilt Cals

clc
%close all
%figure(5)
%clf; hold on;
level_cal = 85/3600.; % degrees per tick

rps_tilt_cals_all = {};
isaac_tilt_cals_all = {};

% First cal, indoors
rps_tilt_ticks = [0 1 2 1 0 0];
rps_tilt_meter = [0.113 0.210 0.302 0.204 0.111 0.110]; % in Volts
rps_tilt_cal = polyfit(rps_tilt_meter,level_cal*rps_tilt_ticks,1);
fprintf('RPS 03 Dec 2021: %f deg/V %+0.2fdeg\n', rps_tilt_cal(1),rps_tilt_cal(2))
%plot(rps_tilt_meter,level_cal*rps_tilt_ticks,'b.')


rps_tilt_ticks = [0 2 1 0 -1 -2 -1];
rps_tilt_meter = [0.1360 0.2610 0.1840 0.1160 0.0345 -0.0445 0.0390];
rps_tilt_cal = polyfit(rps_tilt_meter,level_cal*rps_tilt_ticks,1);
fprintf('RPS 16 Dec 2021: %f deg/V %+0.2fdeg\n', rps_tilt_cal(1),rps_tilt_cal(2))
%plot(rps_tilt_meter,level_cal*rps_tilt_ticks,'bx')

isaac_tilt_ticks = [0 2 0 1 -1 -2 0];
isaac_tilt_meter = [0.553 0.681 0.558 0.615 0.479 0.407 0.569]+0*(0.644-0.553); % in Volts
isaac_tilt_cal = polyfit(isaac_tilt_meter,level_cal*isaac_tilt_ticks,1);
fprintf('ISAAC 9 Dec 2021: %f deg/V %+0.2fdeg\n', isaac_tilt_cal(1),isaac_tilt_cal(2))
%plot(isaac_tilt_meter,level_cal*isaac_tilt_ticks,'r.')
isaac_tilt_ticks = [0 -1 -2 -2 -2 -1 0 1 2 1 0];
isaac_tilt_meter = [0.5585 0.4890 0.3885 0.3925 0.3975 0.4725 0.5485 0.6180 0.7180 0.6380 0.5485];
isaac_tilt_cal = polyfit(isaac_tilt_meter,level_cal*isaac_tilt_ticks,1);
fprintf('ISAAC 16 Dec 2021: %f deg/V %+0.2fdeg\n', isaac_tilt_cal(1),isaac_tilt_cal(2))
%plot(isaac_tilt_meter,level_cal*isaac_tilt_ticks,'rx')

rps_tilt_cals_all{end+1} = rps_tilt_cal;
isaac_tilt_cals_all{end+1} = isaac_tilt_cal;

% RPS - outdoor
inc_tilt_meter = [0.247,0.405,0.550,0.184];
inc_readout = [-0.04 -0.01 +0.03 -0.06]+0.055;
rps_tilt_cal = polyfit(inc_tilt_meter,inc_readout,1);
fprintf('RPS 24 Dec 2021: %f deg/V %+0.2fdeg\n', rps_tilt_cal(1),rps_tilt_cal(2))
%plot(inc_tilt_meter,inc_readout,'b^')

% isaac
inc_tilt_meter = [0.602,0.276,0.250,0.632,0.315,0.515];
inc_readout = -1*[-0.03 +0.06 +0.08 -0.03 +0.04 -0.03]+0.055;
isaac_tilt_cal = polyfit(inc_tilt_meter,inc_readout,1);
fprintf('ISAAC 24 Dec 2021: %f deg/V %+0.2fdeg\n', isaac_tilt_cal(1),isaac_tilt_cal(2))
%plot(inc_tilt_meter,inc_readout,'r^')
%xlabel('Tilt meter [V]')
%ylabel('Cal''d Tilt [Deg]')
%grid on
%legend()

rps_tilt_cals_all{end+1} = rps_tilt_cal;
isaac_tilt_cals_all{end+1} = isaac_tilt_cal;

% Jun 21 RPS Cal?
rps_tilt_ticks = [0 -1 -2 -1 0 1 2 1 0];
rps_tilt_meter = [0.0413 -0.0368 -0.1160 -0.0368 0.0467 0.1074 0.1929 0.1166 0.04684];
rps_tilt_cal = polyfit(rps_tilt_meter,level_cal*rps_tilt_ticks,1);
%fprintf('RPS 21 Jun 2022: %f deg/V %+0.2fdeg\n', rps_tilt_cal(1),rps_tilt_cal(2))
%plot(rps_tilt_meter,level_cal*rps_tilt_ticks,'bo')
%grid on
%rps_tilt_cals_all{end+1} = rps_tilt_cal;

%clf; hold on;
% Jun 25 calibrations
isaac_tilt_ticks = [0 -1 -2 -1 0 2 2];
isaac_tilt_meter = [0.4718 0.4054 0.314 0.389 0.458 0.545 0.5714];
isaac_inc_readout = [-0.05 -0.04 -0.09 -0.06 -0.03 0.01 0.01];
isaac_tilt_cal = polyfit(isaac_tilt_meter,level_cal*isaac_tilt_ticks,1);
fprintf('ISAAC 25 Jun 2022: %f deg/V %+0.2fdeg\n', isaac_tilt_cal(1),isaac_tilt_cal(2))
%plot(isaac_tilt_meter,level_cal*isaac_tilt_ticks,'r+')

rps_tilt_ticks = [0 1 2 1 0 -1 -2 -2 -1 0 0 -2 -2 -1 -1 1 1 2 2];
rps_tilt_meter = [0.0278 0.119 0.2049 0.1208 0.051 -0.039 -0.118 -0.153 -0.035 0.068 0.026 -0.0987 -0.132 -0.051 -0.019 0.109 0.141 0.1875 0.216];
rps_inc_readout = [-0.02 -0.01 +0.02 -0.02 -0.05 -0.06 -0.09 -0.09 -0.07 -0.02 -0.02 -0.1 -0.1 -0.06 -0.07];
rps_tilt_cal = polyfit(rps_tilt_meter,level_cal*rps_tilt_ticks,1);
fprintf('RPS 25 Jun 2022: %f deg/V %+0.2fdeg\n', rps_tilt_cal(1),rps_tilt_cal(2))
%plot(rps_tilt_meter,level_cal*rps_tilt_ticks,'b+')
%grid on
rps_tilt_cals_all{end+1} = rps_tilt_cal;
isaac_tilt_cals_all{end+1} = isaac_tilt_cal;

% Load and fit the data.
%clc
if ~ispc
addpath('~/Documents/kovac_lab/rps_benchtests/rsynced_from_cannon/beammap/')
addpath('~/Documents/kovac_lab/rps_benchtests/rsynced_from_cannon/util/')
addpath('~/Documents/kovac_lab/rps_benchtests/postings/2023mmdd_isaac_data/scripts/')
gitdir = fullfile('/','home','annie','Documents','kovac_lab','rps_benchtests');
%direct = fullfile(gitdir,'rps_cal','modified_data', 'even_angles');
%direct = fullfile(gitdir,'rps_cal','modified_data', 'odd_angles');
figdir = fullfile(gitdir,'postings', '2023mmdd_isaac_data', 'figs');
else
    addpath('z:/pipeline/util/')
    addpath('z:/pipeline/beammap/')
    gitdir = fullfile('c:/','Users','James','Documents','GitHub');
    figdir = fullfile(gitdir,'postings','2023mmdd_isaac_data','figs');
end

direct = fullfile(gitdir,'rps_cal','data');
[prefix, labs] = get_pref_and_labs(); % filenames and labels are moved to the bottom now...

if ~exist('meanguess','var')
    meanguess = 0;
end

% Create the struct that all the fit data (fd) will go into.
fd = struct();
flds = {'schnum','rownum','phi','eps','N1','N2','A','time','tilt','istilt','temp','istemp','istilttemp','wgt','fitind','phi_isaac','obsnum'};
vals = {'prefind','di','parm(1)','parm(2)','parm(3)','parm(4)','parm(5)','nanmean(TIME)',...
    'nanmean(T)','nanmean(T2)','nanmean(TEMP)','nanmean(ISTEMP)','nanmean(ISTILTTEMP)','wgtind','fitind','nanmean(a_isaac)','obsnum'};
for fldind = 1:length(flds)
    fd.(flds{fldind}) = [];
end


[reses, angs, modcurves, modstds, ts,tilts,istilts] = deal({});
rsmax = [];
thresh = 0.7;
weighttitles = {'Uniform',sprintf('Amp<%0.1f',thresh),'1/Amp'};
weightnames = {'uniform','threshold','one_on_r'};
cmlines = colormap('lines');

fitnames = {'fmin','lsq','complex'};
plotmodcurve = 0;
obsnum = 1;
dists = [ones(1,27), 20*0.0254, ones(1,3)*40*0.0254,ones(1,5)*20*0.0254,ones(1,7)*37*0.0254, ones(1,11)*26*0.0254];

% Loop over the file names and fit.
for prefind = 40:length(prefix)%[6 7 9:15, 20:42 44]%20:42

    % Tilt Calibrations
    if ismember(prefind,[1:18])
        rpscal = rps_tilt_cals_all{end-2};
        isaaccal = isaac_tilt_cals_all{end-2};
    elseif ismember(prefind,[19:21 24:26])
        rpscal = rps_tilt_cals_all{end-2};
        isaaccal = isaac_tilt_cals_all{end-2};
    elseif ismember(prefind,[22:23])
        rpscal = rps_tilt_cals_all{end-2};
        isaaccal = isaac_tilt_cals_all{end-2};
    elseif ismember(prefind,[27:length(prefix)])
        rpscal = rps_tilt_cals_all{end};
        isaaccal = isaac_tilt_cals_all{end};
    end
    isaaccal(1) = 0.3;
    %rpscal(1) = 0.3;
    
    d = dir(fullfile(direct,[prefix{prefind} '_scan*']));
    data = {};
    for di = 1:length(d)
        if ~ismember(prefind, [8:21 22:26])
            data{di} = load_lj_data(fullfile(direct,d(di).name));
        else
            data{di} = load_lj_data(fullfile(direct,[prefix{prefind} ...
                sprintf('_scan_%i_loop_1.csv',di-1)]));
        end
    end

    fittype = 'lsq';
    fitind = find(ismember(fitnames,'lsq'));
    % Change the weighting if we want.
    for wgtind = 1%1:3

        % Complex fit requires more parameters (probably don't use)
        if strcmp(fittype,'complex')
            parms = NaN(length(data),7);
        else
            parms = NaN(length(data),5);
        end
        [times, chis, itempstd, itempmn] = deal(NaN(length(data)));
        [mod_curve, res] = deal(NaN(length(unique(data{1}.Angle)),length(data)));

        % Load the data
        for di = 1:length(data)
            a = unique(data{di}.Angle);
            [R, Rs, T, T2, TIME, TEMP, ISTILTTEMP, ISTEMP] = deal([]);
            for ai = 1:length(a)
                aind = data{di}.Angle == a(ai);
                R = [R; mean((sqrt(data{di}.X(aind).^2+data{di}.Y(aind).^2)))]; % could replace with mean(X) if we know we're peaked up 
                Rs = [Rs; std((sqrt(data{di}.X(aind).^2+data{di}.Y(aind).^2)))];
                T = [T; median(polyval(rpscal,data{di}.Tilt(aind)))];
                TIME = [TIME; mean(data{di}.Time(aind))];
                TEMP = [TEMP; mean(data{di}.Temp(aind))];
                if prefind>7
                    ISTEMP = [ISTEMP; data{di}.ISTemp(aind)];
                    ISTILTTEMP = [ISTILTTEMP; data{di}.ISTiltTemp(aind)];
                end

                try
                    T2 = [T2; median(polyval(isaaccal,data{di}.ISTilt(aind)))];
                catch
                    T2 = 0;
                end
                if ai==1 & di==1
                    ts{end+1} = sqrt(data{di}.X(aind).^2+data{di}.Y(aind).^2);
                end
            end

            %R = R./nanmax(R);
            %Rs = Rs./R;
            rsmax(end+1) = max(Rs);
            %
            switch wgtind
                case 1
                    w = 1;
                case 2
                    w = R<0.7*max(R);
                case 3
                    w = 1./(R+1)-min(1./(R+1));
            end

            itempstd(di) = nanstd(ISTILTTEMP);
            itempmn(di) = nanmean(ISTILTTEMP);
            if ismember(prefind, [6:18])
                homeangle = 0.77;
            elseif ismember(prefind,44)
                % Wasn't paying attention to homing position when I put
                % the RPS back together in Jul 2023
                homeangle = 5.454545*2;
            elseif ismember(prefind,[20,22:length(prefix)])
                homeangle = 5.454545;
            else
                homeangle = 0;
            end

            % We're measuring the angle of the ISAAC WRT gravity:
            % Clockwise looking at ISAAC is positive
            % RPS has opposite parity
            rps_tilt = T;
            rps_grid = homeangle;
            a = a + rps_tilt - rps_grid;

            % ISAAC grid and tiltmeter are set up the same as the RPS
            % So the ISAAC angle in the RPS coordinate system is
            % flipped sign.
            isaac_tilt = T2;
            isaac_grid = 0.19;
            a_isaac = -1*mean(isaac_tilt - isaac_grid);
            mod_curve(:,di) = R;

            mxfev = 100000;
            mxiter = 100000;
            options = optimset('TolFun',1e-10,'MaxIter',mxiter,'MaxFunEvals',mxfev,'Display','off');
            %lb = [-20 -0.5 -10 -10 0]; % angle, xpol, collimation x2, amplitude
            %ub = [20 0.5 10 10 1e6];

            lb = [-180 -0.5 -10 -10 0]; % angle, xpol, collimation x2, amplitude
            ub = [180 0.5 10 10 1e6];

            guess = [1e-6,1e-6,0,0,nanmax(R)/2]; % old fminsearch didn't leave zero if given zero as a guess
            %fittype = 'lsq';

            % A million ways to fit:
            switch fittype
                case 'basic' % Don't fit collimation params
                    lb = lb([1 2 end]);
                    ub = ub([1 2 end]);
                    modfunc = @(p,ang) p(3)*(cosd(2*(ang-p(1)))-(p(2)+1)/(p(2)-1));

                    parm = lsqcurvefit(modfunc,[0,0,1]+1e-6,a,R,lb,ub,options);
                    res(:,di) = R-modfunc(parm,a);
                    parm = [parm(1) parm(2) 0 0 parm(3)];
                case 'fmin' % Fminsearch is garbage.
                    modfunc = @(ang,e,A,N1,N2) A/2*(cosd(2*(a-ang))-(e+1)/(e-1)).*(N1*cosd(a)+N2*sind(a)+1);
                    chifun = @(x) sum((R-modfunc(x(1),x(2),x(3),x(4),x(5))).^2);

                    chifun = @(x) sum(w.*(R-rps_get_mod_model(x,a)).^2);

                    parm = fminsearch(chifun,guess,options);
                    res(:,di) = R-modfunc(parm(1),parm(2),parm(3),parm(4),parm(5));

                case 'lsq0' % LSQ works better for some reason. And doesn't require compiling like matmin.
                    parm = lsqcurvefit(@rps_get_mod_model,guess,a,R,lb,ub,options);
                    res(:,di) = R-rps_get_mod_model(parm,a);

                case 'lsq' % LSQ but choose your own chi-sq. Allows weighting.
                    chifunc = @(p) (w.*(R-rps_get_mod_model(p,a)));
                    parm = lsqnonlin(chifunc,guess,lb,ub,options);
                    res(:,di) = R-rps_get_mod_model(parm,a);
                case 'complex' % A mod curve model with complicated geometry considerations.
                    MIRROR = false;
                    PLOT = false;
                    source = struct();
                    source.azimuth = 0;
                    source.elevation = 0;
                    source.distance = dists(prefind); % meters
                    source.height = source.distance*tand(source.elevation); % meters;
                    source.type = 'other';
                    B3 = [1,0,0];
                    B1 = [0,0,1];
                    A = [0,0,0];

                    if 0
                        lb = [-10 -0.5 -10 -10 -10 -10 0 -0.5 -0.5];
                        ub = [10 0.5 10 10 10 10 10 0.5 0.5];
                        guess = [0,0,0,0,0,0,1,0,0];

                        modfunc = @(p) rps_get_mod_model_vectors(p(1:7),a,[0,p(8),p(9)],B3,B1,source,false,false);
                        chifunc = @(p) w.*(R-modfunc(p)');
                    else
                        lb = [-10 -0.5 -10 -10 -10 -10 0];
                        ub = [10 0.5 10 10 10 10 10];
                        guess = [0,0,0,0,0,0,1];
                        chifunc = @(p) ( 1.*(R-rps_get_mod_model_vectors(p,a',A,B3,B1,source,MIRROR,PLOT)') );
                    end
                    try

                        parm = lsqnonlin(chifunc,guess+1E-6,lb,ub,options);
                    catch
                        parm = NaN(1,7);
                    end
                    res(:,di) = chifunc(parm);
                    parm(1) = -parm(1);
            end

            parms(di,1:length(parm)) = parm;
            times(di) = mean(TIME);

            % Do a shitty cut.... do we still need this?
            if abs(R(1)-R((length(R)+1)/2))<0.2
                for fldind = 1:length(flds)
                    eval(sprintf('fd.%s(end+1) = %s;',flds{fldind},vals{fldind}))
                end
                angs{end+1} = a;
                reses{end+1} = R-rps_get_mod_model(parm,a);
                modcurves{end+1} = R;
                modstds{end+1} = Rs;
                tilts{end+1} = T;
                istilts{end+1} = T2;
            else
                parms(di,1:length(parm)) = NaN(1,length(parm));
            end

        end

    end

    fprintf('%i: %s\nAngle: %0.3f +/- %0.3f\n',prefind,labs{prefind},nanmean(parms(:,1))-a_isaac,nanstd(parms(:,1)))
    %fprintf('Eff: %0.4f +/- %0.4f\n',nanmean(parms(:,2)),nanstd(parms(:,2)))
    obsnum = obsnum+1;
end

% Analysis done.
% Plot all the things!
%% Plot a history of measurements
% Only uses measurements 6 to 42
% specifically: [6 7 9:15, 20:42 44]
cmlines = colormap('lines');

fig = figure(2+wgtind);
fig.Position(3:4) = [1900 900]*1;
clf; hold on;
ind = true(size(fd.schnum));

t = tiledlayout(1,1);
t.Padding = 'none';
t.TileSpacing = 'compact';
nexttile
hold on
clear s;
scheds = unique(fd.obsnum);%[12:15, 20:42]; 
[phi, pstd] = deal(NaN(1,length(scheds))); 
for schind = 1:length(scheds)
    idx = fd.obsnum==scheds(schind);
    phi(schind) = mean(fd.phi(idx)-fd.phi_isaac(idx));
    pstd(schind) = std(fd.phi(idx)-fd.phi_isaac(idx));
end
s(2) = plot(fd.obsnum(ind),-1*fd.istilt(ind),'x','Color',cmlines(3,:),'MarkerSize',16,'LineWidth',1.5);
s(3) = plot(fd.obsnum(ind),fd.tilt(ind),'.','MarkerSize',16,'Color',cmlines(2,:));
%s(1) = scatter(fd.schnum(ind),fd.phi(ind)-fd.phi_isaac(ind),14,cmlines(1,:));%,'filled');
s(1) = errorbar(scheds,phi,pstd,'.','Color',cmlines(1,:),'LineWidth',1,'MarkerSize',16);

%plot(fd.schnum(ind),fd.phi_isaac(ind),'^','MarkerSize',10)
plot([1 1]*4.5+5,[-1 0.5],'k--','LineWidth',2)
plot([1 1]*11.5+5,[-1 1],'k--','LineWidth',2)
ylim([-0.6 1])
grid on
xlim([min(scheds)-1,max(scheds)+1])

% Show where we homed during the harvard measurements
plot([11.8 18.2]+5, [1 1]*-0.5,':','LineWidth',2,'Color',[1 1 1]*0.5)
text(12.75+5,-0.45,'Homed once at beginning of dataset','FontSize',13)
plot([18.8 27.2]+5, [1 1]*-0.5,':','LineWidth',2,'Color',[1 1 1]*0.5)
text(20.95+5,-0.45,'Homed before every measurement','FontSize',13)


ax = gca();
legend(s,{'\phi_{fit} minus \phi_{actual}','ISAAC Tilt Angle','RPS Tilt Angle'},...
    'location','northwest','FontSize',13)
text(4.5+2,0.65,'Pole Measurements','FontSize',16)
text(1.25+2,0.45,'No Fine Homing Sw.','FontSize',13)
text(6.5+5,0.45,'Added Fine Homing Sw.','FontSize',13)
text(16.5+5,0.85,'2022 Measurements at Harvard','FontSize',16)
ylabel('Angle [Degrees]')
xlabs = {...
    'Indoor, Homed Each',...    
    'Indoor, Homed Once',...
    'Outdoor,Homed Once',...
    'Outdoor,Homed Once',...
    'Outdoor,Homed Once',...
    'Outdoor,Homed Once',...
    'Outdoor,Homed Each',...
    'Outdoor,Homed Each',...
    'Indoor, Homed Each',...
    'Indoor, Pre-Obs1',...
    'Indoor, Pre-Obs2',...
    'Outdoor, Pre-Obs',...
    'Outdoor, Post-Obs',...
    'Indoor, Post-Obs1',...
    'Indoor, Post-Obs2',...
    'Indoor, Post-Obs3',...
    'Dist 0.5m' ...
    'Dist 1m',...
    'Dist 1m',...
    'Dist 1m, new tilt',...
    'Dist 1m, old tilt',...
    'Dist 0.5m',...
    'Dist 0.5m',...
    'Dist 0.5m',...
    'Dist 0.5m',...
    'Dist 0.5m',...
    'D 0.5m, Y-off 2cm',...
    'D 1m, Y-off 2cm',...
    'D 1m, Y-off 2cm',...
    'D 1m, Y-off 2cm',...
    'D 1m, Y-off 2cm',...
    'D 1m, Y-off 2cm, Baffled',...
    };
ax.XTickLabelRotation = -60;
ax.XTick = linspace(min(scheds)-1,max(scheds)+1,length(scheds)+2);
ax.XTickLabel(2:end-1) = xlabs;
ax.XTickLabel([1, end]) = {'',''};

saveas(fig,'C:\Users\James\Documents\GitHub\postings\2023mmdd_isaac_data\figs\angle_fit_vs_actual_plot','png')

%% Look the other jig data for diagnostics

vals = {fd.temp*1000, fd.istemp*1000, fd.istilttemp*100+273.15};
%vals = {fd.N1,fd.N2}
%vals = {fd.A};
fig = figure(10);
clf; hold on;
for valind = 1:length(vals)
    ind = true(size(fd.schnum)) & fd.wgt == wgtind;
    scatter(fd.schnum(ind),vals{valind}(ind),'.');
    %plot(fd.phi(ind),vals{valind}(ind),'.');
    % ind = ismember(fd.schnum,[28:31 37:39]) & fd.wgt == wgtind;
    % plot(fd.schnum(ind),vals{valind}(ind),'.');
    % ind = ~ismember(fd.schnum,[28:31 37:39]) & fd.wgt == wgtind;
    % plot(fd.schnum(ind),vals{valind}(ind),'.');
end
%ylim([-1 7])
grid on

%% Look at time between measurements
% You can tell when we're homing and when we're not.
dt = NaN(1,max(fd.obsnum));
for obsind = 5:max(fd.obsnum)
    idx = find(fd.obsnum==obsind);
    [s si] = sort(fd.time(idx));
    dt(obsind) = median(diff(fd.time(idx(si))));
end

fig = figure(91824);
plot(dt)

%% Sigma angle vs Attenuation

att = [0,30,40,55];
fit_a = [-2.245, -2.238,-2.170,-1.595]+2.245;
sa = [0.0214,0.024,0.242,47.75];
d0 = 1;
d = sqrt(d0^2./10.^(-att/10));
fig = figure(5);
fig.Position = [2500 100 430 300];
clf; hold on;
ax1 = gca;
%hold on;
errorbar(ax1,att,fit_a,sa,'Color','k')
plot(ax1,att,fit_a,'k.')
xlabel('Attenuation [dB]')
ylabel('Angle Estimation [Degrees]')
grid on
ylim(ax1,[-1 1]*0.8)
xlim(ax1,[-5 60])
ax2 = axes('Position', ax1.Position, 'XAxisLocation', 'top', 'YColor', 'none','Color','none');
hold on;
errorbar(ax2,d,fit_a,sa,'Color','k')
plot(ax2,d,fit_a,'k.')
ylim(ax2,[-1 1]*0.8)
ax2.XScale = 'log';
xlim(ax2,sqrt(d0^2./10.^(-[-5 60]/10)))
xlabel(ax2,'Effective Distance [meters]','Position',[23.7138 0.65 -1]);
ax2.XColor = 'b';

%% B3 Obs
labs = {'B3 Obs, with homing',...
    'B3 Obs, with homing, reseat connector',...
    'B3 Obs, no homing',...
    'B3 Obs, no homing'};

ch = {'696 (Pol A)', '697 (Pol B)'};

resolution = 15;
lims = {[87 90.5], [-4, -1.5]};
for fileind = 1:4
    load(fullfile('dev',['rps_motor_check_', sprintf('%1d.mat',fileind)]))

    fig = figure(9+fileind);
    set(fig,'Position',[200 300 800 400])
    clf; hold on;
    for chind = 1:2
        ang = atand(tand(parms(:,1,chind)));
        subplot(1,2,chind)
        edges = lims{chind}(1):diff(lims{chind})/resolution:lims{chind}(2);
        N = histc(ang,edges);
        bar(edges,N,'histc')
        xlim(lims{chind})
        ylim([0,14])
        xlabel('Estimated Angle (^o)')
        ylabel('N')
        title({['Ch: ',ch{chind}];sprintf('Mean: %0.2f STD: %03f',nanmean(ang),nanstd(ang))})
        grid on
    end
    sgtitle(labs{fileind})
    saveas(fig,fullfile(figdir,sprintf('b3_obs_%01d.png',fileind)))
end

%% Grid angles from mill

load('z:/dev/grid_angles_from_mill.mat')
fig = figure(14);
fig.Position = [400 100 560 420];
resolution = 25;
lims = [0 1.1];
edges = lims(1):diff(lims)/resolution:lims(2);
N = histc(grid_angle,edges);
bar(edges,N,'histc')
xlim(lims)
ylim([0,20])
grid on
xlabel('Grid Horizontal')
ylabel('N')
title('Grid Calibration on knee mill')
text(0.05,33/2,sprintf('angle= %0.2f +/- %0.3f',nanmean(grid_angle),nanstd(grid_angle)))
text(0.05,30/2,sprintf('N Samples= %02d',length(grid_angle)))
%saveas(fig,fullfile(figdir,'isaac_cal_mill.png'))


%% B3 Obs
labs = {'B3 Obs, with homing',...
    'B3 Obs, with homing, reseat connector',...
    'B3 Obs, no homing',...
    'B3 Obs, no homing'};

ch = {'696 (Pol A)', '697 (Pol B)'};

resolution = 15;
lims = {[87 90.5], [-4, -1.5]};
for fileind = 1:4
    load(fullfile('dev',['rps_motor_check_', sprintf('%1d.mat',fileind)]))

    fig = figure(9+fileind);
    set(fig,'Position',[200 300 800 400])
    clf; hold on;
    for chind = 1:2
        ang = atand(tand(parms(:,1,chind)));
        subplot(1,2,chind)
        edges = lims{chind}(1):diff(lims{chind})/resolution:lims{chind}(2);
        N = histc(ang,edges);
        plot(atand(tand(parms(:,1,chind))))
        %bar(edges,N,'histc')
        %xlim(lims{chind})
        %ylim([0,14])
        xlabel('Estimated Angle (^o)')
        ylabel('N')
        title({['Ch: ',ch{chind}];sprintf('Mean: %0.2f STD: %03f',nanmean(ang),nanstd(ang))})
        grid on
    end
    sgtitle(labs{fileind})
    %    saveas(fig,fullfile(figdir,sprintf('b3_obs_%01d.png',fileind)))
end

%% Just look at timestreams

for prefind = 8%26%6:length(prefix)
    fig = figure(6);
    clf; hold on;
    d = dir(fullfile(direct,[prefix{prefind} '*']));

    data = {};
    for di = 1:length(d)
        if 0%prefind <= 22
            %disp(d(di).name)
            data{di} = load_lj_data(fullfile(direct,d(di).name));
        else
            %disp([prefix{prefind} sprintf('_scan_%i_loop_1.csv',di-1)])
            data{di} = load_lj_data(fullfile(direct,[prefix{prefind} ...
                sprintf('_scan_%i_loop_1.csv',di-1)]));
        end
    end

    for di = 1:length(data)
        val = 'Tilt';
        plot(data{di}.Time,data{di}.(val))


    end
    grid on

end

%% Isaac with/without homing and with/without backlash

fig = figure(116468);
fig.Position(3:4) = [440 470];
clf; hold on;
t = tiledlayout('flow');
t.Padding = 'compact';
t.TileSpacing = 'compact';

vals = [12, 20, 27, 34];
vals = [20,12,34,27];
valttls = {'No Backlash, No Homing','No Backlash, With Homing','With Backlash, No homing','With Backlash, With Homing'};
edges = linspace(-1,1,20)*0.2;


for valind = 1:length(vals)
nexttile()
idx = fd.schnum == vals(valind);
V = fd.phi(idx)-fd.phi_isaac(idx);
V = V-nanmean(V);
N = histc(V,edges);
S = nanstd(V);
L = length(V);
bar(edges,N,'histc');
grid on;
ttl = '';
ylbl = '';
if valind == 2
    ttl = 'No Homing';
    
elseif valind == 1
    ttl = 'With Homing';
    ylbl = 'No Backlash';
elseif valind ==3
    ylbl = 'With Backlash';
end

ylabel(ylbl)
title({ttl,...
    sprintf('S Dev: %0.3f | N: %i',S,L)...
    })
pbaspect([1 1 1])
if valind>2
xlabel('Pol Angle (mean-subtracted) [Deg]')
end
end

fname = 'isaac_homing_vs_backlash.png';
exportgraphics(fig,fullfile(figdir,fname),'Resolution',1200)


%% Look at mod curves. Again?
% No. Look at the actual timestreams.
scheds = [22 23]; % pre-obs
scheds = [12 20:27 30:35]; % in-lab D=0.5m
scheds = [12 20 21];% 27 28];

fig = figure(1);
fig.Position(3:4) = [1400 680];
clc
clf; hold on;
clear z
for obsind = scheds
    idx = find(fd.schnum==obsind);
    
    mc = [];
    ms = [];
    a = [];
    for ind = 1:length(idx)
        a(:, end+1) = angs{idx(ind)};
        mc(:,end+1) = modcurves{idx(ind)};%;./fd.A(idx(ind))/2;
        ms(:,end+1) = modstds{idx(ind)}; 
    end
    a = nanmean(a,2);
    %aind = inrange(abs(a),85,95);
    m1 = nanmedian(mc,2);
    m2 = nanmedian(ms,2);
    
    %s = nanstd(mc,[],2);
    subplot(2,1,1)
    hold on
    plot(NaN,NaN,'b')
    plot(NaN,NaN,'r')
    plot(NaN,NaN,'g')
    plot(NaN,NaN,'m')
    
    if obsind <20
        plot(a,m1,'b.-');
    elseif obsind <27
        plot(a,m1,'r.-');
    elseif obsind <34
        plot(a,m1,'g.-');
    else
        plot(a,m1,'m.-');
    end
    
    subplot(2,1,2)
    hold on
    plot(NaN,NaN,'b')
    plot(NaN,NaN,'r')
    plot(NaN,NaN,'g')
    plot(NaN,NaN,'m')
    
    if obsind <20
        plot(a,m2,'b.-');
    elseif obsind <27
        plot(a,m2,'r.-');
    elseif obsind <34
        plot(a,m2,'g.-');
    else
        plot(a,m2,'m.-');
    end
    
end
subplot(2,1,1)
grid on
%ylim([0 1.1])
xlim([-180 180])
legend({'Pole, No Homing','At Pole,W/ Homing','At Harvard, No Homing','At Harvard, With Homing'},'Location','northeastoutside')
xlabel('Command Angle')
ylabel({'Mod Curve Amplitudes (mV?)'})

subplot(2,1,2)
grid on
%ylim([0 0.012])
xlim([-180 180])
legend({'Pole, No Homing','At Pole,W/ Homing','At Harvard, No Homing','At Harvard, With Homing'},'Location','northeastoutside')
xlabel('Command Angle')
ylabel({'Mod Curve Timestream','S.Dev.'})

fname = 'isaac_modcurves_and_uncert_nonorm.png';
%exportgraphics(fig,fullfile(figdir,fname),'Resolution',1200)

%% Make a pager for histograms

fig = figure(723459);
    fig.Position(3:4) = [1000 400];
    clf; hold on;
t = tiledlayout(1,3);
t.Padding = 'compact';
t.TileSpacing = 'compact';

ttls = {'Old','New','New, No LNA'};
unqsch = unique(fd.schnum);
scheds = [34, 62, 46];
for schind = 1:length(scheds)
    idx = fd.schnum==scheds(schind);
    fdsch = structcut(fd,idx);
    dP = fdsch.phi-fdsch.phi_isaac;
    M = nanmean(dP);
    S = nanstd(dP);
    L = length(find(~isnan(dP)));
    
    nexttile()
    hold on
    lims = [-1 1]*0.4;
    edges = linspace(lims(1),lims(2),30);
    N = histc(dP,edges);
    bar(edges,N,'histc');
    %hist(dP)
    grid on
    xlim(lims)
    xlabel('$\phi_{fit}-\phi_{meas}$ [Deg]','FontSize',18)
    ylabel('N')
    t = datestr(nanmean(fdsch.time)/24/3600+datenum('1970-Jan-01:00:00:00','yyyy-mmm-dd:HH:MM:SS'));
    title({ttls{schind}, ...
        sprintf('%s UTC',t),...
        sprintf('M= %0.3f $|$ SDev= %0.3f $|$ N= %i',M,S,L)}, ...
        'FontSize',10)
    pbaspect([1 1 1])
    
end
fname = 'angle_hist_LNA_compare';
saveas(fig,fullfile(figdir,'hists',fname),'png')


%% Make a pager for histograms

unqsch = unique(fd.schnum);
for schind = 1:length(unqsch)
    idx = fd.schnum==unqsch(schind);
    fdsch = structcut(fd,idx);
    dP = fdsch.phi-fdsch.phi_isaac;
    M = nanmean(dP);
    S = nanstd(dP);
    L = length(find(~isnan(dP)));
    fig = figure(7385029);
    fig.Position(3:4) = [500 450];
    clf; hold on;

    hist(dP)
    grid on
    xlabel('$\phi_{fit}-\phi_{meas}$ [Deg]','FontSize',18)
    ylabel('N')
    t = datestr(nanmean(fdsch.time)/24/3600+datenum('1970-Jan-01:00:00:00','yyyy-mmm-dd:HH:MM:SS'));
    title({strrep(labs{unqsch(schind)},'_','\_'), ...
        sprintf('%s UTC',t),...
        sprintf('M= %0.3f $|$ SDev= %0.3f $|$ N= %i',M,S,L)})
    pbaspect([1 1 1])
    fname = sprintf('angle_hists_%i',unqsch(schind));
    saveas(fig,fullfile(figdir,'hists',fname),'png')
end


%% Make a pager for Mod Curve Residuals

unqsch = unique(fd.schnum);
for schind = 1:length(unqsch)
    idx = find(fd.schnum==unqsch(schind));
    fdsch = structcut(fd,idx);
    fig = figure(2875829);
    fig.Position(3:4) = [600 350];
    clf; hold on;
    
    for modind = 1:length(idx)
        plot(angs{idx(modind)},reses{idx(modind)}/fdsch.A(modind),'.-','MarkerSize',14)
        %plot(angs{idx(modind)},reses{idx(modind)},'.-','MarkerSize',14)
    end
    grid on
    ylabel('$A_{meas}-A_{model}$ (Volts)','FontSize',14)
    xlabel('Stage Angle WRT Gravity')
    t = datestr(nanmean(fdsch.time)/24/3600+datenum('1970-Jan-01:00:00:00','yyyy-mmm-dd:HH:MM:SS'));
    ttl = labs{unqsch(schind)};
    ttl = strrep(ttl,'_','\_');
    ttl = strrep(ttl,'%','\%');
    ttl = strrep(ttl,'<','$<$');
    title({ttl, ...
        sprintf('%s UTC',t)})
    ylim([-1 1]*0.1)
    fname = sprintf('reses_%i',unqsch(schind));
    saveas(fig,fullfile(figdir,'reses',fname),'png')

end

%% Make a pager for Mod Curves

unqsch = unique(fd.schnum);
for schind = 1:length(unqsch)
    idx = find(fd.schnum==unqsch(schind));
    fdsch = structcut(fd,idx);
    fig = figure(2875828);
    fig.Position(3:4) = [600 350];
    clf; hold on;
    
    for modind = 1:length(idx)
        plot(angs{idx(modind)},modcurves{idx(modind)},'.-','MarkerSize',14)
    end
    ylim([0 1]*4.5)
    %ylim([0 1]*0.5)
    grid on
    ylabel('$A_{meas}$ (Volts)','FontSize',14)
    xlabel('Stage Angle WRT Gravity')
    t = datestr(nanmean(fdsch.time)/24/3600+datenum('1970-Jan-01:00:00:00','yyyy-mmm-dd:HH:MM:SS'));
    ttl = labs{unqsch(schind)};
    ttl = strrep(ttl,'_','\_');
    ttl = strrep(ttl,'%','\%');
    ttl = strrep(ttl,'<','$<$');
    title({ttl, ...
        sprintf('%s UTC',t)})
    fname = sprintf('mod_curves_%i',unqsch(schind));
    saveas(fig,fullfile(figdir,'mod_curves',fname),'png')

end

%% Make a pager for diff mod curve

unqsch = unique(fd.schnum);
for schind = 1:length(unqsch)
    idx = find(fd.schnum==unqsch(schind));
    fdsch = structcut(fd,idx);
    fig = figure(2875829);
    fig.Position(3:4) = [600 350];
    clf; hold on;
    
    for modind = 1:length(idx)
        m0 = sort(modcurves{idx(modind)});
        plot(1:12,diff(modcurves{idx(modind)})/mean(m0((end-2):end)),'.-','MarkerSize',14)
    end
    ylim([-1 1]*0.5*1.1)
    grid on
    ylabel('$\Delta A_{meas}/A_0$ (Volts)','FontSize',14)
    %xlabel('Stage Angle WRT Gravity')
    t = datestr(nanmean(fdsch.time)/24/3600+datenum('1970-Jan-01:00:00:00','yyyy-mmm-dd:HH:MM:SS'));
    ttl = labs{unqsch(schind)};
    ttl = strrep(ttl,'_','\_');
    ttl = strrep(ttl,'%','\%');
    title({ttl, ...
        sprintf('%s UTC',t)})
    fname = sprintf('diff_mod_curves_%i',unqsch(schind));
    saveas(fig,fullfile(figdir,'mod_curves',fname),'png')

end

%% declare plot ylims
lims = [-0.65 0.4];
%% Plot dPhi vs alignment offset
clc
offs = [0 1 2 0 1 2 -1 -2 -1 0 0.5 0.25 -0.25 -0.5 -0.125];
schnums = [47:61];
%cmlines = distinguishable_colors(length(schnums));

fig = figure(173476);
fig.Position(3:4) = [700 250];
clf; hold on;

[offs_all dp_all,dp_mn] = deal([]);
for schind = 1:length(schnums)
    idx = find(fd.schnum==schnums(schind));
    O = repmat(offs(schind),1,length(idx));
    dP = fd.phi(idx)-fd.phi_isaac(idx);
    plot(O,dP,'.','Color',cmlines(schind,:),'MarkerSize',10)
    %plot(repmat(offs(schind),1,length(idx)),fd.tilt(idx),'x','Color',cmlines(schind,:),'MarkerSize',10)
    %plot(repmat(offs(schind),1,length(idx)),fd.istilt(idx),'^','Color',cmlines(schind,:),'MarkerSize',10)
    offs_all = [offs_all O];
    dp_all = [dp_all dP];
    dp_mn(schind) = nanmean(dP);
end
%plot(offs,dp_mn,'Color',cmlines(2,:))
grid on
xlim([-1 1]*2.5)
ylim(lims);
ylabel('$\phi_{fit}-\phi_{meas}$ [Deg]','FontSize',18)
xlabel('RPS Azimuthal Alignment Offset [Degrees]')
title({'Angle Bias vs. Alignment Offset',''})
fname = 'dp_vs_alignment';
saveas(fig,fullfile(figdir,fname),'png')

%% Plot dPhi vs alignment offset -- with constant amplitude
clc
offs = [0 1 2 1 0 -1 -2];
schnums = [63:69];
%cmlines = distinguishable_colors(length(schnums));
cmlines = colormap('lines');
fig = figure(174357);
fig.Position(3:4) = [700 250];
clf; hold on;

[offs_all dp_all,dp_mn] = deal([]);
for schind = 1:length(schnums)
    idx = find(fd.schnum==schnums(schind));
    O = repmat(offs(schind),1,length(idx));
    dP = fd.phi(idx)-fd.phi_isaac(idx);
    plot(O,dP,'.','Color',cmlines(schind,:),'MarkerSize',14)
    %plot(repmat(offs(schind),1,length(idx)),fd.tilt(idx),'x','Color',cmlines(schind,:),'MarkerSize',10)
    %plot(repmat(offs(schind),1,length(idx)),fd.istilt(idx),'^','Color',cmlines(schind,:),'MarkerSize',10)
    offs_all = [offs_all O];
    dp_all = [dp_all dP];
    dp_mn(schind) = nanmean(dP);
end
plot(offs,dp_mn,'Color',cmlines(2,:))
grid on
xlim([-1 1]*2.5)
ylim(lims);
ylabel('$\phi_{fit}-\phi_{meas}$ [Deg]','FontSize',18)
xlabel('RPS Azimuthal Alignment Offset [Degrees]')
title({'Angle Bias vs. Alignment Offset','Constant Peak Amplitude'})
fname = 'dp_vs_alignment_const_amp';
saveas(fig,fullfile(figdir,fname),'png')

%% Plot dPhi vs alignment offset -- peaked up distance 1m
clc
offs = [0 0.5 1 1.5 2 1 0 -0.5 -1 -2];
schnums = [70:79];
%cmlines = distinguishable_colors(length(schnums));
cmlines = colormap('lines');
fig = figure(173477);
fig.Position(3:4) = [700 250];
clf; hold on;

[offs_all dp_all,dp_mn] = deal([]);
for schind = 1:length(schnums)
    idx = find(fd.schnum==schnums(schind));
    O = repmat(offs(schind),1,length(idx));
    dP = fd.phi(idx)-fd.phi_isaac(idx);
    plot(O,dP,'.','Color',cmlines(1,:),'MarkerSize',14)
    %plot(repmat(offs(schind),1,length(idx)),fd.tilt(idx),'x','Color',cmlines(schind,:),'MarkerSize',10)
    %plot(repmat(offs(schind),1,length(idx)),fd.istilt(idx),'^','Color',cmlines(schind,:),'MarkerSize',10)
    offs_all = [offs_all O];
    dp_all = [dp_all dP];
    dp_mn(schind) = nanmean(dP);
end
plot(offs_all,dp_all,'Color',cmlines(2,:))
grid on
xlim([-1 1]*2.5)
ylim(lims);
ylabel('$\phi_{fit}-\phi_{meas}$ [Deg]','FontSize',18)
xlabel('RPS Azimuthal Alignment Offset [Degrees]')
title({'Angle Bias vs. Alignment Offset, Distance of $\sim1$m','RPS/ISAAC Peaked up to $<1\%$ of max amp'})
fname = 'dp_vs_alignment_new_dist';
%saveas(fig,fullfile(figdir,fname),'png')

%% Plot dPhi vs alignment offset -- peaked up distance 1m -- ISAAC Moving
clc
offs = [0 0.5 1 1.5 2 -2 -1.5 -1 -0.5 0];
schnums = [80:89];
%cmlines = distinguishable_colors(length(schnums));
cmlines = colormap('lines');
fig = figure(173477);
fig.Position(3:4) = [700 250];
clf; hold on;

[offs_all dp_all,dp_mn] = deal([]);
for schind = 1:length(schnums)
    idx = find(fd.schnum==schnums(schind));
    O = repmat(offs(schind),1,length(idx));
    dP = fd.phi(idx)-fd.phi_isaac(idx);
    plot(O,dP,'.','Color',cmlines(1,:),'MarkerSize',14)
    %plot(repmat(offs(schind),1,length(idx)),fd.tilt(idx),'x','Color',cmlines(schind,:),'MarkerSize',10)
    %plot(repmat(offs(schind),1,length(idx)),fd.istilt(idx),'^','Color',cmlines(schind,:),'MarkerSize',10)
    offs_all = [offs_all O];
    dp_all = [dp_all dP];
    dp_mn(schind) = nanmean(dP);
end
plot(offs_all,dp_all,'Color',cmlines(2,:))
grid on
xlim([-1 1]*2.5)
ylim(lims);
ylabel('$\phi_{fit}-\phi_{meas}$ [Deg]','FontSize',18)
xlabel('ISAAC Azimuthal Alignment Offset [Degrees]')
title({'Angle Bias vs. Alignment Offset, Distance of $\sim1$m','ISAAC Moving'})
fname = 'dp_vs_alignment_isaac_move';
saveas(fig,fullfile(figdir,fname),'png')

%% Dphi vs eccosorb test

schnums = [90:97];
idx = ismember(fd.schnum,schnums);
dP = fd.phi(idx)-fd.phi_isaac(idx);

fig = figure(3897532);
fig.Position(3:4) = [700 250];
clf; hold on;
plot(1:length(dP),dP,'.','MarkerSize',14)
plot(1:length(dP),dP)
ylim(lims);
xlim([0 9])
grid on
xlabel('Measurement Number')
ylabel('$\phi_{fit}-\phi_{meas}$ [Deg]','FontSize',18)
title({'Angle Bias vs. Eccosorb',''})
fname = 'dp_vs_alignment_eccosorb';
saveas(fig,fullfile(figdir,fname),'png')

%% Plot dPhi vs alignment offset -- Horn Flipped 180 deg
clc
offs = [0 0.5 1 1.5 2 -2 -1.5 -1 -0.5 0];
schnums = [98:107];
%cmlines = distinguishable_colors(length(schnums));
cmlines = colormap('lines');
fig = figure(173477);
fig.Position(3:4) = [700 250];
clf; hold on;

[offs_all dp_all,dp_mn] = deal([]);
for schind = 1:length(schnums)
    idx = find(fd.schnum==schnums(schind));
    O = repmat(offs(schind),1,length(idx));
    dP = fd.phi(idx)-fd.phi_isaac(idx);
    plot(O,dP,'.','Color',cmlines(1,:),'MarkerSize',14)
    %plot(repmat(offs(schind),1,length(idx)),fd.tilt(idx),'x','Color',cmlines(schind,:),'MarkerSize',10)
    %plot(repmat(offs(schind),1,length(idx)),fd.istilt(idx),'^','Color',cmlines(schind,:),'MarkerSize',10)
    offs_all = [offs_all O];
    dp_all = [dp_all dP];
    dp_mn(schind) = nanmean(dP);
end
plot(offs_all,dp_all,'Color',cmlines(2,:))
grid on
xlim([-1 1]*2.5)
ylim(lims);
ylabel('$\phi_{fit}-\phi_{meas}$ [Deg]','FontSize',18)
xlabel('RPS Azimuthal Alignment Offset [Degrees]')
title({'Angle Bias vs. Alignment Offset, Distance of $\sim1$m','ISAAC Horn Flipped 180deg'})
fname = 'dp_vs_alignment_horn_flipped';
saveas(fig,fullfile(figdir,fname),'png')

%% Plot dPhi vs alignment offset -- New, small horn
clc
offs = [0 0.5 1 1.5 2 -2:0.5:2 0];
schnums = [108:122];
%cmlines = distinguishable_colors(length(schnums));
cmlines = colormap('lines');
fig = figure(173478);
fig.Position(3:4) = [700 250];
clf; hold on;

[offs_all dp_all,dp_mn] = deal([]);
for schind = 1:length(schnums)
    idx = find(fd.schnum==schnums(schind));
    O = repmat(offs(schind),1,length(idx));
    dP = fd.phi(idx)-fd.phi_isaac(idx);
    plot(O,dP,'.','Color',cmlines(1,:),'MarkerSize',14)
    %plot(repmat(offs(schind),1,length(idx)),fd.tilt(idx),'x','Color',cmlines(schind,:),'MarkerSize',10)
    %plot(repmat(offs(schind),1,length(idx)),fd.istilt(idx),'^','Color',cmlines(schind,:),'MarkerSize',10)
    offs_all = [offs_all O];
    dp_all = [dp_all dP];
    dp_mn(schind) = nanmean(dP);
end
plot(offs_all,dp_all,'Color',cmlines(2,:))
grid on
xlim([-1 1]*2.5)
ylim(lims);
ylabel('$\phi_{fit}-\phi_{meas}$ [Deg]','FontSize',18)
xlabel('RPS Azimuthal Alignment Offset [Degrees]')
title({'Angle Bias vs. Alignment Offset, Distance of $\sim1$m','ISAAC Horn Replaced With Small Horn'})
fname = 'dp_vs_alignment_small_horn';
saveas(fig,fullfile(figdir,fname),'png')

%% Plot dPhi vs alignment offset -- Small horn, flipped
clc
offs = [0 0.5 1 1.5 2 -2:0.5:2 0];
schnums = [123:137];
%cmlines = distinguishable_colors(length(schnums));
cmlines = colormap('lines');
fig = figure(173478);
fig.Position(3:4) = [700 250];
clf; hold on;

[offs_all dp_all,dp_mn] = deal([]);
for schind = 1:length(schnums)
    idx = find(fd.schnum==schnums(schind));
    O = repmat(offs(schind),1,length(idx));
    dP = fd.phi(idx)-fd.phi_isaac(idx);
    plot(O,dP,'.','Color',cmlines(1,:),'MarkerSize',14)
    %plot(repmat(offs(schind),1,length(idx)),fd.tilt(idx),'x','Color',cmlines(schind,:),'MarkerSize',10)
    %plot(repmat(offs(schind),1,length(idx)),fd.istilt(idx),'^','Color',cmlines(schind,:),'MarkerSize',10)
    offs_all = [offs_all O];
    dp_all = [dp_all dP];
    dp_mn(schind) = nanmean(dP);
end
plot(offs_all,dp_all,'Color',cmlines(2,:))
grid on
xlim([-1 1]*2.5)
ylim(lims);
ylabel('$\phi_{fit}-\phi_{meas}$ [Deg]','FontSize',18)
xlabel('RPS Azimuthal Alignment Offset [Degrees]')
title({'Angle Bias vs. Alignment Offset, Distance of $\sim1$m','ISAAC Horn Replaced With Small Horn, flipped'})
fname = 'dp_vs_alignment_small_horn';
saveas(fig,fullfile(figdir,fname),'png')

%% Plot dPhi vs alignment offset -- Small horn, flipped, rotating ISAAC
clc
offs = [0 0.5 1 1.5 2 -2:0.5:2 0 -0.5 -1 -1.5 -2 0];
schnums = [138:157];
%cmlines = distinguishable_colors(length(schnums));
cmlines = colormap('lines');
fig = figure(173477);
fig.Position(3:4) = [700 250];
clf; hold on;

[offs_all dp_all,dp_mn] = deal([]);
for schind = 1:length(schnums)
    idx = find(fd.schnum==schnums(schind));
    O = repmat(offs(schind),1,length(idx));
    dP = fd.phi(idx)-fd.phi_isaac(idx);
    plot(O,dP,'.','Color',cmlines(1,:),'MarkerSize',14)
    %plot(repmat(offs(schind),1,length(idx)),fd.tilt(idx),'x','Color',cmlines(schind,:),'MarkerSize',10)
    %plot(repmat(offs(schind),1,length(idx)),fd.istilt(idx),'^','Color',cmlines(schind,:),'MarkerSize',10)
    offs_all = [offs_all O];
    dp_all = [dp_all dP];
    dp_mn(schind) = nanmean(dP);
end
plot(offs_all,dp_all,'Color',cmlines(2,:))
grid on
xlim([-1 1]*2.5)
ylim(lims);
ylabel('$\phi_{fit}-\phi_{meas}$ [Deg]','FontSize',18)
xlabel('ISAAC Azimuthal Alignment Offset [Degrees]')
title({'Angle Bias vs. Alignment Offset, Distance of $\sim1$m','15dB horn, flipped, rotating ISAAC'})
fname = 'dp_vs_alignment_small_horn_rot_isaac_allangles';
saveas(fig,fullfile(figdir,fname),'png')

%% Plot translation vs alignment offset -- Small horn, flipped, moving ISAAC
clc
offs = [0.0 -0.65 -1.33 -2.0 -1.34 -0.67 0 0.67 1.33 1.99 0];
schnums = [158:168];
%cmlines = distinguishable_colors(length(schnums));
cmlines = colormap('lines');
fig = figure(173477);
fig.Position(3:4) = [700 250];
clf; hold on;

[offs_all dp_all,dp_mn] = deal([]);
for schind = 1:length(schnums)
    idx = find(fd.schnum==schnums(schind));
    O = repmat(offs(schind),1,length(idx));
    dP = fd.phi(idx)-fd.phi_isaac(idx);
    plot(O,dP,'.','Color',cmlines(1,:),'MarkerSize',14)
    %plot(repmat(offs(schind),1,length(idx)),fd.tilt(idx),'x','Color',cmlines(schind,:),'MarkerSize',10)
    %plot(repmat(offs(schind),1,length(idx)),fd.istilt(idx),'^','Color',cmlines(schind,:),'MarkerSize',10)
    offs_all = [offs_all O];
    dp_all = [dp_all dP];
    dp_mn(schind) = nanmean(dP);
end
plot(offs_all,dp_all,'Color',cmlines(2,:))
grid on
xlim([-1 1]*2.5)
ylim(lims);
ylabel('$\phi_{fit}-\phi_{meas}$ [Deg]','FontSize',18)
xlabel('Equivalent RPS angle [Degrees]')
title({'Angle Bias vs. Translation Offset, Distance of $\sim1$m','15dB horn, flipped, moving ISAAC'})
fname = 'dp_vs_alignment_moving_isaac';
saveas(fig,fullfile(figdir,fname),'png')

%% Plot dPhi vs alignment offset -- Removed RPS shroud
clc
offs = [0 0.5 1 1.5 2 -2:0.5:2 0];
schnums = [169:183];
%cmlines = distinguishable_colors(length(schnums));
cmlines = colormap('lines');
fig = figure(173478);
fig.Position(3:4) = [700 250];
clf; hold on;

[offs_all dp_all,dp_mn] = deal([]);
for schind = 1:length(schnums)
    idx = find(fd.schnum==schnums(schind));
    O = repmat(offs(schind),1,length(idx));
    dP = fd.phi(idx)-fd.phi_isaac(idx);
    plot(O,dP,'.','Color',cmlines(1,:),'MarkerSize',14)
    %plot(repmat(offs(schind),1,length(idx)),fd.tilt(idx),'x','Color',cmlines(schind,:),'MarkerSize',10)
    %plot(repmat(offs(schind),1,length(idx)),fd.istilt(idx),'^','Color',cmlines(schind,:),'MarkerSize',10)
    offs_all = [offs_all O];
    dp_all = [dp_all dP];
    dp_mn(schind) = nanmean(dP);
end
plot(offs_all,dp_all,'Color',cmlines(2,:))
grid on
xlim([-1 1]*2.5)
ylim(lims-6.9)
ylabel('$\phi_{fit}-\phi_{meas}$ [Deg]','FontSize',18)
xlabel('RPS Azimuthal Alignment Offset [Degrees]')
title({'Angle Bias vs. Alignment Offset, Distance of $\sim1$m','RPS shroud and wire grid removed'})
fname = 'dp_vs_alignment_no_shroud';
saveas(fig,fullfile(figdir,fname),'png')

%% Plot dPhi vs alignment offset -- Removed RPS shroud again
clc
offs = [0 0.5 1 1.5 2 -2:0.5:2 0];
schnums = [184:198];
%cmlines = distinguishable_colors(length(schnums));
cmlines = colormap('lines');
fig = figure(173478);
fig.Position(3:4) = [700 250];
clf; hold on;

[offs_all dp_all,dp_mn] = deal([]);
for schind = 1:length(schnums)
    idx = find(fd.schnum==schnums(schind));
    O = repmat(offs(schind),1,length(idx));
    dP = fd.phi(idx)-fd.phi_isaac(idx);
    plot(O,dP,'.','Color',cmlines(1,:),'MarkerSize',14)
    %plot(repmat(offs(schind),1,length(idx)),fd.tilt(idx),'x','Color',cmlines(schind,:),'MarkerSize',10)
    %plot(repmat(offs(schind),1,length(idx)),fd.istilt(idx),'^','Color',cmlines(schind,:),'MarkerSize',10)
    offs_all = [offs_all O];
    dp_all = [dp_all dP];
    dp_mn(schind) = nanmean(dP);
end
plot(offs_all,dp_all,'Color',cmlines(2,:))
grid on
xlim([-1 1]*2.5)
ylim(lims-6.9)
ylabel('$\phi_{fit}-\phi_{meas}$ [Deg]','FontSize',18)
xlabel('RPS Azimuthal Alignment Offset [Degrees]')
title({'Angle Bias vs. Alignment Offset, Distance of $\sim1$m','RPS shroud and wire grid removed'})
fname = 'dp_vs_alignment_no_shroud_pt2';
saveas(fig,fullfile(figdir,fname),'png')

%% Plot dPhi vs alignment offset -- Moved to 53in
clc
offs = [0 0.5 1 1.5 2 -2:0.5:2 0];
schnums = [199:213];
%cmlines = distinguishable_colors(length(schnums));
cmlines = colormap('lines');
fig = figure(173478);
fig.Position(3:4) = [700 250];
clf; hold on;

[offs_all dp_all,dp_mn] = deal([]);
for schind = 1:length(schnums)
    idx = find(fd.schnum==schnums(schind));
    O = repmat(offs(schind),1,length(idx));
    dP = fd.phi(idx)-fd.phi_isaac(idx);
    plot(O,dP,'.','Color',cmlines(1,:),'MarkerSize',14)
    %plot(repmat(offs(schind),1,length(idx)),fd.tilt(idx),'x','Color',cmlines(schind,:),'MarkerSize',10)
    %plot(repmat(offs(schind),1,length(idx)),fd.istilt(idx),'^','Color',cmlines(schind,:),'MarkerSize',10)
    offs_all = [offs_all O];
    dp_all = [dp_all dP];
    dp_mn(schind) = nanmean(dP);
end
plot(offs_all,dp_all,'Color',cmlines(2,:))
grid on
xlim([-1 1]*2.5)
ylim(lims);
ylabel('$\phi_{fit}-\phi_{meas}$ [Deg]','FontSize',18)
xlabel('RPS Azimuthal Alignment Offset [Degrees]')
title({'Angle Bias vs. Alignment Offset, Distance of $\sim1.3$m','15dB horn, flipped'})
fname = 'dp_vs_alignment_1p3meters';
saveas(fig,fullfile(figdir,fname),'png')

%% Plot dPhi vs alignment offset -- Rolled ISAAC by 90deg
clc
offs = [0 0.5 1 1.5 2 -2:0.5:2 0];
schnums = [214:228];
%cmlines = distinguishable_colors(length(schnums));
cmlines = colormap('lines');
fig = figure(173478);
fig.Position(3:4) = [700 250];
clf; hold on;

[offs_all dp_all,dp_mn] = deal([]);
for schind = 1:length(schnums)
    idx = find(fd.schnum==schnums(schind));
    O = repmat(offs(schind),1,length(idx));
    dP = fd.phi(idx)-fd.phi_isaac(idx);
    plot(O,dP,'.','Color',cmlines(1,:),'MarkerSize',14)
    %plot(repmat(offs(schind),1,length(idx)),fd.tilt(idx),'x','Color',cmlines(schind,:),'MarkerSize',10)
    %plot(repmat(offs(schind),1,length(idx)),fd.istilt(idx),'^','Color',cmlines(schind,:),'MarkerSize',10)
    offs_all = [offs_all O];
    dp_all = [dp_all dP];
    dp_mn(schind) = nanmean(dP);
end
plot(offs_all,dp_all,'Color',cmlines(2,:))
grid on
xlim([-1 1]*2.5)
ylim(lims+87.5);
ylabel('$\phi_{fit}-\phi_{meas}$ [Deg]','FontSize',18)
xlabel('RPS Azimuthal Alignment Offset [Degrees]')
title({'Angle Bias vs. Alignment Offset, Distance of $\sim1$m','ISAAC flipped 90deg'})
fname = 'dp_vs_alignment_isaac_roll90';
saveas(fig,fullfile(figdir,fname),'png')

%% Plot dPhi vs alignment offset -- Rolled ISAAC by 90deg
clc
offs = [0 0.5 1 1.5 2 -2:0.5:2 0];
schnums = [229:243];
%cmlines = distinguishable_colors(length(schnums));
cmlines = colormap('lines');
fig = figure(173478);
fig.Position(3:4) = [700 250];
clf; hold on;

[offs_all dp_all,dp_mn] = deal([]);
for schind = 1:length(schnums)
    idx = find(fd.schnum==schnums(schind));
    O = repmat(offs(schind),1,length(idx));
    dP = fd.phi(idx)-fd.phi_isaac(idx);
    plot(O,dP,'.','Color',cmlines(1,:),'MarkerSize',14)
    %plot(repmat(offs(schind),1,length(idx)),fd.tilt(idx),'x','Color',cmlines(schind,:),'MarkerSize',10)
    %plot(repmat(offs(schind),1,length(idx)),fd.istilt(idx),'^','Color',cmlines(schind,:),'MarkerSize',10)
    offs_all = [offs_all O];
    dp_all = [dp_all dP];
    dp_mn(schind) = nanmean(dP);
end
plot(offs_all,dp_all,'Color',cmlines(2,:))
grid on
xlim([-1 1]*2.5)
ylim(lims);
ylabel('$\phi_{fit}-\phi_{meas}$ [Deg]','FontSize',18)
xlabel('RPS Azimuthal Alignment Offset [Degrees]')
title({'Angle Bias vs. Alignment Offset, Distance of $\sim2.4$m','15dB horn, flipped'})
fname = 'dp_vs_alignment_2p4meters';
saveas(fig,fullfile(figdir,fname),'png')

%% Combine Mod curves and Compare as-fit residuals 

clc
unqsch = unique(fd.schnum);
for schind = 1:length(unqsch)
    idx = find(fd.schnum==unqsch(schind));
    fdsch = structcut(fd,idx);
    dP = fd.phi(idx)-fd.phi_isaac(idx);
    fig = figure(2875828);
    fig.Position(3:4) = [800 550];
    clf; hold on;
    t = tiledlayout(2,1);
    t.Padding = 'compact';
    t.TileSpacing = 'compact';
    
    nexttile(1); hold on;
    for modind = 1:length(idx)
        plot(angs{idx(modind)},modcurves{idx(modind)},'.-','MarkerSize',14)
    end
    ylim([0 1]*4.5)
    %ylim([0 1]*0.5)
    grid on
    ylabel('$A_{meas}$ (Volts)','FontSize',14)
    %xlabel('Stage Angle WRT Gravity')
    t = datestr(nanmean(fdsch.time)/24/3600+datenum('1970-Jan-01:00:00:00','yyyy-mmm-dd:HH:MM:SS'));
    ttl = labs{unqsch(schind)};
    ttl = strrep(ttl,'_','\_');
    ttl = strrep(ttl,'%','\%');
    ttl = strrep(ttl,'<','$<$');
    title({ttl, ...
        sprintf('%s UTC',t)})

    nexttile(2); hold on;
    for modind = 1:length(idx)
        parm0 = [fdsch.phi_isaac(modind) fdsch.eps(modind) fdsch.N1(modind) fdsch.N2(modind) fdsch.A(modind)];
        idealres = rps_get_mod_model(parm0,angs{idx(modind)})...
            - rps_get_mod_model(parm0-[dP(modind) 0 0 0 0],angs{idx(modind)});
        plot(angs{idx(modind)},idealres/fdsch.A(modind),'.-','MarkerSize',14,'Color',cmlines(3,:))
        fixedres = modcurves{idx(modind)} - rps_get_mod_model(parm0,angs{idx(modind)});
        plot(angs{idx(modind)},fixedres/fdsch.A(modind),'.-','MarkerSize',14,'Color',cmlines(2,:))
        plot(angs{idx(modind)},reses{idx(modind)}/fdsch.A(modind),'.-','MarkerSize',14,'Color',cmlines(1,:))
    end
    ylim([-1 1]*0.03)
    %ylim([0 1]*0.5)
    legend({'Ideal Offset','Angle Fixed','Real Residuals'},'Location','northeastoutside')
    grid on
    ylabel('$(A_{meas}-A_{model})/A_{max}$','FontSize',10)
    xlabel('Stage Angle WRT Gravity')
    
    fname = sprintf('mod_curve_and_residuals_%i',unqsch(schind));
    saveas(fig,fullfile(figdir,'mod_curves',fname),'png')

end


%% Make plots for tilt meters
clc
unqsch = unique(fd.schnum);
for schind = 1:length(unqsch)
    idx = find(fd.schnum==unqsch(schind));
    fdsch = structcut(fd,idx);
    fig = figure(2877389);
    fig.Position(3:4) = [600 350];
    clf; hold on;
    t = tiledlayout(2,1);
    t.Padding = 'compact';
    t.TileSpacing = 'compact';

    nexttile(1)
    for modind = 1:length(idx)
        plot(angs{idx(modind)},tilts{idx(modind)},'.-','MarkerSize',14)
    end
    %ylim([0 1]*4.5)
    %ylim([0 1]*0.5)
    grid on
    ylabel('RPS Tilt','FontSize',14)
    %xlabel('Stage Angle WRT Gravity')
    t = datestr(nanmean(fdsch.time)/24/3600+datenum('1970-Jan-01:00:00:00','yyyy-mmm-dd:HH:MM:SS'));
    ttl = labs{unqsch(schind)};
    ttl = strrep(ttl,'_','\_');
    ttl = strrep(ttl,'%','\%');
    ttl = strrep(ttl,'<','$<$');
    title({ttl, ...
        sprintf('%s UTC',t)})
    
    nexttile(2)
    for modind = 1:length(idx)
        plot(angs{idx(modind)},istilts{idx(modind)},'.-','MarkerSize',14)
    end
    %ylim([0 1]*4.5)
    %ylim([0 1]*0.5)
    grid on
    ylabel('ISAAC Tilt','FontSize',14)
    xlabel('Stage Angle WRT Gravity')
    
    fname = sprintf('tilts_%i',unqsch(schind));
    saveas(fig,fullfile(figdir,'tilts',fname),'png')

end

%end % End of Main function
%%
function [prefix, labs] = get_pref_and_labs()

prefix = {...
    'isaac_cal_1_with_home_';... %0dB
    'isaac_cal_1_no_home_';... %0dB
    'isaac_cal_2_';... %-30dB
    'isaac_cal_3_';... %-45dB
    'isaac_cal_4_';... %-55dB
    'isaac_cal_1_new_home';... % New homing, Indoor
    'isaac_cal_5_no_home';... % No Homing, Indoor
    'isaac_cal_6_no_home';... % placed outside
    'isaac_cal_7_no_home';... % Moved to same level
    'isaac_cal_8_no_home';... % outside relevel
    'isaac_cal_9_no_home';... % Outside after grid check
    'isaac_cal_10_no_home';... % Homed the RPS. Checking again.
    'isaac_cal_2_with_home';... % Cal with constant homing
    'isaac_cal_3_with_home';... % After bringing it back down and fixing recessed pin.
    'isaac_cal_4_with_home_indoor';... % Indoor test again
    'isaac_cal_5_with_home_indoor';... % Indoor test after installing heater
    'isaac_cal_6_with_home_indoor';... % Indoor test because I'm paranoid
    'isaac_cal_7_with_home_indoor';... % Another indoor test because everything is still shitty.
    'isaac_cal_8_with_home_indoor';... % Added a 'fine' homing switch to the stepper motor. Skip this one.
    'isaac_cal_1_fine_home_check';... % Edited the homing to avoid intermittent homing position.
    'isaac_cal_2_fine_home_check';... % Checking if I'm an idiot...
    'isaac_cal_1_fine_home_check_outdoor';... % Outdoor Check!
    'isaac_cal_1_fine_home_check_outdoor_post';... % Outdoor Check, post obs
    'isaac_cal_1_fine_home_check_indoor_post';... % Indoor Check, post obs
    'isaac_cal_2_fine_home_check_indoor_post';... % Indoor Check, post obs, ISAAC refrozen
    'isaac_cal_3_fine_home_check_indoor_post';... % Indoor Check, post obs, RPS refrozen
    'isaac_cal_jig_nonhome_1';... % RPS jig check 25 Jun 2022 @ 20"
    'isaac_cal_jig_nonhome_2';... % RPS jig check 25 Jun 2022 @ 40"
    'isaac_cal_jig_nonhome_3';... % RPS jig check 28 Jun 2022 @ 40" no changes
    'isaac_cal_jig_nonhome_4';... % RPS jig check 28 Jun 2022 @ 40", changed tilt
    'isaac_cal_jig_nonhome_5';... % RPS jig check 28 Jun 2022 @ 40", back to orig. tilt
    'isaac_cal_jig_nonhome_6';... % RPS jig check 28 Jun 2022 @ 20",
    'isaac_cal_jig_nonhome_7';... % RPS jig check 28 Jun 2022 @ 20", repeat of previous
    'isaac_cal_jig_nonhome_8';... % RPS jig check 29 Jun 2022 @ 20", This was actually homed!
    'isaac_cal_jig_withhome_2';... % 30 Jun 2022 @ 20"; homed. Repeat of last
    'isaac_cal_jig_withhome_3';... % 30 Jun 2022 @ 20"; homed. Repeat of last
    'isaac_cal_jig_latmove_1in_1';... %30 Jun 2022 @ 20"; homed. Offset 1in.
    'isaac_cal_jig_withhome_4';... % 29 Jul 2022 @ 37"; homed.
    'isaac_cal_jig_withhome_5';... % 31 Jul 2022 @ 37"; homed. Repeat of last
    'isaac_cal_jig_withhome_6';... % 31 Jul 2022 @ 37"; homed. Zeroed RPS to -0.03
    'isaac_cal_jig_withhome_7';... % 31 Jul 2022 @ 37"; homed. Zeroed ISAAC
    'isaac_cal_jig_withhome_8';... % 31 Jul 2022 @ 37"; homed. Added baffling.
    'isaac_cal_jig_offsets_1';... 28 Jul 2022 @ 40; homed. Ella moving around at 40"
    'isaac_cal_jig_july23_nohome_1';... 03 July 2023 @ 26". Removed Isaac LNA, added absorber.
    'isaac_cal_jig_july23_nohome_2';... 04 July 2023 @ 26". Removed Isaac LNA, added absorber.
    'isaac_cal_jig_july23_with_home_1';... 04 July 2023 @ 26". Homing
    'isaac_cal_jig_july23_with_home_aligncheck_0_1';... % Moving through a range of alignments on the RPS
    'isaac_cal_jig_july23_with_home_aligncheck_1deg_1';...
    'isaac_cal_jig_july23_with_home_aligncheck_2deg_1';...
    'isaac_cal_jig_july23_with_home_aligncheck_0deg_2';...
    'isaac_cal_jig_july23_with_home_aligncheck_1deg_2';...
    'isaac_cal_jig_july23_with_home_aligncheck_2deg_2';...
    'isaac_cal_jig_july23_with_home_aligncheck_m1deg_1';...
    'isaac_cal_jig_july23_with_home_aligncheck_m2deg_1';...
    'isaac_cal_jig_july23_with_home_aligncheck_m1deg_2';...
    'isaac_cal_jig_july23_with_home_aligncheck_0deg_3';...
    'isaac_cal_jig_july23_with_home_aligncheck_0p5deg_1';...
    'isaac_cal_jig_july23_with_home_aligncheck_0p25deg_1';...
    'isaac_cal_jig_july23_with_home_aligncheck_m0p25deg_1';...
    'isaac_cal_jig_july23_with_home_aligncheck_m0p5deg_1';...
    'isaac_cal_jig_july23_with_home_aligncheck_m0p125deg_1';...
    'isaac_cal_jig_july23_with_home_with_LNA_1';... Reinstalled the LNA to look at noise again.
    'isaac_cal_jig_july23_constamp_0deg_1';... % Alignment test, but this time I kept a constant peak amplitude
    'isaac_cal_jig_july23_constamp_1deg_1';...
    'isaac_cal_jig_july23_constamp_2deg_1';...
    'isaac_cal_jig_july23_constamp_1deg_2';...
    'isaac_cal_jig_july23_constamp_0deg_2';...
    'isaac_cal_jig_july23_constamp_m1deg_1';...
    'isaac_cal_jig_july23_constamp_m2deg_1';...
    'isaac_cal_jig_july23_with_homing_all_peaked_0deg_1';...
    'isaac_cal_jig_july23_with_homing_all_peaked_0p5deg_1';...
    'isaac_cal_jig_july23_with_homing_all_peaked_1deg_1';...
    'isaac_cal_jig_july23_with_homing_all_peaked_1p5deg_1';...
    'isaac_cal_jig_july23_with_homing_all_peaked_2deg_1';...
    'isaac_cal_jig_july23_with_homing_all_peaked_1deg_2';...
    'isaac_cal_jig_july23_with_homing_all_peaked_0deg_2';...
    'isaac_cal_jig_july23_with_homing_all_peaked_m0p5deg_1';...
    'isaac_cal_jig_july23_with_homing_all_peaked_m1deg_1';...
    'isaac_cal_jig_july23_with_homing_all_peaked_m2deg_1';...
    'isaac_cal_jig_july23_with_homing_isaac_align_0deg_1';...
    'isaac_cal_jig_july23_with_homing_isaac_align_0p5deg_1';...
    'isaac_cal_jig_july23_with_homing_isaac_align_1deg_1';...
    'isaac_cal_jig_july23_with_homing_isaac_align_1p5deg_1';...
    'isaac_cal_jig_july23_with_homing_isaac_align_2deg_1';...
    'isaac_cal_jig_july23_with_homing_isaac_align_m2deg_1';...
    'isaac_cal_jig_july23_with_homing_isaac_align_m1p5deg_1';...
    'isaac_cal_jig_july23_with_homing_isaac_align_m1deg_1';...
    'isaac_cal_jig_july23_with_homing_isaac_align_m0p5deg_1';...
    'isaac_cal_jig_july23_with_homing_isaac_align_0deg_2';...
    'isaac_cal_jig_july23_with_homing_eccosorb_no_absorb_1';...
    'isaac_cal_jig_july23_with_homing_eccosorb_behind_isaac_1';...
    'isaac_cal_jig_july23_with_homing_eccosorb_in_window_1';...
    'isaac_cal_jig_july23_with_homing_eccosorb_on_table_1';...
    'isaac_cal_jig_july23_with_homing_eccosorb_table_rps_up_1';...
    'isaac_cal_jig_july23_with_homing_eccosorb_table_isaac_up_1';...
    'isaac_cal_jig_july23_with_homing_eccosorb_fully_blocked_1';...
    'isaac_cal_jig_july23_with_homing_eccosorb_both_sides_1';...
    'isaac_cal_jig_july23_with_homing_horn_flipped_0deg_1';...
    'isaac_cal_jig_july23_with_homing_horn_flipped_0p5deg_1';...
    'isaac_cal_jig_july23_with_homing_horn_flipped_1deg_1';...
    'isaac_cal_jig_july23_with_homing_horn_flipped_1p5deg_1';...
    'isaac_cal_jig_july23_with_homing_horn_flipped_2deg_1';...
    'isaac_cal_jig_july23_with_homing_horn_flipped_m2deg_1';...
    'isaac_cal_jig_july23_with_homing_horn_flipped_m1p5deg_1';...
    'isaac_cal_jig_july23_with_homing_horn_flipped_m1deg_1';...
    'isaac_cal_jig_july23_with_homing_horn_flipped_m0p5deg_1';...
    'isaac_cal_jig_july23_with_homing_horn_flipped_0deg_2';...
    'isaac_cal_jig_july23_with_homing_small_horn_0deg_1';...
    'isaac_cal_jig_july23_with_homing_small_horn_0p5deg_1';...
    'isaac_cal_jig_july23_with_homing_small_horn_1deg_1';...
    'isaac_cal_jig_july23_with_homing_small_horn_1p5deg_1';...
    'isaac_cal_jig_july23_with_homing_small_horn_2deg_1';...
    'isaac_cal_jig_july23_with_homing_small_horn_m2deg_1';...
    'isaac_cal_jig_july23_with_homing_small_horn_m1p5deg_1';...
    'isaac_cal_jig_july23_with_homing_small_horn_m1deg_1';...
    'isaac_cal_jig_july23_with_homing_small_horn_m0p5deg_1';...
    'isaac_cal_jig_july23_with_homing_small_horn_0deg_2';...
    'isaac_cal_jig_july23_with_homing_small_horn_0p5deg_2';...
    'isaac_cal_jig_july23_with_homing_small_horn_1deg_2';...
    'isaac_cal_jig_july23_with_homing_small_horn_1p5deg_2';...
    'isaac_cal_jig_july23_with_homing_small_horn_2deg_2';...
    'isaac_cal_jig_july23_with_homing_small_horn_0deg_3';...
    'isaac_cal_jig_july23_with_homing_small_flipped_0deg_1';...
    'isaac_cal_jig_july23_with_homing_small_flipped_0p5deg_1';...
    'isaac_cal_jig_july23_with_homing_small_flipped_1deg_1';...
    'isaac_cal_jig_july23_with_homing_small_flipped_1p5deg_1';...
    'isaac_cal_jig_july23_with_homing_small_flipped_2deg_1';...
    'isaac_cal_jig_july23_with_homing_small_flipped_m2deg_1';...
    'isaac_cal_jig_july23_with_homing_small_flipped_m1p5deg_1';...
    'isaac_cal_jig_july23_with_homing_small_flipped_m1deg_1';...
    'isaac_cal_jig_july23_with_homing_small_flipped_m0p5deg_1';...
    'isaac_cal_jig_july23_with_homing_small_flipped_0deg_2';...
    'isaac_cal_jig_july23_with_homing_small_flipped_0p5deg_2';...
    'isaac_cal_jig_july23_with_homing_small_flipped_1deg_2';...
    'isaac_cal_jig_july23_with_homing_small_flipped_1p5deg_2';...
    'isaac_cal_jig_july23_with_homing_small_flipped_2deg_2';...
    'isaac_cal_jig_july23_with_homing_small_flipped_0deg_3';...
    'isaac_cal_jig_july23_with_homing_small_flipped_isaac_0deg_1';...
    'isaac_cal_jig_july23_with_homing_small_flipped_isaac_0p5deg_1';...
    'isaac_cal_jig_july23_with_homing_small_flipped_isaac_1deg_1';...
    'isaac_cal_jig_july23_with_homing_small_flipped_isaac_1p5deg_1';...
    'isaac_cal_jig_july23_with_homing_small_flipped_isaac_2deg_1';...
    'isaac_cal_jig_july23_with_homing_small_flipped_isaac_m2deg_1';...
    'isaac_cal_jig_july23_with_homing_small_flipped_isaac_m1p5deg_1';...
    'isaac_cal_jig_july23_with_homing_small_flipped_isaac_m1deg_1';...
    'isaac_cal_jig_july23_with_homing_small_flipped_isaac_m0p5deg_1';...
    'isaac_cal_jig_july23_with_homing_small_flipped_isaac_0deg_2';...
    'isaac_cal_jig_july23_with_homing_small_flipped_isaac_0p5deg_2';...
    'isaac_cal_jig_july23_with_homing_small_flipped_isaac_1deg_2';...
    'isaac_cal_jig_july23_with_homing_small_flipped_isaac_1p5deg_2';...
    'isaac_cal_jig_july23_with_homing_small_flipped_isaac_2deg_2';...
    'isaac_cal_jig_july23_with_homing_small_flipped_isaac_0deg_3';...
    'isaac_cal_jig_july23_with_homing_small_flipped_isaac_m0p5deg_2';...
    'isaac_cal_jig_july23_with_homing_small_flipped_isaac_m1deg_2';...
    'isaac_cal_jig_july23_with_homing_small_flipped_isaac_m1p5deg_2';...
    'isaac_cal_jig_july23_with_homing_small_flipped_isaac_m2deg_2';...
    'isaac_cal_jig_july23_with_homing_small_flipped_isaac_0deg_4';...
    'isaac_cal_jig_july23_with_homing_translate_isaac_0in_1';...
    'isaac_cal_jig_july23_with_homing_translate_isaac_0p5in_1';...
    'isaac_cal_jig_july23_with_homing_translate_isaac_1in_1';...
    'isaac_cal_jig_july23_with_homing_translate_isaac_1p5in_1';...
    'isaac_cal_jig_july23_with_homing_translate_isaac_1in_2';...
    'isaac_cal_jig_july23_with_homing_translate_isaac_0p5in_2';...
    'isaac_cal_jig_july23_with_homing_translate_isaac_0in_2';...
    'isaac_cal_jig_july23_with_homing_translate_isaac_m0p5in_1';...
    'isaac_cal_jig_july23_with_homing_translate_isaac_m1in_1';...
    'isaac_cal_jig_july23_with_homing_translate_isaac_m1p5in_1';...
    'isaac_cal_jig_july23_with_homing_translate_isaac_0in_3';...
    'isaac_cal_jig_july23_with_homing_no_shroud_0deg_1';...
    'isaac_cal_jig_july23_with_homing_no_shroud_0p5deg_1';...
    'isaac_cal_jig_july23_with_homing_no_shroud_1deg_1';...
    'isaac_cal_jig_july23_with_homing_no_shroud_1p5deg_1';...
    'isaac_cal_jig_july23_with_homing_no_shroud_2deg_1';...
    'isaac_cal_jig_july23_with_homing_no_shroud_m2deg_1';...
    'isaac_cal_jig_july23_with_homing_no_shroud_m1p5deg_1';...
    'isaac_cal_jig_july23_with_homing_no_shroud_m1deg_1';...
    'isaac_cal_jig_july23_with_homing_no_shroud_m0p5deg_1';...
    'isaac_cal_jig_july23_with_homing_no_shroud_0deg_2';...
    'isaac_cal_jig_july23_with_homing_no_shroud_0p5deg_2';...
    'isaac_cal_jig_july23_with_homing_no_shroud_1deg_2';...
    'isaac_cal_jig_july23_with_homing_no_shroud_1p5deg_2';...
    'isaac_cal_jig_july23_with_homing_no_shroud_2deg_2';...
    'isaac_cal_jig_july23_with_homing_no_shroud_0deg_3';...
    'isaac_cal_jig_july23_with_homing_no_shroud_pt2_0deg_1';...
    'isaac_cal_jig_july23_with_homing_no_shroud_pt2_0p5deg_1';...
    'isaac_cal_jig_july23_with_homing_no_shroud_pt2_1deg_1';...
    'isaac_cal_jig_july23_with_homing_no_shroud_pt2_1p5deg_1';...
    'isaac_cal_jig_july23_with_homing_no_shroud_pt2_2deg_1';...
    'isaac_cal_jig_july23_with_homing_no_shroud_pt2_m2deg_1';...
    'isaac_cal_jig_july23_with_homing_no_shroud_pt2_m1p5deg_1';...
    'isaac_cal_jig_july23_with_homing_no_shroud_pt2_m1deg_1';...
    'isaac_cal_jig_july23_with_homing_no_shroud_pt2_m0p5deg_1';...
    'isaac_cal_jig_july23_with_homing_no_shroud_pt2_0deg_2';...
    'isaac_cal_jig_july23_with_homing_no_shroud_pt2_0p5deg_2';...
    'isaac_cal_jig_july23_with_homing_no_shroud_pt2_1deg_2';...
    'isaac_cal_jig_july23_with_homing_no_shroud_pt2_1p5deg_2';...
    'isaac_cal_jig_july23_with_homing_no_shroud_pt2_2deg_2';...
    'isaac_cal_jig_july23_with_homing_no_shroud_pt2_0deg_3';...
    'isaac_cal_jig_july23_with_homing_1p3meters_0deg_1';...
    'isaac_cal_jig_july23_with_homing_1p3meters_0p5deg_1';...
    'isaac_cal_jig_july23_with_homing_1p3meters_1deg_1';...
    'isaac_cal_jig_july23_with_homing_1p3meters_1p5deg_1';...
    'isaac_cal_jig_july23_with_homing_1p3meters_2deg_1';...
    'isaac_cal_jig_july23_with_homing_1p3meters_m2deg_1';...
    'isaac_cal_jig_july23_with_homing_1p3meters_m1p5deg_1';...
    'isaac_cal_jig_july23_with_homing_1p3meters_m1deg_1';...
    'isaac_cal_jig_july23_with_homing_1p3meters_m0p5deg_1';...
    'isaac_cal_jig_july23_with_homing_1p3meters_0deg_2';...
    'isaac_cal_jig_july23_with_homing_1p3meters_0p5deg_2';...
    'isaac_cal_jig_july23_with_homing_1p3meters_1deg_2';...
    'isaac_cal_jig_july23_with_homing_1p3meters_1p5deg_2';...
    'isaac_cal_jig_july23_with_homing_1p3meters_2deg_2';...
    'isaac_cal_jig_july23_with_homing_1p3meters_0deg_3';...
    'isaac_cal_jig_july23_with_homing_isaac_roll90_0deg_1';...
    'isaac_cal_jig_july23_with_homing_isaac_roll90_0p5deg_1';...
    'isaac_cal_jig_july23_with_homing_isaac_roll90_1deg_1';...
    'isaac_cal_jig_july23_with_homing_isaac_roll90_1p5deg_1';...
    'isaac_cal_jig_july23_with_homing_isaac_roll90_2deg_1';...
    'isaac_cal_jig_july23_with_homing_isaac_roll90_m2deg_1';...
    'isaac_cal_jig_july23_with_homing_isaac_roll90_m1p5deg_1';...
    'isaac_cal_jig_july23_with_homing_isaac_roll90_m1deg_1';...
    'isaac_cal_jig_july23_with_homing_isaac_roll90_m0p5deg_1';...
    'isaac_cal_jig_july23_with_homing_isaac_roll90_0deg_2';...
    'isaac_cal_jig_july23_with_homing_isaac_roll90_0p5deg_2';...
    'isaac_cal_jig_july23_with_homing_isaac_roll90_1deg_2';...
    'isaac_cal_jig_july23_with_homing_isaac_roll90_1p5deg_2';...
    'isaac_cal_jig_july23_with_homing_isaac_roll90_2deg_2';...
    'isaac_cal_jig_july23_with_homing_isaac_roll90_0deg_3';...
    'isaac_cal_jig_july23_with_homing_2p4meters_0deg_1';...
    'isaac_cal_jig_july23_with_homing_2p4meters_0p5deg_1';...
    'isaac_cal_jig_july23_with_homing_2p4meters_1deg_1';...
    'isaac_cal_jig_july23_with_homing_2p4meters_1p5deg_1';...
    'isaac_cal_jig_july23_with_homing_2p4meters_2deg_1';...
    'isaac_cal_jig_july23_with_homing_2p4meters_m2deg_1';...
    'isaac_cal_jig_july23_with_homing_2p4meters_m1p5deg_1';...
    'isaac_cal_jig_july23_with_homing_2p4meters_m1deg_1';...
    'isaac_cal_jig_july23_with_homing_2p4meters_m0p5deg_1';...
    'isaac_cal_jig_july23_with_homing_2p4meters_0deg_2';...
    'isaac_cal_jig_july23_with_homing_2p4meters_0p5deg_2';...
    'isaac_cal_jig_july23_with_homing_2p4meters_1deg_2';...
    'isaac_cal_jig_july23_with_homing_2p4meters_1p5deg_2';...
    'isaac_cal_jig_july23_with_homing_2p4meters_2deg_2';...
    'isaac_cal_jig_july23_with_homing_2p4meters_0deg_3';...
    };

labs = {...;
    '0dB Attenuation, With Homing';...
    '0dB Attenuation, No Homing';...
    '30dB Attenuation, No Homing';...
    '45dB Attenuation, No Homing';...
    '55dB Attenuation, No Homing';...
    'Indoor test with grids, with homing';...
    'Indoor test with grids, No homing';...
    'Outdoor test, No homing';... 8
    'Outdoor test, Same platform, No Homing';...
    'Outdoor test, Same platform, Releveled';... 10
    'Outdoor test, After ISAAC grid check';...
    'Outdoor test, Homed RPS once';... 12
    'Outdoor test, Homed RPS every scan';...
    'Outdoor test, Homed RPS every scan';... 14
    'Indoor Test, Homed RPS every scan';...
    'Indoor Test, Homed RPS every scan, with stage heater';...16
    'Indoor Test, Homed RPS every scan, with stage heater';...
    'Indoor Test, Homed RPS every scan, with stage heater';...%18
    'Indoor Test, Added Fine Homing Switch';...
    'Indoor Test, Added Fine Homing Switch, Wrong Homing Position';...
    'Indoor Test, Added Fine Homing Switch';...%21
    'Outdoor Test, Pre RPS Observation, Homed Every Scan';...
    'Outdoor Test, Post RPS Observation, Homed Every Scan';...
    'Indoor Test, Post RPS Observation';... % 24
    'Indoor Test, Post RPS Observation, ISAAC Refrozen';...
    'Indoor Test, Post RPS Observation, RPS Refrozen';...
    'On Alignment Jig, No Homing, Dist = 20"';... % 27
    'On Alignment Jig, No Homing, Dist = 40"';...
    'On Alignment Jig, No Homing, Dist = 40"';...
    'On Alignment Jig, No Homing, Dist = 40", Tilt -0.2';... %30
    'On Alignment Jig, No Homing, Dist = 40", Tilt -0.4';... %31
    'On Alignment Jig, No Homing, Dist = 20"';... % 32
    'On Alignment Jig, No Homing, Dist = 20", repeat';...
    'On Alignment Jig, Homing, Dist = 20"';...
    'On Alignment Jig, Homing, Dist = 20"';... % 35
    'On Alignment Jig, Homing, Dist = 20"';... % 36
    'On Alignment Jig, Homing, Dist = 20", offset 1"';... % 37
    'On Alignment Jig, Homing, Dist = 37"';... % 38
    'On Alignment Jig, Homing, Dist = 37", repeat';... % 39
    'On Alignment Jig, Homing, Dist = 37", zeroed RPS';... % 40
    'On Alignment Jig, Homing, Dist = 37", zeroed isaac';... % 41
    'On Alignment Jig, Homing, Dist = 37", added baffling';... % 42
    'On Alignment Jig, Homing, Dist = 40", offsets';... % 43
    'On Alignment Jig, No Homing, Dist = 26", No LNA';... % 44
    'On Alignment Jig, No Homing, Dist = 26", No LNA, 5.45 homing';... % 45
    'On Alignment Jig, Homing, Dist = 26"';... % 46
    'On Alignment Jig, Homing, Dist = 26", align-check +0deg';... % 47
    'On Alignment Jig, Homing, Dist = 26", align-check +1deg';... % 48
    'On Alignment Jig, Homing, Dist = 26", align-check +2deg';... % 49
    'On Alignment Jig, Homing, Dist = 26", align-check +0deg';... % 50
    'On Alignment Jig, Homing, Dist = 26", align-check +1deg';... % 51
    'On Alignment Jig, Homing, Dist = 26", align-check +2deg';... % 52
    'On Alignment Jig, Homing, Dist = 26", align-check -1deg';... % 53
    'On Alignment Jig, Homing, Dist = 26", align-check -2deg';... % 54
    'On Alignment Jig, Homing, Dist = 26", align-check -1deg';... % 55
    'On Alignment Jig, Homing, Dist = 26", align-check +0deg';... % 56
    'On Alignment Jig, Homing, Dist = 26", align-check +0.5deg';... % 57
    'On Alignment Jig, Homing, Dist = 26", align-check +0.25deg';... % 58
    'On Alignment Jig, Homing, Dist = 26", align-check -0.25deg';... % 59
    'On Alignment Jig, Homing, Dist = 26", align-check -0.5deg';... % 60
    'On Alignment Jig, Homing, Dist = 26", align-check -0.125deg';... % 61
    'On Alignment Jig, Homing, Dist = 26", Noise Check'; ... % 62
    'On Alignment Jig, Homing, Dist = 26", cont-amp +0deg';... % 63
    'On Alignment Jig, Homing, Dist = 26", cont-amp +1deg';... % 64
    'On Alignment Jig, Homing, Dist = 26", cont-amp +2deg';... % 65
    'On Alignment Jig, Homing, Dist = 26", cont-amp +1deg';... % 66
    'On Alignment Jig, Homing, Dist = 26", cont-amp +0deg';... % 67
    'On Alignment Jig, Homing, Dist = 26", cont-amp -1deg';... % 68
    'On Alignment Jig, Homing, Dist = 26", cont-amp -2deg';... % 69
    'On Alignment Jig, Homing, Dist = 43", +0deg';... % 70
    'On Alignment Jig, Homing, Dist = 43", +0.5deg';... % 71
    'On Alignment Jig, Homing, Dist = 43", +1deg';... % 72
    'On Alignment Jig, Homing, Dist = 43", +1.5deg';... % 73
    'On Alignment Jig, Homing, Dist = 43", +2deg';... % 74
    'On Alignment Jig, Homing, Dist = 43", +1deg';... % 75
    'On Alignment Jig, Homing, Dist = 43", +0deg';... % 76
    'On Alignment Jig, Homing, Dist = 43", -0.5deg';... % 77
    'On Alignment Jig, Homing, Dist = 43", -1deg';... % 78
    'On Alignment Jig, Homing, Dist = 43", -2deg';... % 79
    'On Alignment Jig, Homing, Dist = 43", 0deg ISAAC Offset';... % 80
    'On Alignment Jig, Homing, Dist = 43", 0.5deg ISAAC Offset';... % 81
    'On Alignment Jig, Homing, Dist = 43", 1deg ISAAC Offset';... % 82
    'On Alignment Jig, Homing, Dist = 43", 1.5deg ISAAC Offset';... % 83
    'On Alignment Jig, Homing, Dist = 43", 2deg ISAAC Offset';... % 84
    'On Alignment Jig, Homing, Dist = 43", -2deg ISAAC Offset';... % 85
    'On Alignment Jig, Homing, Dist = 43", -1.5deg ISAAC Offset';... % 86
    'On Alignment Jig, Homing, Dist = 43", -1deg ISAAC Offset';... % 87
    'On Alignment Jig, Homing, Dist = 43", -0.5deg ISAAC Offset';... % 88
    'On Alignment Jig, Homing, Dist = 43", 0deg ISAAC Offset';... % 89
    'On Alignment Jig, Homing, Dist = 43", Absorber test: no absorber';... 90
    'On Alignment Jig, Homing, Dist = 43", Absorber test: behind isaac';... 91
    'On Alignment Jig, Homing, Dist = 43", Absorber test, in window';... 92
    'On Alignment Jig, Homing, Dist = 43", Absorber test, on table';... 93
    'On Alignment Jig, Homing, Dist = 43", Absorber test, table rps up';... 94
    'On Alignment Jig, Homing, Dist = 43", Absorber test, table isaac up';... 95
    'On Alignment Jig, Homing, Dist = 43", Absorber test, fully blocked';... 96
    'On Alignment Jig, Homing, Dist = 43", Absorber test, both sides';... 97
    'On Jig, Homing, Dist = 43", Horn Flipped 180, RPS align +0.0deg';... 98
    'On Jig, Homing, Dist = 43", Horn Flipped 180, RPS align +0.5deg';... 99
    'On Jig, Homing, Dist = 43", Horn Flipped 180, RPS align +1.0deg';... 100
    'On Jig, Homing, Dist = 43", Horn Flipped 180, RPS align +1.5deg';... 101
    'On Jig, Homing, Dist = 43", Horn Flipped 180, RPS align +2.0deg';... 102
    'On Jig, Homing, Dist = 43", Horn Flipped 180, RPS align -2.0deg';... 103
    'On Jig, Homing, Dist = 43", Horn Flipped 180, RPS align -1.5deg';... 104
    'On Jig, Homing, Dist = 43", Horn Flipped 180, RPS align -1.0deg';... 105
    'On Jig, Homing, Dist = 43", Horn Flipped 180, RPS align -0.5deg';... 106
    'On Jig, Homing, Dist = 43", Horn Flipped 180, RPS align +0.0deg';... 107
    'On Jig, Homing, Dist = 43", Small Horn, RPS align +0.0deg';... 108
    'On Jig, Homing, Dist = 43", Small Horn, RPS align +0.5deg';... 109
    'On Jig, Homing, Dist = 43", Small Horn, RPS align +1.0deg';... 110
    'On Jig, Homing, Dist = 43", Small Horn, RPS align +1.5deg';... 111
    'On Jig, Homing, Dist = 43", Small Horn, RPS align +2.0deg';... 112
    'On Jig, Homing, Dist = 43", Small Horn, RPS align -2.0deg';... 113
    'On Jig, Homing, Dist = 43", Small Horn, RPS align -1.5deg';... 114
    'On Jig, Homing, Dist = 43", Small Horn, RPS align -1.0deg';... 115
    'On Jig, Homing, Dist = 43", Small Horn, RPS align -0.5deg';... 116
    'On Jig, Homing, Dist = 43", Small Horn, RPS align +0.0deg';... 117
    'On Jig, Homing, Dist = 43", Small Horn, RPS align +0.5deg';... 118
    'On Jig, Homing, Dist = 43", Small Horn, RPS align +1.0deg';... 119
    'On Jig, Homing, Dist = 43", Small Horn, RPS align +1.5deg';... 120
    'On Jig, Homing, Dist = 43", Small Horn, RPS align +2.0deg';... 121
    'On Jig, Homing, Dist = 43", Small Horn, RPS align +0.0deg';... 122
    'On Jig, Homing, Dist = 43", Small Horn Flipped, RPS align +0.0deg';... 123
    'On Jig, Homing, Dist = 43", Small Horn Flipped, RPS align +0.5deg';... 124
    'On Jig, Homing, Dist = 43", Small Horn Flipped, RPS align +1.0deg';... 125
    'On Jig, Homing, Dist = 43", Small Horn Flipped, RPS align +1.5deg';... 126
    'On Jig, Homing, Dist = 43", Small Horn Flipped, RPS align +2.0deg';... 127
    'On Jig, Homing, Dist = 43", Small Horn Flipped, RPS align -2.0deg';... 128
    'On Jig, Homing, Dist = 43", Small Horn Flipped, RPS align -1.5deg';... 129
    'On Jig, Homing, Dist = 43", Small Horn Flipped, RPS align -1.0deg';... 130
    'On Jig, Homing, Dist = 43", Small Horn Flipped, RPS align -0.5deg';... 131
    'On Jig, Homing, Dist = 43", Small Horn Flipped, RPS align +0.0deg';... 132
    'On Jig, Homing, Dist = 43", Small Horn Flipped, RPS align +0.5deg';... 133
    'On Jig, Homing, Dist = 43", Small Horn Flipped, RPS align +1.0deg';... 134
    'On Jig, Homing, Dist = 43", Small Horn Flipped, RPS align +1.5deg';... 135
    'On Jig, Homing, Dist = 43", Small Horn Flipped, RPS align +2.0deg';... 136
    'On Jig, Homing, Dist = 43", Small Horn Flipped, RPS align +0.0deg';... 137
    'On Jig, Homing, Dist = 43", Small Horn Flipped, ISAAC align +0.0deg';... 138
    'On Jig, Homing, Dist = 43", Small Horn Flipped, ISAAC align +0.5deg';... 139
    'On Jig, Homing, Dist = 43", Small Horn Flipped, ISAAC align +1.0deg';... 140
    'On Jig, Homing, Dist = 43", Small Horn Flipped, ISAAC align +1.5deg';... 141
    'On Jig, Homing, Dist = 43", Small Horn Flipped, ISAAC align +2.0deg';... 142
    'On Jig, Homing, Dist = 43", Small Horn Flipped, ISAAC align -2.0deg';... 143
    'On Jig, Homing, Dist = 43", Small Horn Flipped, ISAAC align -1.5deg';... 144
    'On Jig, Homing, Dist = 43", Small Horn Flipped, ISAAC align -1.0deg';... 145
    'On Jig, Homing, Dist = 43", Small Horn Flipped, ISAAC align -0.5deg';... 146
    'On Jig, Homing, Dist = 43", Small Horn Flipped, ISAAC align +0.0deg';... 147
    'On Jig, Homing, Dist = 43", Small Horn Flipped, ISAAC align +0.5deg';... 148
    'On Jig, Homing, Dist = 43", Small Horn Flipped, ISAAC align +1.0deg';... 149
    'On Jig, Homing, Dist = 43", Small Horn Flipped, ISAAC align +1.5deg';... 150
    'On Jig, Homing, Dist = 43", Small Horn Flipped, ISAAC align +2.0deg';... 151
    'On Jig, Homing, Dist = 43", Small Horn Flipped, ISAAC align +0.0deg';... 152
    'On Jig, Homing, Dist = 43", Small Horn Flipped, ISAAC align -0.5deg';... 153
    'On Jig, Homing, Dist = 43", Small Horn Flipped, ISAAC align -1.0deg';... 154
    'On Jig, Homing, Dist = 43", Small Horn Flipped, ISAAC align -1.5deg';... 155
    'On Jig, Homing, Dist = 43", Small Horn Flipped, ISAAC align -2.0deg';... 156
    'On Jig, Homing, Dist = 43", Small Horn Flipped, ISAAC align +0.0deg';... 157
    'On Jig, Homing, Dist = 43", Small Horn Flipped, Moved ISAAC +0.0in';... 158
    'On Jig, Homing, Dist = 43", Small Horn Flipped, Moved ISAAC +0.5in';... 159
    'On Jig, Homing, Dist = 43", Small Horn Flipped, Moved ISAAC +1.0in';... 160
    'On Jig, Homing, Dist = 43", Small Horn Flipped, Moved ISAAC +1.5in';... 161
    'On Jig, Homing, Dist = 43", Small Horn Flipped, Moved ISAAC +1.0in';... 162
    'On Jig, Homing, Dist = 43", Small Horn Flipped, Moved ISAAC +0.5in';... 163
    'On Jig, Homing, Dist = 43", Small Horn Flipped, Moved ISAAC +0.0in';... 164
    'On Jig, Homing, Dist = 43", Small Horn Flipped, Moved ISAAC -0.5in';... 165
    'On Jig, Homing, Dist = 43", Small Horn Flipped, Moved ISAAC -1.0in';... 166
    'On Jig, Homing, Dist = 43", Small Horn Flipped, Moved ISAAC -1.5in';... 167
    'On Jig, Homing, Dist = 43", Small Horn Flipped, Moved ISAAC +0.0in';... 168
    'On Jig, Homing, Dist = 43", RPS Shroud Removed, RPS align +0.0deg';... 169
    'On Jig, Homing, Dist = 43", RPS Shroud Removed, RPS align +0.5deg';... 170
    'On Jig, Homing, Dist = 43", RPS Shroud Removed, RPS align +1.0deg';... 171
    'On Jig, Homing, Dist = 43", RPS Shroud Removed, RPS align +1.5deg';... 172
    'On Jig, Homing, Dist = 43", RPS Shroud Removed, RPS align +2.0deg';... 173
    'On Jig, Homing, Dist = 43", RPS Shroud Removed, RPS align -2.0deg';... 174
    'On Jig, Homing, Dist = 43", RPS Shroud Removed, RPS align -1.5deg';... 175
    'On Jig, Homing, Dist = 43", RPS Shroud Removed, RPS align -1.0deg';... 176
    'On Jig, Homing, Dist = 43", RPS Shroud Removed, RPS align -0.5deg';... 177
    'On Jig, Homing, Dist = 43", RPS Shroud Removed, RPS align +0.0deg';... 178
    'On Jig, Homing, Dist = 43", RPS Shroud Removed, RPS align +0.5deg';... 179
    'On Jig, Homing, Dist = 43", RPS Shroud Removed, RPS align +1.0deg';... 180
    'On Jig, Homing, Dist = 43", RPS Shroud Removed, RPS align +1.5deg';... 181
    'On Jig, Homing, Dist = 43", RPS Shroud Removed, RPS align +2.0deg';... 182
    'On Jig, Homing, Dist = 43", RPS Shroud Removed, RPS align +0.0deg';... 183 -------
    'On Jig, Homing, Dist = 43", RPS Shroud Removed, round 2, RPS align +0.0deg';... 184
    'On Jig, Homing, Dist = 43", RPS Shroud Removed, round 2, RPS align +0.5deg';... 185
    'On Jig, Homing, Dist = 43", RPS Shroud Removed, round 2, RPS align +1.0deg';... 186
    'On Jig, Homing, Dist = 43", RPS Shroud Removed, round 2, RPS align +1.5deg';... 187
    'On Jig, Homing, Dist = 43", RPS Shroud Removed, round 2, RPS align +2.0deg';... 188
    'On Jig, Homing, Dist = 43", RPS Shroud Removed, round 2, RPS align -2.0deg';... 189
    'On Jig, Homing, Dist = 43", RPS Shroud Removed, round 2, RPS align -1.5deg';... 190
    'On Jig, Homing, Dist = 43", RPS Shroud Removed, round 2, RPS align -1.0deg';... 191
    'On Jig, Homing, Dist = 43", RPS Shroud Removed, round 2, RPS align -0.5deg';... 192
    'On Jig, Homing, Dist = 43", RPS Shroud Removed, round 2, RPS align +0.0deg';... 193
    'On Jig, Homing, Dist = 43", RPS Shroud Removed, round 2, RPS align +0.5deg';... 194
    'On Jig, Homing, Dist = 43", RPS Shroud Removed, round 2, RPS align +1.0deg';... 195
    'On Jig, Homing, Dist = 43", RPS Shroud Removed, round 2, RPS align +1.5deg';... 196
    'On Jig, Homing, Dist = 43", RPS Shroud Removed, round 2, RPS align +2.0deg';... 197
    'On Jig, Homing, Dist = 43", RPS Shroud Removed, round 2, RPS align +0.0deg';... 198
    'On Jig, Homing, Dist = 53", Small Horn Flipped, RPS align +0.0deg';... 199
    'On Jig, Homing, Dist = 53", Small Horn Flipped, RPS align +0.5deg';... 200
    'On Jig, Homing, Dist = 53", Small Horn Flipped, RPS align +1.0deg';... 201
    'On Jig, Homing, Dist = 53", Small Horn Flipped, RPS align +1.5deg';... 202
    'On Jig, Homing, Dist = 53", Small Horn Flipped, RPS align +2.0deg';... 203
    'On Jig, Homing, Dist = 53", Small Horn Flipped, RPS align -2.0deg';... 204
    'On Jig, Homing, Dist = 53", Small Horn Flipped, RPS align -1.5deg';... 205
    'On Jig, Homing, Dist = 53", Small Horn Flipped, RPS align -1.0deg';... 206
    'On Jig, Homing, Dist = 53", Small Horn Flipped, RPS align -0.5deg';... 207
    'On Jig, Homing, Dist = 53", Small Horn Flipped, RPS align +0.0deg';... 208
    'On Jig, Homing, Dist = 53", Small Horn Flipped, RPS align +0.5deg';... 209
    'On Jig, Homing, Dist = 53", Small Horn Flipped, RPS align +1.0deg';... 210
    'On Jig, Homing, Dist = 53", Small Horn Flipped, RPS align +1.5deg';... 211
    'On Jig, Homing, Dist = 53", Small Horn Flipped, RPS align +2.0deg';... 212
    'On Jig, Homing, Dist = 53", Small Horn Flipped, RPS align +0.0deg';... 213
    'On Jig, Homing, Dist = 43", ISAAC roll 90deg, RPS align +0.0deg';... 214
    'On Jig, Homing, Dist = 43", ISAAC roll 90deg, RPS align +0.5deg';... 215
    'On Jig, Homing, Dist = 43", ISAAC roll 90deg, RPS align +1.0deg';... 216
    'On Jig, Homing, Dist = 43", ISAAC roll 90deg, RPS align +1.5deg';... 217
    'On Jig, Homing, Dist = 43", ISAAC roll 90deg, RPS align +2.0deg';... 218
    'On Jig, Homing, Dist = 43", ISAAC roll 90deg, RPS align -2.0deg';... 219
    'On Jig, Homing, Dist = 43", ISAAC roll 90deg, RPS align -1.5deg';... 220
    'On Jig, Homing, Dist = 43", ISAAC roll 90deg, RPS align -1.0deg';... 221
    'On Jig, Homing, Dist = 43", ISAAC roll 90deg, RPS align -0.5deg';... 222
    'On Jig, Homing, Dist = 43", ISAAC roll 90deg, RPS align +0.0deg';... 223
    'On Jig, Homing, Dist = 43", ISAAC roll 90deg, RPS align +0.5deg';... 224
    'On Jig, Homing, Dist = 43", ISAAC roll 90deg, RPS align +1.0deg';... 225
    'On Jig, Homing, Dist = 43", ISAAC roll 90deg, RPS align +1.5deg';... 226
    'On Jig, Homing, Dist = 43", ISAAC roll 90deg, RPS align +2.0deg';... 227
    'On Jig, Homing, Dist = 43", ISAAC roll 90deg, RPS align +0.0deg';... 228 
    'On Jig, Homing, Dist = 93", Small Horn Flipped, RPS align +0.0deg';... 229
    'On Jig, Homing, Dist = 93", Small Horn Flipped, RPS align +0.5deg';... 220
    'On Jig, Homing, Dist = 93", Small Horn Flipped, RPS align +1.0deg';... 231
    'On Jig, Homing, Dist = 93", Small Horn Flipped, RPS align +1.5deg';... 232
    'On Jig, Homing, Dist = 93", Small Horn Flipped, RPS align +2.0deg';... 233
    'On Jig, Homing, Dist = 93", Small Horn Flipped, RPS align -2.0deg';... 234
    'On Jig, Homing, Dist = 93", Small Horn Flipped, RPS align -1.5deg';... 235
    'On Jig, Homing, Dist = 93", Small Horn Flipped, RPS align -1.0deg';... 236
    'On Jig, Homing, Dist = 93", Small Horn Flipped, RPS align -0.5deg';... 237
    'On Jig, Homing, Dist = 93", Small Horn Flipped, RPS align +0.0deg';... 238
    'On Jig, Homing, Dist = 93", Small Horn Flipped, RPS align +0.5deg';... 239
    'On Jig, Homing, Dist = 93", Small Horn Flipped, RPS align +1.0deg';... 240
    'On Jig, Homing, Dist = 93", Small Horn Flipped, RPS align +1.5deg';... 241
    'On Jig, Homing, Dist = 93", Small Horn Flipped, RPS align +2.0deg';... 242
    'On Jig, Homing, Dist = 93", Small Horn Flipped, RPS align +0.0deg';... 243
    };

%end