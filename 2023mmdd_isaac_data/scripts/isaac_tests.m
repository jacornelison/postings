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


%%
clc
addpath('z:/pipeline/beammap/')
gitdir = fullfile('C:','Users','James','Documents','');
%figdir = fullfile(gitdir,'GitHub','postings','2021mmdd_isaac_test','figs');
figdir = fullfile(gitdir,'GitHub','postings','2023mmdd_isaac_data','figs');
direct = fullfile(gitdir,'GitHub','rps_cal','data','');
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
    'isaac_cal_8_with_home_indoor';... % Added a 'fine' homing switch to the stepper motor.
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
    'isaac_cal_jig_withhome_4';... % 29 Jul 2022 @ 37"; homed.
    'isaac_cal_jig_withhome_5';... % 31 Jul 2022 @ 37"; homed. Repeat of last
    'isaac_cal_jig_withhome_6';... % 31 Jul 2022 @ 37"; homed. Zeroed RPS to -0.03
    'isaac_cal_jig_withhome_7';... % 31 Jul 2022 @ 37"; homed. Zeroed ISAAC
    'isaac_cal_jig_withhome_8';... % 31 Jul 2022 @ 37"; homed. Added baffling.
    };

labs = {...;
    '0dB Attenuation, With Homing';...
    '0dB Attenuation, No Homing';...
    '30dB Attenuation, No Homing';...
    '45dB Attenuation, No Homing';...
    '55dB Attenuation, No Homing';...
    'Indoor test with grids, with homing';...
    'Indoor test with grids, No homing';...
    'Outdoor test, No homing';...
    'Outdoor test, Same platform, No Homing';...
    'Outdoor test, Same platform, Releveled';...
    'Outdoor test, After ISAAC grid check';...
    'Outdoor test, Homed RPS once';...
    'Outdoor test, Homed RPS every scan';...
    'Outdoor test, Homed RPS every scan';...
    'Indoor Test, Homed RPS every scan';...
    'Indoor Test, Homed RPS every scan, with stage heater';...
    'Indoor Test, Homed RPS every scan, with stage heater';...
    'Indoor Test, Homed RPS every scan, with stage heater';...
    'Indoor Test, Added Fine Homing Switch';...
    'Indoor Test, Added Fine Homing Switch, Wrong Homing Position';...
    'Indoor Test, Added Fine Homing Switch';...
    'Outdoor Test, With Fine Homing Sw, Homed Every Scan';...
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
    'On Alignment Jig, Homing, Dist = 37"';... % 37
    'On Alignment Jig, Homing, Dist = 37", repeat';... % 38
    'On Alignment Jig, Homing, Dist = 37", zeroed RPS';... % 39
    'On Alignment Jig, Homing, Dist = 37", zeroed isaac';... % 40
    'On Alignment Jig, Homing, Dist = 37", added baffling';... % 41
    };
if ~exist('meanguess','var')
    meanguess = 0;
end

fd = struct();
flds = {'schnum','rownum','phi','eps','N1','N2','A','time','tilt','istilt','temp','istemp','istilttemp','wgt','fitind','phi_isaac'};
vals = {'prefind','di','parm(1)','parm(2)','parm(3)','parm(4)','parm(5)','nanmean(TIME)',...
    'nanmean(T)','nanmean(T2)','nanmean(TEMP)','nanmean(ISTEMP)','nanmean(ISTILTTEMP)','wgtind','fitind','nanmean(a_isaac)'};
for fldind = 1:length(flds)
    fd.(flds{fldind}) = [];
end

reses = {};
angs = {};

thresh = 0.7;
weighttitles = {'Uniform',sprintf('Amp<%0.1f',thresh),'1/Amp'};
weightnames = {'uniform','threshold','one_on_r'};

fitnames = {'fmin','lsq','complex'};
plotmodcurve = 0;

dists = [ones(1,27), 20*0.0254, ones(1,3)*40*0.0254,ones(1,5)*20*0.0254,ones(1,5)*37*0.0254];
for prefind = 21:41


    if ismember(prefind,[1:20 23:25])
        rpscal = rps_tilt_cals_all{end};
        isaaccal = isaac_tilt_cals_all{end};
    elseif ismember(prefind,[21:22])
        rpscal = rps_tilt_cals_all{end};
        isaaccal = isaac_tilt_cals_all{end};
    elseif ismember(prefind,[26:length(prefix)])
        rpscal = rps_tilt_cals_all{end};
        isaaccal = isaac_tilt_cals_all{end};
        %isaaccal = [0.28 -0.1778];
    end

    d = dir(fullfile(direct,[prefix{prefind} '_scan*']));
    disp(labs{prefind})
    data = {};
    for di = 1:length(d)
        if ~ismember(prefind, [8:21 22:26])
            %disp(d(di).name)
            data{di} = load_lj_data(fullfile(direct,d(di).name));
        else
            %disp([prefix{prefind} sprintf('_scan_%i_loop_1.csv',di-1)])
            data{di} = load_lj_data(fullfile(direct,[prefix{prefind} ...
                sprintf('_scan_%i_loop_1.csv',di-1)]));
        end
    end
    for fitind = 2%1:length(fitnames)
        %fittype=fitnames{fitind};
        fittype = 'lsq0';
        for wgtind = 1%1:3

            lock_cal = 50/10; %mV per V
            if strcmp(fittype,'complex')
                parms = NaN(length(data),7);
            else
            parms = NaN(length(data),5);
            end
            [times, chis, itempstd, itempmn] = deal(NaN(length(data)));
            [mod_curve, res] = deal(NaN(length(unique(data{1}.Angle)),length(data)));


            if plotmodcurve
                fig = figure(prefind);
                fig.Position(3:4) = [1000,300];
                clf; hold on;
            end
            for di = 1:length(data)
                a = unique(data{di}.Angle);
                [R, Rs, T, T2, TIME, TEMP, ISTILTTEMP, ISTEMP] = deal([]);
                for ai = 1:length(a)
                    aind = data{di}.Angle == a(ai);
                    R = [R; mean(sqrt(data{di}.X(aind).^2+data{di}.Y(aind).^2))];
                    Rs = [Rs; std(sqrt(data{di}.X(aind).^2+data{di}.Y(aind).^2))];
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

                end
                R = R./max(R);

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
                elseif ismember(prefind,[12,20,22:length(prefix)])
                    homeangle = 5.454545;
                else
                    homeangle = 0;
                end


                % We're measuring the angle of the ISAAC WRT gravity:
                % Clockwise looking at ISAAC is positive
                % RPS has opposite parity
                % Angle_out = RPS_angle - RPS tilt - RPS grid horz. + ISAAC Tilt + ISAAC grid horz.
                rps_tilt = T;
                rps_grid = homeangle;
                isaac_tilt = T2;
                isaac_grid = 0.17;
                a = a + rps_tilt - rps_grid;% + isaac_tilt + isaac_grid;
                a_isaac = mean(isaac_grid - isaac_tilt);
                mod_curve(:,di) = R;

                mxfev = 100000;
                mxiter = 100000;
                options = optimset('TolFun',1e-10,'MaxIter',mxiter,'MaxFunEvals',mxfev,'Display','off');
                lb = [-20 -0.5 -10 -10 0];
                ub = [20 0.5 10 10 10];
                meanguess = 0;%0;
                for nind = 1%1:3

                    guess = [meanguess,0.002,0,0,max(R)/2]+1e-6;
                    %fittype = 'lsq';
                    switch fittype
                        case 'basic'
                            modfunc = @(ang,e,A) A/2*(cosd(2*(a-ang))-(e+1)/(e-1))+(0*cosd(a-T)+0*sind(a-T));
                            chifun = @(x) sum((R-modfunc(x(1),x(2),x(3))).^2);
                            parm = fminsearch(chifun,[0.75,0,max(R)]);
                            %plot(a,R-modfunc(parm(1),parm(2),parm(3)))
                            res(:,di) = R-modfunc(parm(1),parm(2),parm(3));
                        case 'fmin'

                            modfunc = @(ang,e,A,N1,N2) A/2*(cosd(2*(a-ang))-(e+1)/(e-1)).*(N1*cosd(a)+N2*sind(a)+1);
                            chifun = @(x) sum((R-modfunc(x(1),x(2),x(3),x(4),x(5))).^2);

                            chifun = @(x) sum(w.*(R-rps_get_mod_model(x,a)).^2);

                            parm = fminsearch(chifun,guess,options);
                            res(:,di) = R-modfunc(parm(1),parm(2),parm(3),parm(4),parm(5));

                        case 'lsq0'
                            lb = [-20 -0.5 -10 -10 0];
                            ub = [20 0.5 10 10 1e6];
                            parm = lsqcurvefit(@rps_get_mod_model,guess,a,R,lb,ub,options);
                            res(:,di) = R-rps_get_mod_model(parm,a);
                        case 'lsq'

                            chifunc = @(p) w.*(R-rps_get_mod_model(p,a));
                            parm = lsqnonlin(chifunc,guess,lb,ub,options);
                            res(:,di) = R-rps_get_mod_model(parm,a);
                        case 'complex'
                            
                            
                            
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
                    meanguess = nanmean(parms(:,1));
                    if plotmodcurve
                        
                        fig = figure(prefind);
                        %plot(a,R)
                        ind1 = 1:7;
                        ind2 = 7:13;
                        plot(res(:,di))
                        ylim([-1 1]*0.01)
                        %plot(R(ind1)-R(ind2));
                        %ylim([-1 1]*0.04)
                        %xlim([-1 1]*190)
                        %ylim([-0.1,6])
                        xlabel('Angle (^o)')
                        ylabel('Lock-In X (Volts)')
                        grid on
                        title(labs{prefind})

                    end
                    %plot(a,res(:,di))

                end

                if abs(R(1)-R((length(R)+1)/2))<0.05
                    for fldind = 1:length(flds)
                        eval(sprintf('fd.%s(end+1) = %s;',flds{fldind},vals{fldind}))
                    end
                    angs{end+1} = a;
                    parm(1) = 0;
                    reses{end+1} = R-rps_get_mod_model(parm,a);
                else
                    parms(di,1:length(parm)) = NaN(1,length(parm));
                end

            end

%             % cuts
%             switch prefind
%                 case 19
%                     ind = [4:6 8:size(parms,1)];
%                 case 23
%                     ind = 2:size(parms,1);
%                 case 26
%                     ind = [1, 3:size(parms,1)];
%                 otherwise
%                     ind = 1:size(parms,1);
%             end
%             parms = parms(ind,:);

            if 0
                lims = {[-0.5 1.5] [0 1]*0.06};
                offs = {0 0};
                titles = {'\phi [Degrees]','\epsilon'};
                fig = figure(2);
                fig.Position(3:4) = [560*2 420];
                clf; hold on;

                for parmind = 1:2
                    subplot(1,2,parmind)
                    resolution = 20;
                    edges = lims{parmind}(1):diff(lims{parmind})/resolution:lims{parmind}(2)+offs{parmind};
                    N = histc(parms(:,parmind),edges);
                    bar(edges,N,'histc')
                    xlim(lims{parmind}+offs{parmind})
                    xlabel(titles{parmind});
                    ylabel('N')
                    title({labs{prefind},...
                        sprintf('Mean=%0.3f STD=%0.3f N=%02d, %s',nanmean(parms(:,parmind)),nanstd(parms(:,parmind)),length(data),weighttitles{wgtind})})
                    grid on
                end

                saveas(fig,fullfile(figdir,sprintf('isaac_cal_%03d_%s_%s.png',prefind,weightnames{wgtind},fittype)))
            end
        end
    end
    fprintf('Angle: %0.3f +/- %0.3f\n',nanmean(parms(:,1))-a_isaac,nanstd(parms(:,1)))
    fprintf('Eff: %0.4f +/- %0.4f\n',nanmean(parms(:,2)),nanstd(parms(:,2)))
end


%% Look the jig data

for wgtind = 1%1:3
    fig = figure(2+wgtind);
    clf; hold on;
    % ind = ismember(fd.schnum,[28:31 37:39]) & fd.wgt == wgtind;
    % plot(fd.schnum(ind),fd.phi(ind)+1*5.454545-1*fd.tilt(ind)-1*fd.istilt(ind),'.')
    % ind = ~ismember(fd.schnum,[28:31 37:39]) & fd.wgt == wgtind;
    % plot(fd.schnum(ind),fd.phi(ind)+1*5.454545-1*fd.tilt(ind)-1*fd.istilt(ind),'.')
    ind = true(size(fd.schnum)) & fd.wgt == wgtind;
    plot(fd.schnum(ind),fd.phi(ind)-fd.phi_isaac,'.');
    plot(fd.schnum(ind),fd.istilt(ind),'.')
    plot(fd.schnum(ind),fd.tilt(ind),'.')
    ylim([-1 1])
    grid on
end

%% Look the jig data

vals = {fd.temp*1000, fd.istemp*1000, fd.istilttemp*100+273.15};
vals = {fd.N1,fd.N2}
%vals = {fd.A};
%fig = figure(2+wgtind);
clf; hold on;
for valind = 1:length(vals)
    ind = true(size(fd.schnum)) & fd.wgt == wgtind;
    plot(fd.schnum(ind),vals{valind}(ind),'.');
    %plot(fd.phi(ind),vals{valind}(ind),'.');
    % ind = ismember(fd.schnum,[28:31 37:39]) & fd.wgt == wgtind;
    % plot(fd.schnum(ind),vals{valind}(ind),'.');
    % ind = ~ismember(fd.schnum,[28:31 37:39]) & fd.wgt == wgtind;
    % plot(fd.schnum(ind),vals{valind}(ind),'.');
end
%ylim([-1 7])
grid on

%% look at residuals

ind = find(fd.schnum==37 & fd.wgt==3);
figure(6);
clf; hold on;
for rowind = ind
    plot(angs{rowind},reses{rowind})

end

%% Tilt / Temp Plots

plotflag = 1;
temp = [];
figure(4)
clf; hold on;
for di =1:length(data)
    if plotflag==1
        plot(data{di}.Time,polyval(rps_tilt_cal,data{di}.Tilt),'b')
        plot(data{di}.Time,polyval(isaac_tilt_cal,data{di}.ISTilt),'r')
    elseif plotflag==2
        plot(data{di}.Time,data{di}.ISTemp,'r')
        plot(data{di}.Time,data{di}.Temp,'b')
    elseif plotflag==3
        plot(data{di}.Time,data{di}.Tilt,'b')
        plot(data{di}.Time,data{di}.ISTilt,'r')
    end
    xlabel('Time')
    grid on
    %temp = [temp; data{di}.Temp];
end


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

