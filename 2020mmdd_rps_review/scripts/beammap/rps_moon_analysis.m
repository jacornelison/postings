%% Start:
% This code acts as a preliminary analysis of moon observations that were
% made in early 2017. We have measurements at 4 decks, but I've only
% analyzed DK = 0 right now.
%
% We also have data for all of tiles 9 through 13, but I only use 12
% channels for the prelim analysis. When I have time, I'll incorporate more
% channels.
%
% Once I derived the mirror parameters from the moon obs, I use RPS data to
% determine the parameters of the source assuming the mirror parameters are
% 100% correct. The timestreams are much smaller so I use all Pol A tiles
% from 9 to 13.
%
% Current best-fit params as of 19 Mar 2019:
% 
rpsopt.mirror.tilt = 44.50514;
rpsopt.mirror.roll =  0.09811;

rpsopt.source.azimuth = -177.6500060189; % +/- 4.5e-10
rpsopt.source.el = 2.3129580213; % +/- 1.6e-10
rpsopt.source.distance = 195.5; % +/- 1
rpsopt.source.height = rpsopt.source.distance*tand(rpsopt.source.el);


%% Load moon data

dirname = '~jcornelison/rps_data/2016_moon/';
[p, p_ind] = get_array_info(20191101,'obs','obs','obs');
prx = 2 * sind(p.r / 2) .* cosd(p.theta) * 180 / pi;
pry = 2 * sind(p.r / 2) .* sind(p.theta) * 180 / pi;

% We'll need all three deck angles. But for now, we're just doing the
% first.


if 0
    % From obs log
    t1 = '20170121 22:00:00';
    t2 = '20170121 23:30:00';
    
    [msg, t] = search_run_log('moonraster',t1,t2,1,'log/');
    
    t1 = t(1);
    t2 = t(2);
    
    % Downselect to however many channels we want.
    % There's a lot of data, so this takes a while to load.
    % Picking 5 known good pairs from T11 for now.
    ch = intersect(find(ismember(p.tile,9:13)),p_ind.rgl100);
    
    % We should include some tiles that are farther out.
    til = [09 11 11 11 11 13];
    col = [03 07 07 02 05 05]; %
    row = [03 01 08 06 05 05]; %
    ch = [];
    for i = 1:length(col)
        ch = [ch find(p.tile == til(i) & p.det_col==col(i) & p.det_row == row(i))];
    end
    ch = reshape(ch,1,[]);
    
    ch = intersect(find(ismember(p.tile,9:13)),p_ind.rgl100);
    ch = sort(ch);
    
    
    disp('Selecting Registers')
    % Select registers to read using load_arc.
    %   * array contains feature bits
    %   * antenna0 contains telescope encoders, timestamps
    %   * mceX.data.fb for selected channels
    %   * mceX.cc, mceX.rc1, mceX.frame are necessary for deconvolution
    reg = {'array', 'antenna0'};
    
    expt = get_experiment_name();
    
    if strcmp(expt,'keck')
        for rx=unique(p.rx')
            mceselect = sprintf('mce%i', rx);
            mce_channels = ch(floor((ch - 1) / 528) == rx) - 528 * rx;
            if length(mce_channels) > 0
                chselect = [mceselect sprintf('.data.fb(%i', mce_channels(1))];
                for i=2:length(mce_channels)
                    chselect = [chselect ',' sprintf('%i', mce_channels(i))];
                end
                chselect = [chselect ')'];
            else
                chselect = [mceselect '.data.fb()'];
            end
            reg{length(reg) + 1} = chselect;
            reg{length(reg) + 1} = [mceselect '.cc'];
            reg{length(reg) + 1} = [mceselect '.rc1'];
            reg{length(reg) + 1} = [mceselect '.frame'];
        end
    elseif strcmp(expt,'bicep3')
        for i=unique(p.mce')
            ind = intersect(ch,find(p.mce==i))';
            mceselect = sprintf('mce%i',i);
            chselect = [mceselect '.data.fb('];
            chselect = [chselect sprintf('%i,',1+p.gcp(ind))];
            chselect = [chselect(1:end-1) ')'];
            reg{length(reg) + 1} = chselect;
            reg{length(reg) + 1} = [mceselect '.cc'];
            reg{length(reg) + 1} = [mceselect '.rc1'];
            reg{length(reg) + 1} = [mceselect '.frame'];
        end
    else
        disp('Unrecognized experiment name. We can only reduce keck and bicep3 data right now.')
        return
    end
    
    
    d = load_arc('arc',t1,t2,reg);
    
    d = make_utc_single_col(d);
    
    % Combine all mce into mce0 field of the data structure.
    if ~isfield(d, 'mce0')
        d.mce0 = [];
    end
    for m = {'mce1', 'mce2', 'mce3', 'mce4'}
        if isfield(d, m)
            d.mce0 = structcat(2, [d.mce0, getfield(d, m{1})]);
            d = rmfield(d, m);
        end
    end
    
    % Ran into an issue where the DK encoder wasn't being properly offset in
    % the pointing model. So set the encoder cals and offsets to the 'Normal'
    % values and let keck_beam_map_pointing do its job.
    
    % Schedule encoder_zeros -52.621, -7.0, 2.626
    % Normal encoder_zeros 127.37899, -82.41876, -2.626 # 2016 Feb 2 - BICEP3 from 20160122
    % Schdule encoder_cals: 2304000, -2304000, -574400, ant=all # az,el,dk counts per turn, reverse el&dk direction for ffflat
    % Normal encoder_cals: 2304000, 2304000, 574400
    
    d.antenna0.tracker.encoder_off(:,3) = ...
        -1 * d.antenna0.tracker.encoder_off(:,3);
    d.antenna0.tracker.encoder_mul(1,:) = [2304000, 2304000, 574400];
    d.antenna0.tracker.encoder_off = repmat(([127.37899, -82.41876, -2.626])*3.6e6 ,size(d.antenna0.tracker.encoder_off,1),1);
    
    pm = get_pointing_model(d.antenna0.time.utcfast(1,1), 0, d);
    pm.az_tilt_ha = 0;
    pm.az_tilt_lat = 0;
    pm.el_tilt = 0;
    d = invpointing_model(d, pm);
    d.ch = ch;
    
    save([dirname 'moondata_dk000'],'d','-v7.3')
    
else
    
    load([dirname 'moondata_dk000'])
    ch = d.ch;
end

md = []; %moon data
md.tod = {};
md.az     = {};
md.el     = {};
md.dk     = {};
md.dk0    = [];
md.ch     = [];


% Manually find good channels that don't have flux jumps, noise, etc... 
good_chans = [134 135 136 137 140 141 226 227 160 161 179 180 ...
    668 669 672 673 684 685 686 687 692 693 696 697 708 709 714 715 ...
    928 929 934 935 940 941 970 971 978 979];
%good_chans = find(ismember(p.tile,13));


for i = 1:length(ch)
    
    if any(ismember(ch(i),good_chans))
        md.tod{end+1} = d.mce0.data.fb(:,ch==ch(i));
        md.az{end+1} = d.pointing.hor.az;
        md.el{end+1} = d.pointing.hor.el;
        md.dk{end+1} = d.pointing.hor.dk;
        md.dk0(end+1) = nanmean(d.pointing.hor.dk);
        md.ch(end+1) = ch(i);
    end
end

% Grab cmb-derived beam centroids
% This is the 'data' for matmin.m
xy0 = [prx(md.ch); pry(md.ch)];


%% Set up pointing model parameters

moonopt = [];
moonopt.mount.aperture_offr = 0;
moonopt.mount.aperture_offz = 0.9970;
moonopt.mount.dk_offy = 0;
moonopt.mount.dk_offx = 0;
moonopt.mount.el_tilt = 0;
moonopt.mount.el_offx = 0;
moonopt.mount.el_offz = 0;
moonopt.mount.az_tilt_ha = 0;
moonopt.mount.az_tilt_lat = 0;

moonopt.mirror.height = 0.9540;
moonopt.mirror.tilt = 44.35;
moonopt.mirror.roll = -1;

% Boresight 'fpu_data'
% keck_beam_map_pointing, by default, gives you detector centered
% coordinates (x'/y') for every detector in p = get_array_info.
% Substituting 'bs' for 'p' will return a single timestream in
% boresight-centered coordinates (x/y);

% This is useful because A. less data to deal with and 
% B. beam centers centers will be located at a given detector's r/theta
% instead of zero which is useful for taking beam-center residuals.
bs.r = 0;
bs.theta = 0;
bs.chi = 0;
bs.chi_thetaref = 0;
bs.drumangle = 0;
mount.aperture_offr = 0;

%% Let's figure out the source position
% Get moon ephemeris data from jpl.horizons database
f = fopen('~jcornelison/rps_data/2016_moon/moon_ephem_2017_dk000.csv','r');
moondata = textscan(f,'%s %n %n %n %n %s','Delimiter',',');
fclose(f);

% moondata comes in RA and DEC. We need to convert this to AZ/EL

[yr,mo,da,hr,mi] = datevec(moondata{1},'yyyy-mmm-dd:HH:MM');
sec = zeros(size(mi));
tmoon = date2mjd(yr,mo,da,hr,mi,sec);
[az, el] = radec2azel(moondata{2}+180,moondata{3},tmoon,-89.9911,-44.65);

% Interpolate source position for timestream
t = d.antenna0.time.utcfast;
az = interp1(tmoon,az,t);
el = interp1(tmoon,el,t);
moondist = interp1(tmoon,moondata{4}*1000,t);

% then convert El to source height off using the distance data supplied in column moondata{4}
moonopt.source.azimuth = az;
moonopt.source.el = el;
moonopt.source.distance = moondist.*cosd(el);
moonopt.source.height = moondist.*sind(el);

% Didn't we already calculate this?
[p, p_ind] = get_array_info(20170131,'obs','obs','obs');
%[p, p_ind] = get_array_info(20170131,'ideal');
prx = 2 * sind(p.r / 2) .* cosd(p.theta) * 180 / pi;
pry = 2 * sind(p.r / 2) .* sind(p.theta) * 180 / pi;

% Params from fit: tilt = 44.4735955 +/- 1.3880e-7, roll = 9.94091e-2 +/- 5.5535e-7  
moonopt.mirror.tilt = 44.51078;
moonopt.mirror.roll = 0.11346;

[r, theta, psi] = ...
    keck_beam_map_pointing(d.pointing.hor.az, d.pointing.hor.el, ...
    d.pointing.hor.dk, moonopt.mount, moonopt.mirror, moonopt.source, bs);

x = 2 * sind(r / 2) .* cosd(theta) * 180 / pi;
y = 2 * sind(r / 2) .* sind(theta) * 180 / pi;

good_chans = [134 135 136 137 140 141 226 227 160 161 179 180 ...
    668 669 672 673 684 685 686 687 692 693 696 697 708 709 714 715 ...
    928 929 934 935 940 941 970 971 978 979];

% I use these to check the current parameters.
[mx,my] = deal([]);
fig = figure(1);
hold off
for i = 1:length(md.ch)
    %if ~any(ismember(md.ch(i),[705 710])) 
    data = md.tod{i};
    ind = abs(x-prx(md.ch(i))) < 2 & abs(y-pry(md.ch(i))) < 2;
    xi = x(ind);
    yi = y(ind);
    data = data(ind);
    [m,mind] = max(data);
    guess = [m prx(md.ch(i)) pry(md.ch(i)) p.fwhm_maj(md.ch(i))^2 p.fwhm_maj(md.ch(i))^2 0 mean(data)];
    fun = @(z)sum((data-egauss2(z,xi,yi)).^2);    
    %fun = @(z)chimin(z,tod,x,y);
    
    [T, param] = evalc('fminsearch(fun,guess)');
    
    mx(end+1) = param(2);
    my(end+1) = param(3);
    
    %mx(end+1) = xi(mind);
    %my(end+1) = yi(mind);
    
    if 0 
    hold off
    plot3(xi,yi,data)
    hold on
    plot3([1 1]*prx(md.ch(i)),[1 1]*pry(md.ch(i)),[0,3000],'r')
    %plot3([1 1]*xi(mind),[1 1]*yi(mind),[0,3000],'g')
    plot3([1 1]*param(2),[1 1]*param(3),[0,3000],'g')
    
    title(sprintf('channel %d, tile %d, row %d, col %d',md.ch(i),...
        p.tile(md.ch(i)), p.det_row(md.ch(i)), p.det_col(md.ch(i))))
        %pause(0.5)
    zzz = input(sprintf('dx %d dy %d',prx(md.ch(i))-param(2),pry(md.ch(i))-param(3)));
    end

end

if 1 
fig = figure(2);
clf(fig)
hold off
plot(x,y)
hold on
%plot(mirror_model(1:end/2),mirror_model((end/2+1):end),'kx')
plot(mx,my,'gx')
plot(prx(md.ch),pry(md.ch),'rx')
plot(prx,pry,'rx')
plot(mx,my,'gx')
xlim([-20 20])
ylim([-20 20])
xlabel('x (^o)')
ylabel('y (^o)')
grid on
end


%% Fit stuff
% I fit for the mirror parameters that minimize residuals between CMB- and
% moon derived beam centers. This currently requires fitting each beam
% every time a parameter changes. Thus, there's three beamfitting options:

% 'max' assumes the max value in the timestreams is the beamcenter. This is
% really quick, but gives wrong answers if the real beam is far away from
% your guess.

% 'matmin' does a chi-square minimization on each beam using matmin.m. Has
% a lot of text output and I'm guessing matmin doesn't like to be nested as
% this option will fail if you try to do a chi-square min on the mirror
% params. Use only if you're doing one-offs or a grid search on the mirror
% params.

% 'fminsearch' uses fminsearch.m for the beam fit. Slowest of the options,
% but not by much and doesn't fail if you try to do a minimization on the 
% mirror params.

fp.free = [1 1];
fp.ub = [44.7 1];
fp.lb = [44.3 -1];

guess = [44.51078 0.11346];

%moonopt.fittype = 'max';
%moonopt.fittype = 'matmin';
moonopt.fittype = 'fminsearch';

[param, err, gof] = matmin('chisq',guess, fp,'rps_get_moon_model',xy0,ones(size(d))*1e-4,md,moonopt,p);

%% Chi2 grid map

ub = [44.7 1];
lb = [44.3 -1];
n = 5;

param1 = lb(1):(ub(1)-lb(1))/(n-1):ub(1);
param2 = lb(2):(ub(2)-lb(2))/(n-1):ub(2);

mirrchi = NaN(length(param1),length(param2));
%mirrchi = normrnd(0,1,length(param1),length(param2));
%mirrchi = rand(length(param1),length(param2));

moonopt.fittype = 'max';
for i = 1:length(param1)
    for j = 1:length(param2)
        mirror_model = rps_get_moon_model([param1(i) param2(j)],md,moonopt,p);
        mirrchi(i,j) = nansum((xy0-mirror_model).^2);
    end
end


fig = figure(1);
clf(fig)
imagesc(param1,param2,mirrchi);%,[0,400]);
axis square
colorbar()

tic
moonopt.fittype = 'fminsearch';
[mirror_model, mx, my] = rps_get_moon_model([44.51078 0.11346],md,moonopt,p);
toc
nansum((xy0-mirror_model).^2)

[P1, P2] = meshgrid(param1,param2);
P1 = reshape(P1,1,[]);
P2 = reshape(P2,1,[]);
chimap = reshape(mirrchi,1,[]);


%% Fitting source params
% Orig. source params: az: -177.38, el: 1.58, dist: 210
% Fit params from a single row: az: -177.38, el: 2.30, dist: 210

ch = [148 150 245 247 544 546 548 641 643 759 767 777 785 883 891 899 917 931 975]';
data = rps_get_timestream(sch,1,1,p,rpsopt,ch)

md = []; %moon data
md.tod = {};
md.az     = {};
md.el     = {};
md.dk     = {};
md.dk0    = [];
md.ch     = [];


for j = 1%:3
    for i = 1:length(ch)
        md.tod{end+1} = data.mce0.data.fb(:,ch==ch(i));
        md.az{end+1} = data.pointing.hor.az;
        md.el{end+1} = data.pointing.hor.el;
        md.dk{end+1} = data.pointing.hor.dk;
        md.dk0(end+1) = nanmean(data.pointing.hor.dk);
        md.ch(end+1) = ch(i);
        
    end
end


% Grab cmb-derived beam centroids
xy0 = [prx(md.ch); pry(md.ch)];

% Boresight 'fpu_data'
bs.r = 0;
bs.theta = 0;
bs.chi = 0;
bs.chi_thetaref = 0;
bs.drumangle = 0;

moonopt.mirror.tilt = 44.51078;
moonopt.mirror.roll = 0.11346;

moonopt.source.azimuth = -177.38;
moonopt.source.el = 2.30;
moonopt.source.distance = 210;
moonopt.source.height = moonopt.source.distance*tand(moonopt.source.el);

[r, theta, psi] = ...
    keck_beam_map_pointing(data.pointing.hor.az, data.pointing.hor.el, data.pointing.hor.dk, moonopt.mount, ...
                           moonopt.mirror, moonopt.source, bs);
                       

x = 2 * sind(r / 2) .* cosd(theta) * 180 / pi;
y = 2 * sind(r / 2) .* sind(theta) * 180 / pi;

[mx,my] = deal([]);
fig = figure(2);
hold off
for i = 1:length(ch)
    tod = md.tod{i};
    ind = abs(x-prx(md.ch(i))) < 5 & abs(y-pry(md.ch(i))) < 5;
    xi = x(ind);
    yi = y(ind);
    tod = tod(ind);
    [m,mind] = max(tod);
    mx(end+1) = xi(mind);
    my(end+1) = yi(mind);
%     hold off
%     plot3(xi,yi,data)
%     hold on
%     plot3([1 1]*prx(ch(i)),[1 1]*pry(ch(i)),[0,3000],'r')
%     pause(0.5)
end

if 1 
fig = figure(2);
clf(fig)
hold off
plot(x,y)
hold on
plot(prx(ch),pry(ch),'rx')
plot(mx,my,'gx')
xlim([-20 20])
ylim([1 3])
grid on
end

fig = figure(1);
for i = 1:length(ch)
    clf(fig)
    plot3(x,y,data.mce0.data.fb(:,i))
    grid on
    title(sprintf('CH: %d',data.ch(i)))
    zzz = input('Hit enter to continue');
end

%% Run reduc with preliminary params
% We got the centers pretty close. Now let's run rps_reduc_driver and use
% demodulated data.
% Mirror Params from dk 0: tilt: 44.51078 roll: 0.11346
% Source params from a single row: az: -177.38, el: 2.30, dist: 210

rpsopt.min_dist = 0.25; % Accept only tods with beam centers in LOS;
rpsopt.mirror.tilt = 44.51078;
rpsopt.mirror.roll = 0.11346;

rpsopt.source.azimuth = -177.65858;
rpsopt.source.el = 2.275;
rpsopt.source.distance = 210;
rpsopt.source.height = rpsopt.source.distance*tand(rpsopt.source.el);

dirname = '~jcornelison/rps_data/2016_original/';
load([dirname 'b3rpsfiles'])
rps_reduc_driver(sch(1),p,rpsopt,dirname,(0:7)*13+1)


%% Now with demodulated data

%ch = [148 150 245 247 544 546 548 641 643 759 767 785 883 891 899 917 931]';
%rpstod = rps_read(sch{1}.scans(1),ch,p,rpsopt);


% 

[p, p_ind] = get_array_info(20170131,'obs','obs','obs');
%[p, p_ind] = get_array_info(20170131,'ideal');
prx = 2 * sind(p.r / 2) .* cosd(p.theta) * 180 / pi;
pry = 2 * sind(p.r / 2) .* sind(p.theta) * 180 / pi;

tods = rps_load_tods(sch,1,(0:7)*13+1,dirname);

md = []; %moon data
md.tod = {};
md.az     = {};
md.el     = {};
md.dk     = {};
md.dk0    = [];
md.ch     = [];
md.ind = [];

for j = 1:length(tods)
rpstod = tods{j};
ch = rpstod.ch;

% We're at dk0 so grab A dets (or B dets from MCE0).
ind = find((p.mce==0 & strcmp(p.pol,'A')) | (ismember(p.mce,1:3) & strcmp(p.pol,'B')));
ch = intersect(ch,ind);

cutchans = [152 243 733 741 751 865 873 909 939 949 957 965 542 713 809 ...
    1043 731 739 747 783 863 871 879 915 955 963 154 170 220 230 550 ...
    626 689 1019 1045 154 164 234 616 634 629 637 745 774 781 861 869 ...
    877 737 953 961 156 194 636 681 953 961 1003 604 681 735 743 769 ...
    772 867 875 904 951 959 967 977 1023 177 214 1041];
ch = setdiff(ch,cutchans);

    for i = 1:length(ch)
        md.tod{end+1} = rpstod.todcos(:,rpstod.ch==ch(i));
        md.az{end+1} = rpstod.az;
        md.el{end+1} = rpstod.el;
        md.dk{end+1} = rpstod.dk;
        md.dk0(end+1) = nanmean(rpstod.dk);
        md.ch(end+1) = ch(i);
        md.ind(end+1) = j;
    end
end

% Grab cmb-derived beam centroids
xy0 = [prx(md.ch); pry(md.ch)];

bs.r = 0;
bs.theta = 0;
bs.chi = 0;
bs.chi_thetaref = 0;
bs.drumangle = 0;

moonopt.mirror.tilt = 44.51078;
moonopt.mirror.roll = 0.11346;

moonopt.source.azimuth = -177.65001;
moonopt.source.el = 2.3129580123;
moonopt.source.distance = 195.5;
moonopt.source.height = moonopt.source.distance*tand(moonopt.source.el);

x = {};
y = {};
tic
for i = 1:length(md.ch)
if i == 1 | ((md.ind(i-1)-md.ind(i))<0)
[r, theta, psi] = ...
    keck_beam_map_pointing(md.az{i}, md.el{i}, md.dk{i}, moonopt.mount, ...
                           moonopt.mirror, moonopt.source, bs);
end
x{end+1} = 2 * sind(r / 2) .* cosd(theta) * 180 / pi;
y{end+1} = 2 * sind(r / 2) .* sind(theta) * 180 / pi;
end
toc

[mx,my] = deal([]);
fig = figure(2);
hold off
tic
for i = 1:length(md.ch)
    tod = md.tod{i};
    ind = abs(x{i}-prx(md.ch(i))) < 5 & abs(y{i}-pry(md.ch(i))) < 5;
    xi = x{i}(ind,1);
    yi = y{i}(ind,1);
    tod = tod(ind);
    [m,mind] = max(tod);
    %guess = [m xi(mind) yi(mind) p.fwhm_maj(md.ch(i))^2 p.fwhm_maj(md.ch(i))^2 0];
    guess = [m prx(md.ch(i)) pry(md.ch(i)) p.fwhm_maj(md.ch(i))^2 p.fwhm_maj(md.ch(i))^2 0];
    fun = @(z)sum((tod-egauss2(z,xi,yi)).^2);    
    %fun = @(z)chimin(z,tod,x,y);
    
    [T, param] = evalc('fminsearch(fun,guess)');
    
    mx(end+1) = param(2);
    my(end+1) = param(3);
    %mx(end+1) = xi(mind);
    %my(end+1) = yi(mind);
    
    if 0
    hold off
    plot3(xi,yi,tod)
    hold on
    plot3([1 1]*prx(md.ch(i)),[1 1]*pry(md.ch(i)),[0,200],'r')
    plot3([1 1]*mx(end),[1 1]*my(end),[0,200],'g')
    title(md.ch(i))
    grid on
    zzz = input('');
    end
end
toc

if 1   
fig = figure(2);
clf(fig)
hold off
for i = 1:length(md.ch)
plot(x{i},y{i})
hold on
end
plot(mx,my,'kx')
%plot(source_model(1:end/2),source_model((end/2+1):end),'gx')
plot(prx,pry,'rx')
plot(prx(md.ch),pry(md.ch),'gx')
xlim([-20 20])
ylim([-20 20])
xlabel('x (^o)')
ylabel('y (^o)')
grid on
end

fig = figure(1);
clf(fig)
qscale = 10;
quiver(prx(md.ch),pry(md.ch),(prx(md.ch)-mx')*qscale,(pry(md.ch)-my')*qscale,0,'k')
xlim([-15 15])
ylim([-15 15])
xlabel('x (^o)')
ylabel('y (^o)')
grid on

fig = figure(1);
clf(fig)
edg = -0.1:0.01:0.1;
subplot(1,2,1)
bar(edg,histc([prx(md.ch)]-[mx'],edg),'histc')
xlim([-0.1 0.1])
ylim([0 80])
xlabel('dx (^o)')
ylabel('N')
grid on
subplot(1,2,2)
bar(edg,histc([pry(md.ch)]-[my'],edg),'histc')
xlim([-0.1 0.1])
xlabel('dy (^o)')
ylim([0 80])
grid on


%% Grid fit
% source distance is no closer than 180m and no farther than 210m
% see http://bicep.rc.fas.harvard.edu/bkcmb/analysis_logbook/analysis/20130307_station_layout_2/
% Using a hi-res map, the physical source is 195.5 +/- 1m

ub = [-177.65 2.34 196.5];
lb = [-177.67 2.3 194.5];
n = 10;

param1 = lb(1):(ub(1)-lb(1))/(n-1):ub(1);
param2 = lb(2):(ub(2)-lb(2))/(n-1):ub(2);
%param3 = lb(3):(ub(3)-lb(3))/(n-1):ub(3);

srcchi = NaN(length(param1),length(param2));
moonopt.fittype = 'max';
tic
for i = 1:length(param1)
    for j = 1:length(param2)
        source_model = rps_get_source_model([param1(i) param2(j) 195.5],md,moonopt,p);
        srcchi(i,j) = nansum((xy0-source_model).^2);
    end
end
toc

fig = figure(1);
clf(fig)
imagesc(param2,param1,(srcchi));%,[0,5]);
axis square
colorbar()


moonopt.fittype = 'fminsearch';
tic
[source_model, mx, my] = rps_get_source_model([-177.65001 2.3129580123 195.5],md,moonopt,p);
nansum((xy0-source_model).^2)
toc

[P1, P2] = meshgrid(param1,param2);
P1 = reshape(P1,1,[]);
P2 = reshape(P2,1,[]);
chimap = reshape(mirrchi,1,[]);



%% min-chi-square Fit

% Model = [Az El];
fp.free = [1 1 0];
fp.ub = [-177.65 2.4 210];
fp.lb = [-178.67 2.25 180];

guess = [-177.65858 2.318 195.5];

%moonopt.fittype = 'max';
moonopt.fittype = 'fminsearch';

%tic
%[param, err, gof] = matmin('chisq',guess, fp,'rps_get_source_model',xy0,ones(size(xy0))*1e-4,md,moonopt,p);
%toc

% Using fminsearch took 9 hours to complete.
% Source params from tiles 9 to 13: az: -177.65001, el: 2.3129580123, dist: 195.5



