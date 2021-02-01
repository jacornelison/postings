
[p, p_ind] = get_array_info(datestr(datenum(moonsch.t{1},'yyyy-mmm-dd:HH:MM:SS'),'yyyymmdd'),'obs','obs','obs');
prx = 2 * sind(p.r / 2) .* cosd(p.theta) * 180 / pi;
pry = 2 * sind(p.r / 2) .* sind(p.theta) * 180 / pi;

dirname = ['rps_data/' moonsch.name '/'] ;

load([dirname 'fit_data_temp.mat'])
load([dirname 'fit_data_final.mat'])



% Look at how residuals change when we lock-in the 'drum angle'

fit_ch = [];
for schind = 1:length(md)
    fit_ch = [fit_ch md{schind}.ch];
end

data = [prx(fit_ch); pry(fit_ch)];

% How many parameters?
% [tilt, roll, fpu angle]
fp.free = [1 1 0];
fp.ub = [60 5 5];
fp.lb = [40 -5 -5];

guess = [45.0 0 0-1e-3] + 1e-3;
[param, err, gof] = matmin('chisq',guess, fp,'rps_get_moon_model_v2',data,ones(size(data)),md,moonsch.mount);

model = rps_get_moon_model_v2(param,md,moonsch.mount);

res = reshape(data-model,[],2);

scale = 0.4;
figure(1)
set(gcf,'Position',[100,400,1200,400])
subplot(1,3,1)
plot(res(:,1),res(:,2),'.')
grid on
title({'Best Fit Residuals',...
    'Dk\_offs locked at 0',...
    sprintf('Best Fit Tilt: %2.2f^o Roll: %2.2f^o',param(1),param(2))})
xlim([-1,1]*scale)
ylim([-1,1]*scale)
xlabel('X (deg)')
ylabel('Y (deg)')
axis square

subplot(1,3,2)
plot(fd.resx,fd.resy,'.')
grid on
title({'Best Fit Residuals',...
    'Dk\_offs fit to -1.21',...
    sprintf('Best Fit Tilt: %2.2f^o Roll: %2.2f^o',fd.fitparam(1),fd.fitparam(2))})
xlim([-1,1]*scale)
ylim([-1,1]*scale)
xlabel('X (deg)')
%ylabel('Y (deg)')
axis square

subplot(1,3,3)
plot(res(:,1)-fd.resx,res(:,2)-fd.resy,'.')
grid on
title({'Difference in residuals', 'dk\_locked - dk\_fit'})
xlim([-1,1]*scale)
ylim([-1,1]*scale)
xlabel('X (deg)')
%ylabel('Y (deg)')
axis square


% Okay, now lets get rid of all knowledge of the drum angle and see what
% the fits look like.



fit_ch = [];
for schind = 1:length(md)
    fit_ch = [fit_ch md{schind}.ch];
end

b = find_file_by_date(20170121,'aux_data//beams/beamcen');
brx = 2 * sind(b.r / 2) .* cosd(b.theta) * 180 / pi;
bry = 2 * sind(b.r / 2) .* sind(b.theta) * 180 / pi;

data = [brx(fit_ch); bry(fit_ch)];

% How many parameters?
% [tilt, roll, fpu angle]
fp.free = [1 1 0];
fp.ub = [60 5 5];
fp.lb = [40 -5 -5];

guess = [45.0 0 1.5-1e-3] + 1e-3;
[param, err, gof] = matmin('chisq',guess, fp,'rps_get_moon_model_v2',data,ones(size(data)),md,moonsch.mount);

model = rps_get_moon_model_v2(param,md,moonsch.mount);

res = reshape(data-model,[],2);

scale = 0.4;
figure(1)
set(gcf,'Position',[100,400,1200,400])
subplot(1,3,1)
plot(res(:,1),res(:,2),'.')
grid on
title({'Best Fit Residuals',...
    'Dk\_offs locked at 1.5^o',...
    sprintf('Best Fit Tilt: %2.2f^o Roll: %2.2f^o',param(1),param(2))})
xlim([-1,1]*scale)
ylim([-1,1]*scale)
xlabel('X (deg)')
ylabel('Y (deg)')
axis square

fp.free = [1 1 1];
fp.ub = [60 5 5];
fp.lb = [40 -5 -5];

guess = [45.0 0 0-1e-3] + 1e-3;
[fitparam, err, gof] = matmin('chisq',guess, fp,'rps_get_moon_model_v2',data,ones(size(data)),md,moonsch.mount);

model = rps_get_moon_model_v2(fitparam,md,moonsch.mount);

fitres = reshape(data-model,[],2);

subplot(1,3,2)
plot(fitres(:,1),fitres(:,2),'.')
grid on
title({'Best Fit Residuals',...
    ['Dk\_offs ' sprintf('fit to %2.2f^o',fitparam(3))],...
    sprintf('Best Fit Tilt: %2.2f^o Roll: %2.2f^o',fitparam(1),fitparam(2))})
xlim([-1,1]*scale)
ylim([-1,1]*scale)
xlabel('X (deg)')
%ylabel('Y (deg)')
axis square

subplot(1,3,3)
plot(res(:,1)-fitres(:,1),res(:,2)-fitres(:,2),'.')
grid on
title({'Difference in residuals', 'dk\_locked - dk\_fit'})
xlim([-1,1]*scale)
ylim([-1,1]*scale)
xlabel('X (deg)')
%ylabel('Y (deg)')
axis square




% With the best-fit params, what's the difference between beam centers
% estimated from gaussian fits and beam centers taken from the max
% amplitude?

for schind = 1:3
    
    fname = sprintf('moondata_sch_%03i',schind);
    fprintf(['Loading from:\n' dirname fname '\n'])
    tic
    load([dirname fname])
    toc
    % Select only data with the correct feature bit set.
    [mask_slow, mask_fast] = select_feature_bit(d, 3);
    %mask_fast = true(size(d.mce0.data.fb,1),1);
    
    % Also get rid of skipped samples:
    % They should be zero for all channels, so if the sum of the zeros is equal
    % to the number of channels, the sample is considered skipped.
    ind = find(sum(d.mce0.data.fb == 0,2)==size(d.mce0.data.fb,2));
    mask_fast(ind) = false;
    
    % Calculate downsampled encoder positions.
    az = d.pointing.hor.az(mask_fast);
    el = d.pointing.hor.el(mask_fast);
    dk = d.pointing.hor.dk(mask_fast);
    
    fb = d.mce0.data.fb(mask_fast,:);
    chselect = find(ismember(p_ind.rgl100,md{schind}.ch));
    fb = fb(:,chselect);
    ch = p_ind.rgl100(chselect);
    
    t = d.antenna0.time.utcfast(mask_fast);
    t = (t-t(1))*24*3600; % Start from zero, convert to seconds
    
    % Apply filters to get rid atmospheric effects and noise.
    disp('Applying filters')
    dt = diff(t);
    sampFreq = 1/median(dt(dt>0));
    [b,a] = butter(3,1/300/(sampFreq/2),'high');
    t0 = NaN(size(ch));
    for chind = 1:length(ch)
        tic;
        fb(:,chind) = medfilt1(filtfilt(b,a,fb(:,chind)));
        t0(chind) =  toc;
        fprintf('\rest time remaining: %04.2f minutes',nanmean(t0)*(length(chselect)-chind)/60)
    end
    disp('')
    
    fbstruct{schind}.fb = fb;
    
end


mirror = md{1}.mirror;
bs = md{1}.bs;

mirror.tilt = params(1);
mirror.roll = params(2);
dk_off = params(3);

prx = 2 * sind(p.r / 2) .* cosd(p.theta) * 180 / pi;
pry = 2 * sind(p.r / 2) .* sind(p.theta) * 180 / pi;

for schind = 1:length(md)

    [r, theta, psi] = keck_beam_map_pointing(md{schind}.az,...
        md{schind}.el,...
        md{schind}.dk+dk_off,...
        mount,...
        mirror,...
        md{schind}.source,...
        bs);

    % Convert to x/y using Lambert azimuthal equal area projection
    %[x, y] = pol2cart(theta*pi/180,r);
    md{schind}.x = 2 * sind(r / 2) .* cosd(theta) * 180 / pi;
    md{schind}.y = 2 * sind(r / 2) .* sind(theta) * 180 / pi;
    
end

clr = {[0, 0.4470, 0.7410],...
    [0.8500, 0.3250, 0.0980],...
    [0.4660, 0.6740, 0.1880]*0.8,...
    [0.9290, 0.6940, 0.1250]};
figure(1);
set(gcf,'Position',[100,400,1200,400])
clf; hold on;
for schind = 1:3
    subplot(1,3,schind)
    chind = md{schind}.ch==696;%fit_ch(690);
    if ~isempty(find(chind))
    px = prx(md{schind}.ch(chind));
    py = pry(md{schind}.ch(chind));
    x = md{schind}.x;
    y = md{schind}.y;
    az = md{schind}.az;
    el = md{schind}.el;
    maskind = sqrt((x-px).^2+(y-py).^2) < 1;
    %plot3(x(maskind),y(maskind),fbstruct{schind}.fb(maskind,chind),'Color',clr{schind})
    plot3(el(maskind),az(maskind),fbstruct{schind}.fb(maskind,chind),'Color',clr{schind})
    grid on
    xlabel('x (^o)')
    ylabel('y (^o)')
    end
end




fit_ch = [];
fit_sch = [];
fit_ind = [];
for schind = 1:length(md)
    fit_ch = [fit_ch md{schind}.ch];
    fit_sch = [fit_sch ones(size(md{schind}.ch))*schind];
    fit_ind = [fit_ind 1:length(md{schind}.ch)];
end


tic; [model_beams, val, flag] = rps_get_moon_model_v3(fd.fitparam,md,moonsch.mount,fbstruct,p); toc;
bms_fit = reshape(model_beams,[],2);


data = [prx(fit_ch); pry(fit_ch)];

tic; model = rps_get_moon_model_v2(fd.fitparam,md,moonsch.mount); toc;
res = reshape(data-model,[],2);

res_beams = reshape(data-model_beams,[],2);

res_fits = reshape(model-model_beams,[],2);

scale = 0.15;

figure(3)
set(gcf,'Position',[100,400,1200,400])
subplot(1,3,1)
plot(res_beams(:,1),res_beams(:,2),'.')
xlim([-1,1]*scale)
ylim([-1,1]*scale)
xlabel('X (deg)')
ylabel('Y (deg)')
grid on
axis square

subplot(1,3,2)
plot(res(:,1),res(:,2),'.')
xlim([-1,1]*scale)
ylim([-1,1]*scale)
xlabel('X (deg)')
ylabel('Y (deg)')
grid on
axis square

subplot(1,3,3)
plot(res_fits(:,1),res_fits(:,2),'.')
xlim([-1,1]*scale)
ylim([-1,1]*scale)
xlabel('X (deg)')
ylabel('Y (deg)')
grid on
axis square


figure(4)
scale = 10;
quiver(prx(fit_ch),pry(fit_ch),res_beams(:,1)*scale,res_beams(:,2)*scale,0)


chsub = [];
for i = 1:3
    chsub = [chsub md{i}.ch(1:20)];
end

dsub = [prx(chsub) pry(chsub)];


% Grab an obviously wrong centroid and look at it
%fitind = find(res_beams(:,1)>0.03 & res_beams(:,2)>0.06);
fitind = find(res_fits(:,1)>0.0 & res_fits(:,2)<0);
%fitind = find(sqrt(res_beams(:,1).^2 + res_beams(:,2).^2)>0.06);
schind = fit_sch(fitind);
chind = fit_ind(fitind);

figure(4)
clf
set(gcf,'Position',[100,400,1200,400])
subplot(1,3,1)
hold on
A = fbstruct{schind}.fb(:,chind);
plot(md{schind}.x,A)
plot([1 1]*model_beams(fitind),[min(A) max(A)*1.1],'r') 
plot([1 1]*model(fitind),[min(A) max(A)*1.1],'r') 
grid on
xlabel('x')

subplot(1,3,2)
hold on
A = fbstruct{schind}.fb(:,chind);
plot(md{schind}.y,A)
plot([1 1]*model_beams(length(res_beams)+fitind),[min(A) max(A)*1.1],'r')
plot([1 1]*model(length(res)+fitind),[min(A) max(A)*1.1],'g')
grid on
xlabel('y')


fitind = find(res_fits(:,1)<0.0 & res_fits(:,2)<0);
figure()
hold off; plot(res_fits(:,1),res_fits(:,2),'.')          
hold on; plot(res_fits(fitind,1),res_fits(fitind,2),'r.')
figure()
hold off; plot(bms_fit(:,1),bms_fit(:,2),'.')
hold on; plot(bms_fit(fitind,1),bms_fit(fitind,2),'r.')



% AB offsets?

[p, p_ind] = get_array_info(datestr(datenum(moonsch.t{1},'yyyy-mmm-dd:HH:MM:SS'),'yyyymmdd'),'obs');
prx = 2 * sind(p.r / 2) .* cosd(p.theta) * 180 / pi;
pry = 2 * sind(p.r / 2) .* sind(p.theta) * 180 / pi;


b = find_file_by_date(20170121,'aux_data//beams/beamcen');

ind = p_ind;
if ~isempty(b)
  % Ensure no A/B offsets are added here: to make sure,
  % average the A and B beam centers.
  % lat/lon of detectors, possibly with A/B offsets
  [lat,lon]=reckon(0,0,b.r,b.theta+90);
  % distance and angle of offset
  [R,AZ]=distance(lat(ind.a),lon(ind.a),lat(ind.b),lon(ind.b));
  % step back by half of offset
  [lat(ind.a),lon(ind.a)]=reckon(lat(ind.a),lon(ind.a),R/2,AZ);
  lat(ind.b)=lat(ind.a);
  lon(ind.b)=lon(ind.a);

  % go back to r/theta
  b.r=distance(0,0,lat,lon);
  b.theta=azimuth(0,0,lat,lon)-90;

  % Could do same thing using 'paircenter' function
  % p.r=b.r;
  % p.theta=b.theta;
  % [p.r(ind.a),p.theta(ind.a)]=paircenter(b.r(ind.a),b.theta(ind.a),b.r(ind.b),b.theta(ind.b));
  % p.r(ind.b)=p.r(ind.a);
  % p.theta(ind.b)=p.theta(ind.a);
end

brx = 2 * sind(b.r / 2) .* cosd(b.theta) * 180 / pi;
bry = 2 * sind(b.r / 2) .* sind(b.theta) * 180 / pi;

figure(4)
clf; hold on
for i = 1:2
    if i==1
        ind = p_ind.rgl100a;
    else
        ind = p_ind.rgl100b;
    end
    
    quiver(brx(ind),bry(ind),prx(ind)-brx(ind),pry(ind)-bry(ind),0,'Color',clr{i})


end


%% Look at the moon dataset

figure(2)
set(gcf,'Position',[1000 100 600 600])
clf; hold on;
for i = 1:3
plot(md{i}.x,md{i}.y,'Color',clr{i})
end
plot(brx,bry,'kx')
grid on
xlabel('x (^o)')
ylabel('y (^o)')
title('Moon Observation Coverage Jan 2017')
legend({'DK = 0','DK = 45','DK = 90','CMB-derived pnts'})
axis square


% Bootstrap the beam fitting

freepar.free = [1 1 1 1 1 1 1];
freepar.lb = [-1e4 -1e3 -1e3 0 0 -1 -100];
freepar.ub = [1e4 1e3 1e3 2 2 1 100];


N=100;
bootparm = NaN(N,7);
for schind = 1%:3
    chind = md{schind}.ch==696;%fit_ch(690);
    if ~isempty(find(chind))
        px = prx(md{schind}.ch(chind));
        py = pry(md{schind}.ch(chind));
        x = md{schind}.x;
        y = md{schind}.y;
        
        centind = md{schind}.centind(chind);
        
        
        A = fbstruct{schind}.fb(md{schind}.centind(chind),chind);
        px = prx(md{schind}.ch(chind));
        py = pry(md{schind}.ch(chind));
        guess = [A px py 0.15 0.15 0 0];
        %maskind = sqrt((x-x(centind)).^2+(y-y(centind)).^2) < 1;
        %maskind = sqrt((x-px).^2+(y-py).^2) < 1;
        
        for i = 1:N
            % Bootstrap funny business.
            % Mask data beyond 1 deg of the beam.
            % Then select random points of the mask with replacement.
            maskind = find(sqrt((x-px).^2+(y-py).^2) < 1);
            s = size(maskind);
            bootind = randi([1,max(s)],s(1),s(2));
            maskind = maskind(bootind);
            
            %fun = @(prm)chimin(prm,fbstruct{schind}.fb(maskind,chind),x(maskind),y(maskind));
            %[parm, fval, eflag] = fminsearch(fun,guess);
            [parm, err, gof, stat, cov] = matmin('chisq',guess, freepar,'egauss2',...
                fbstruct{schind}.fb(maskind,chind), ones(size(fbstruct{schind}.fb(maskind,chind))),...
                x(maskind),y(maskind));
            
            bootparm(i,:) = parm;
        end
        
    end
end

% Fix starpoint pointing model

for schind=2:3
    fname = sprintf('moondata_sch_%03i',schind);
    load([dirname fname]);    

    pm = ParameterRead('aux_data/pointingmodels/pointingmodel_complete_v1.csv');
    pm = structcut(pm,moonsch.strpnt);
    d = invpointing_model(d, pm);
    d.ch = ch;
    
    % Select only data with the correct feature bit set.
    [mask_slow, mask_fast] = select_feature_bit(d, 3);
    %mask_fast = true(size(d.mce0.data.fb,1),1);
    
    % Also get rid of skipped samples:
    % They should be zero for all channels, so if the sum of the zeros is equal
    % to the number of channels, the sample is considered skipped.
    ind = find(sum(d.mce0.data.fb == 0,2)==size(d.mce0.data.fb,2));
    mask_fast(ind) = false;
    
    % Calculate downsampled encoder positions.
    az = d.pointing.hor.az(mask_fast);
    el = d.pointing.hor.el(mask_fast);
    dk = d.pointing.hor.dk(mask_fast);
    
    md{schind}.az=az;
    md{schind}.el=el;
    md{schind}.dk=dk;
    
end


