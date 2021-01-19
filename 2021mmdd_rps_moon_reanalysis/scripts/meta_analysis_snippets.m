
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

guess = [45.0 0 -1.5-1e-3] + 1e-3;
[param, err, gof] = matmin('chisq',guess, fp,'rps_get_moon_model_v2',data,ones(size(data)),md,moonsch.mount);

model = rps_get_moon_model_v2(param,md,moonsch.mount);

res = reshape(data-model,[],2);

scale = 0.2;
figure(1)
set(gca,'Position',[100,400,1200,400])
subplot(1,3,1)
plot(res(:,1),res(:,2),'.')
grid on
title({'Best Fit Residuals','Dk\_offs locked at 1.5'})
xlim([-1,1]*scale)
ylim([-1,1]*scale)
xlabel('X (deg)')
ylabel('Y (deg)')
axis square

subplot(1,3,2)
plot(fd.resx,fd.resy,'.')
grid on
title({'Best Fit Residuals','Dk\_offs fitted to 1.21'})
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

[xmodel, ymodel] = deal([]);
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



fit_ch = [];
fit_sch = [];
fit_ind = [];
for schind = 1:length(md)
    fit_ch = [fit_ch md{schind}.ch];
    fit_sch = [fit_sch ones(size(md{schind}.ch))*schind];
    fit_ind = [fit_ind 1:length(md{schind}.ch)];
end

data = [prx(fit_ch); pry(fit_ch)];

tic; model = rps_get_moon_model_v2(fd.fitparam,md,moonsch.mount); toc;
res = reshape(data-model,[],2);

tic; model_beams = reshape(rps_get_moon_model_v3(fd.fitparam,md,moonsch.mount,fbstruct,p),[],2); toc;
model_beams = reshape(model_beams,[],1);

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


chsub = [];
for i = 1:3
    chsub = [chsub md{i}.ch(1:20)];
end

dsub = [prx(chsub) pry(chsub)];
    
   
% Grab an obviously wrong centroid and look at it
%fitind = find(res_beams(:,1)>0.03 & res_beams(:,2)>0.06);
fitind = find(res_fits(:,1)>0.0 & res_fits(:,2)>0.1);
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




