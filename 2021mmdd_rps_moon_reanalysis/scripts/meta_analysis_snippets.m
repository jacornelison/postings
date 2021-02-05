
[p, p_ind] = get_array_info(datestr(datenum(moonsch.t{1},'yyyy-mmm-dd:HH:MM:SS'),'yyyymmdd'),'obs','obs','obs');
prx = 2 * sind(p.r / 2) .* cosd(p.theta) * 180 / pi;
pry = 2 * sind(p.r / 2) .* sind(p.theta) * 180 / pi;

dirname = ['rps_data/' moonsch.name '/'] ;

load([dirname 'fit_data_temp.mat'])
load([dirname 'fit_data_final.mat'])




fit_ch = [];
for schind = 1:length(md)
    fit_ch = [fit_ch md{schind}.ch];
end

data = [prx(fit_ch); pry(fit_ch)];

% How many parameters?
% [tilt, roll, fpu angle]
fp.free = [1 1];
fp.ub = [60 5];
fp.lb = [40 -5];

guess = [45.0 0] + 1e-6;
%[param, err, gof] = matmin('chisq',guess, fp,'rps_get_moon_model_v2',data,ones(size(data)),md,moonsch.mount);

model = rps_get_moon_model_v2(param,md,moonsch.mount);

res = reshape(data-model,[],2);


for schind = 1%2:length(moonsch.dk)
    fname = sprintf('moondata_sch_%03i.mat',schind);
    fprintf(['Loading from:\n' dirname fname '\n'])
    tic
    load(fullfile(dirname,fname))
    toc
    ch = d.ch;
    
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
    
    
    fb = d.mce0.data.fb(mask_fast,ismember(ch,md{schind}.ch));
    t = d.t(mask_fast);
    t = (t-t(1))*24*3600; % Start from zero, convert to seconds
    ch = md{schind}.ch;
    
    % Apply filters to get rid atmospheric effects and noise.
    % High-pass filter for atmo effects.
    % Median filter for noise.
    disp('Applying filters')
    dt = diff(t);
    sampFreq = 1/median(dt(dt>0));
    [b,a] = butter(3,1/300/(sampFreq/2),'high');
    t0 = NaN(size(ch));
    for chind = 1:length(ch)
        tic;
        fb(:,chind) = medfilt1(filtfilt(b,a,fb(:,chind)));
        t0(chind) =  toc;
        fprintf('\rest time remaining: %04.2f minutes',nanmean(t0)*(length(ch)-chind)/60)
    end
    disp('')
    
    fbstruct{schind}.fb = fb;
end

% This takes fucking forever.
%[model, val, flag] = rps_get_moon_model_v3(fd.fitparam,md,mount,fbstruct,p)
mirror = md{1}.mirror;
bs = md{1}.bs;

mirror.tilt = fd.fitparam(1);
mirror.roll = fd.fitparam(2);

for schind = 1:length(md)
    [r, theta, psi] = keck_beam_map_pointing(md{schind}.az,...
        md{schind}.el,...
        md{schind}.dk,...
        moonsch.mount,...
        mirror,...
        md{schind}.source,...
        md{schind}.bs);
    md{schind}.x = 2 * sind(r / 2) .* cosd(theta) * 180 / pi;
    md{schind}.y = 2 * sind(r / 2) .* sind(theta) * 180 / pi;
end

modsplit = reshape(model,[],2);
modelx = modsplit(:,1);
modely = modsplit(:,2);

res_beam = reshape(data-model,[],2);
resx = res_beam(:,1);
resy = res_beam(:,2);
pxy = reshape(fd.data,[],2);
clr = {[0, 0.4470, 0.7410],...
    [0.8500, 0.3250, 0.0980],...
    [0.4660, 0.6740, 0.1880]*0.8,...
    [0.9290, 0.6940, 0.1250]};

plts = {[1];[2];[3];,[2,3]};
pltlabs = {'0','45','90','all'};
legs = {{'0'},{'45'},{'90'},{'0','45','90'}};
for pltind = 4%1:length(plts)
    scale = 10;
    fig = figure(3);
    %fig.Position = [700 400 500 500];
    clf; hold on;
    
    pxy = reshape(fd.data,[],2);
    
    dks = plts{pltind};
    
    for dkind = 1:length(dks)
        chind = ismember(fd.sch,dks(dkind)) & flag'==3;
        quiver(pxy(chind,1),pxy(chind,2),resx(chind)*scale,resy(chind)*scale,0,'Color',clr{dkind})
        %plot(resx(chind),resy(chind),'b.')
    end
    grid on
    %legend({'DK=0','DK=45','DK=90'})
    legend(legs{pltind})
    title('Beam Center Best-Fit Residuals x10')
    xlabel('x (^o)')
    ylabel('y (^o)')
    xlim([-15 15])
    ylim([-15 15])
    figname = ['figs/' 'moonfit_residuals_quiver_zoomed_dk' pltlabs{pltind} '_beamfit'];
    saveas(gcf,figname,'png')
end



