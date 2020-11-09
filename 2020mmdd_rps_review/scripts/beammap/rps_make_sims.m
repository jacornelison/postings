function rpssim = rps_make_sims(sch, nsch, nrows,dirname, in_param, p, rpsopt,chans,iter,DEBUG)
% function rpssim = rps_fit_beam(sch,nsch,nrows,dirname,rpsopt,p,[chans],[DEBUG])
% Estimates beam and polarization parameters from RPS data. This assumes that
% all RPS data in a schedule have been reduced into rpstod.mat files.
% Look into rps_reduc_driver if tods need to be reduced first.
%
% [ Arguments ]
%	sch		Cell array of metadata structures describing RPS schedules returned
%			by rps_log_read.m
%	nsch	sch index used for farming.
%	nrows	elevation offset index. Used for farming.
%	dirname	Directory to which the parameter data will be saved. Looks for a
%			params directory. Enter '' if you don't want to save the file.
%			(Future: makes one if it can't find one)
%	p		Focal plane info: p = get_array_info('yyyymmmdd','obs');
%	rpsopt	Data structure controlling options for RPS analysis
%	[chans]	Optional: Selects only channels specific channels if provided.
%			Channels are indexed at 1. I.E. gcp 0 == channel 1.
%	[DEBUG]	Optional: Returns extra variables for diagnostic use such as
%			data, x, y, model, and initial guess timestreams. File sizes are a
%			lot larger though. Wouldn't recommend saving.
%
%
% [ Output ]
%	rpssim	Data structure containing best-fit paramters and fit diagnostics.
%
%	Fields for rpssim data structure:
%		ch:			Channel list common to all TOD's in an elevation offset.
%		ang_param:	Array of polarization response parameters
%		ang_err: 	1-sigma error returned by matmin.
%		ang_status:	Flag indicating fit status. 3 is convergence.
%		ang_gof:	Goodness-of-fit value returned by .
%		ang_cov:	Covariance matrix returned by matmin.
%		ang_var:	Residual variance between best-fit model and data.
%		dk:			Deck angle for the observation schedule.
%		data_rms:	Estimated rms of timestream data.
%
%	[Optional outputs]
%
%		data:		Concatenated amplitude timestream data across all RPS angles
%		model:		Modeled amplitude timestream data
%		x:			Concatenated x position data
%		y:			Concatenated y position data
%		guess:		First-guess modeled timestream data that's fed to matmin.
%
% Last update: 2018-Jul-02 JAC
%
%
% function rpssim = rps_make_sims(sch,nsch,nrows,dirname,rpsopt,p,[chans],[DEBUG])

if ~exist('DEBUG','var')
    DEBUG = 0;
end

sch = sch(nsch);
rpssim = {};
prx = 2 * sind(p.r / 2) .* cosd(p.theta) * 180 / pi;
pry = 2 * sind(p.r / 2) .* sind(p.theta) * 180 / pi;



% Load in all timestreams for a given El offset in a schedule.
schname = sch{1}.name(end-15:end-4);
tods = {};
for j = 1:sch{1}.nrps
    fname = [dirname 'tods/tod_', schname, sprintf('_%03i.mat', j+(nrows-1)*sch{1}.nrps)];
    try
        load(fname);
        % Sometimes the tod is empty. Skip it.
        if isfield(rpstod,'az')
            tods{end+1} = rpstod; % Allows us to skip faulty timestreams if needed.
        end
    catch ME
        if strcmp(ME.identifier,'MATLAB:load:couldNotReadFile')
            fprintf('Could not load: %s\nSkipping...\n',fname)
        else
            disp(sprintf('Caught exception: %s\nWe''ll try to continue anyway.',ME.identifier))
        end
    end
    
    disp(['loaded ' fname '...'])
end
if isempty(tods)
    fprintf('No data in TODs!\n')
    return
end
dirname = 'rps_data/2018/';
%param_name = [dirname, 'params/param_', schname, sprintf('_%02i.mat',nrows)];
%load(param_name);

% Convert to x/y focal plane coordinates by using the boresight parameters
% r and theta which are identically zero.
bs.r = 0;
bs.theta = 0;
bs.chi = 0;
bs.chi_thetaref = 0;
bs.drumangle = 0;
mount.aperture_offr = 0; % Locate boresight 'aperture' in center of drum.

% Check if any specific channels are wanted.
if exist('chans','var')
    channel = chans;
else
    channel = tods{1}.ch;
    for i = 2:length(tods)
        channel = intersect(channel,tods{i}.ch);
    end
end

% Check if TOD's have any channels at all.
if isempty(channel)
    fprintf('Row %02i doesn''t have any channels. Skipping...\n',nrows)
    return
end
fprintf('Fitting a total of %i channels. This could take awhile.\n', length(channel))

% rpstod.x and rpstod.y are actually det-centered coords (x' and y')
% Change this to boresight-centered coordinates (x and y)
for j = 1:length(tods)
    [r, theta, psi] = ...
        keck_beam_map_pointing(tods{j}.az-180, 90-tods{j}.el, -1*tods{j}.dk, ...
        rpsopt.mount,rpsopt.mirror, rpsopt.source, bs);
    tods{j}.x = 2 * sind(r / 2) .* cosd(theta) * 180 / pi;
    tods{j}.y = 2 * sind(r / 2) .* sind(theta) * 180 / pi;
end

% Loop over channels
for i = 1
    
    ch = channel;
    
    % Build best-guess parameter array.
    % Func form: F_gauss(x,y)*(cos(2*(th+psi))-(eps+1)/(eps-1))*(N!*cos(th+N2)+1)
    % [x0 y0 sigx sigy rho B0 psi xpol N1 N2 A]
    %
    %
    % Guess psi based off Colin's keck_beam_map_pointing model
    %ang_guess = atand(tand(mean(tods{1}.psi(:,tods{1}.ch==ch))));
    %A_guess = ones(1,length(tods));
    for j=1:length(tods)
        ch_ind = find(tods{j}.ch == ch);
        
        A_guess(j) = max(tods{j}.todcos(:,ch_ind));
        rot(j)  = tods{j}.scan.abs_rot-rpsopt.grid.horizontal;
    end
    
    
    
    
    %Generate amp params from input values
    initparam = [prx(ch) pry(ch) p.fwhm_maj(ch)^2 p.fwhm_maj(ch)^2 0 rps_get_mod_model(in_param,rot)];
    ang_model_0 = rps_get_model(initparam,tods,ch,rpsopt);    
    
    guess = initparam;
    
    % Concat timestream data together
    [data,x,y] = rps_concat_data(tods,ch);
    
    % Very rough estimate of data rms noise by masking all data within 10sigma
    % of the beam, then calculate std.
    dist = sqrt((x-guess(1)).^2+(y-guess(2)).^2);
    ind = find(dist>mean(guess(3:4))*10);
    sd = std(data(ind));
    clear data
    %%%%%
    % Create model timestream
    

    
    % Trying matmin at Victor's suggestion
    % Found that it's infinitely better than fminsearch.
    % Build free parameter structure
    % [x0 y0 sigx sigy rho A1 ... AN]
    freepar.free = ones(size(guess));
    freepar.lb = [-30 -30 0 0 -1 repmat(-1e4,size(A_guess))];
    freepar.ub = [30 30 2 2 1 repmat(1e4,size(A_guess))];
    
    % loop over iterations
    rpssim  = [];
    [rpssim.ch, rpssim.berr, rpssim.bparam,rpssim.bgof, rpssim.bstat] = deal([]);
    [rpssim.aerr, rpssim.aparam,rpssim.agof, rpssim.astat, rpssim.amodel] = deal([]);
    [rpssim.rot] = deal([]);
    
    for k = 1:iter
        fprintf('Running sim %i of %i...\n',k,iter)
        % Add White noise noise estimate to model
        wnoise = normrnd(0,sd,size(ang_model_0));
        ang_data = ang_model_0+wnoise;
        
        % If wanted, add amp-dependent noise
        anoise = normrnd(0,0.05,size(ang_model_0)).*ang_data;
        ang_data = ang_data+anoise;
        % Estimate parameters
        [bparam, berr, bgof, bstat, ang_cov] = matmin('chisq',...
            guess, freepar,	'rps_get_model',ang_data, sd, tods, ch, rpsopt);
        
        rpssim.ch(k) = ch;
        rpssim.bparam(k,:) = bparam;
        rpssim.berr(k,:) = berr;
        rpssim.bgof(k) = bgof;
        rpssim.bstat(k) = bstat;
        
        %% Extract best fit model and calculate residuals
        %ang_model = rps_get_model(ang_param,tods,ch,rpsopt);
        %res = ang_model-ang_data;
        %
        %rpssim{end}.ang_var = nanvar(res);
        rpssim.dk(k) = tods{1}.scan.dk_ofs;
        rpssim.data_rms(k) = sd;
        % debugging
        if DEBUG & iter == 1
            rpssim.model = rps_get_model(bparam,tods,ch,rpsopt);
            rpssim.data = ang_data;
            rpssim.x = x;
            rpssim.y = y;
            rpssim.wnoise = wnoise;
            rpssim.anoise = anoise;
        end
        
        %%
        % Fit modulation curve.
        
        % Recycle parameter bounds
        % [psi eps N1 N2 A]
        afreepar.free = [1 1 1 1 1];
        afreepar.lb = [-360 -1 -1e4 -1e4 -1e4];
        afreepar.ub = [360 1 1e4 1e4 1e4];
        
        % First guess.
        ang_guess = atand(tand(mean(tods{1}.psi(:,tods{1}.ch==ch))));
        aguess = [ang_guess 0 0 0 max(bparam(6:end))/2];
        
        % Estimate parameters
        [aparam, aerr, agof, astat, acov] = matmin('chisq',...
            aguess, afreepar,'rps_get_mod_model',bparam(6:end),berr(6:end),rot);
        
        rpssim.aparam(k,:) = aparam;
        rpssim.aerr(k,:) = aerr;
        rpssim.agof(k) = agof;
        rpssim.apte(k) = 1-chi2cdf(agof,length(rot)-length(aguess));
        rpssim.astat(k) = astat;
        %rpssim.acov = acov;
        rpssim.amodel(k,:) = rps_get_mod_model(aparam,rot);
        rpssim.rot(k,:) = rot;
        
        if DEBUG
            rpssim.aguess = aguess;
        end
        
    end
    
end
% Save the output file if we want
%  if exist('dirname','var') & ~isempty(dirname)
%    file_name = [dirname, 'params/param_', schname, sprintf('_%02i.mat',nrows)];
%	disp(['Saving file: ', file_name])
%    save(file_name,'rpssim');
%  end


disp('Complete!')


% Subfunction: rps_concat_data
% ---
% Loops through source angles and
% concats data for this channel.
function [data,x,y] = rps_concat_data(tods,ch)

[data,x,y] = deal([]);

for i=1:length(tods)
    ch_ind = find(tods{i}.ch == ch);
    data = [data; tods{i}.todcos(:,ch_ind)];
    x = [x; tods{i}.x];
    y = [y; tods{i}.y];
end


