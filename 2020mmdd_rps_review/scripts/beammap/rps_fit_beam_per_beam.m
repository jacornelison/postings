function bparam = rps_fit_beam_per_beam(sch, nsch, nrows, dirname, p, rpsopt,chans,DEBUG)
% function bparam = rps_fit_beam(sch,nsch,nrows,dirname,rpsopt,p,[chans],[DEBUG])
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
%	bparam	Data structure containing best-fit paramters and fit diagnostics.
%
%	Fields for bparam data structure:
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
% function bparam = rps_make_sims(sch,nsch,nrows,dirname,rpsopt,p,[chans],[DEBUG])

% Start by re-running reduc on all tods.
if 0    
    for i = 1:sch{nsch}.nrps
        fprintf('Reducing tod %i of %i\n.\n.\n.\n',i,sch{nsch}.nrps)
        rps_reduc(sch,nsch,i+(nrows-1)*sch{nsch}.nrps,'rps_data/2018/',p,rpsopt)
    end
end

if ~exist('DEBUG','var')
    DEBUG = 0;
end

sch = sch(nsch);
bparam = {};
prx = 2 * sind(p.r / 2) .* cosd(p.theta) * 180 / pi;
pry = 2 * sind(p.r / 2) .* sind(p.theta) * 180 / pi;



% Load in all timestreams for a given El offset in a schedule.
schname = sch{1}.name(end-15:end-4);
tods = {};
for j = 1:sch{1}.nrps
    fname = ['rps_data/2018/tods/tod_', schname, sprintf('_%03i.mat', j+(nrows-1)*sch{1}.nrps)];
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


bparam  = [];
[bparam.ch, bparam.err, bparam.param,bparam.gof, bparam.stat] = deal([]);
if DEBUG
    [bparam.model, bparam.data, bparam.x, bparam.y] = deal([]);
end


% Loop over tods
for i = 1:length(tods)
    % rpstod.x and rpstod.y are actually det-centered coords (x' and y')
    % Change this to boresight-centered coordinates (x and y)
    [r, theta, psi] = ...
        keck_beam_map_pointing(tods{j}.az-180, 90-tods{j}.el, -1*tods{j}.dk, ...
        rpsopt.mount,rpsopt.mirror, rpsopt.source, bs);
    tods{j}.x = 2 * sind(r / 2) .* cosd(theta) * 180 / pi;
    tods{j}.y = 2 * sind(r / 2) .* sind(theta) * 180 / pi;
    
    ch = channel;
    
    % Build best-guess parameter array.
    
    ch_ind = find(tods{i}.ch == ch);
    
    % Try instead fitting in just az and el.
    [A_guess, mind] = max(tods{i}.todcos(:,ch_ind));
    
    guess = [A_guess prx(ch) pry(ch) p.fwhm_maj(ch)^2 p.fwhm_maj(ch)^2 0 0];
    %guess = [A_guess tods{i}.az(mind) tods{i}.el(mind) p.fwhm_maj(ch)^2 p.fwhm_maj(ch)^2 0 0];
    
    % Trying matmin at Victor's suggestion
    % Found that it's infinitely better than fminsearch.
    % Build free parameter structure
    % [x0 y0 sigx sigy rho B0 psi xpol N1 N2 A]
    freepar.free = [1 1 1 1 1 1 1];
    freepar.lb = [-1e4 -1e3 -1e3 0 0 -1 -1000];
    freepar.ub = [1e4 1e3 1e3 2 2 1 1000];
    
    % Concat timestream data together
    %[ang_data,x,y] = rps_concat_data(tods,ch);
    x = tods{i}.x;
    y = tods{i}.y;
    data = tods{i}.todcos(:,ch_ind);
    
    % Very rough estimate of data rms noise by masking all data within 10sigma
    % of the beam, then calculate std.
    dist = sqrt((x-guess(2)).^2+(y-guess(3)).^2);
    ind = find(dist>mean(guess(4:5))*10);
    sd = std(data(ind));
    %mn = mean(ang_data(ind));
    % Create model timestream
    
    
    fprintf('Running sim %i of %i...\n',i,length(tods))
    
    % Estimate parameters
    [param, err, gof, stat, cov] = matmin('chisq',...
        guess, freepar,	'egauss2', data, sd, x, y);
    
    bparam.ch(i) = ch;
    bparam.param(i,:) = param;
    bparam.err(i,:) = err;
    bparam.gof(i) = gof;
    bparam.pte(i) = chi2cdf(gof,length(data)-length(guess));
    bparam.stat(i) = stat;
    bparam.dk(i) = tods{1}.scan.dk_ofs;
    bparam.data_rms(i) = sd;
    
    
    % debugging
    if DEBUG
        bparam.model = [bparam.model; egauss2(param,x,y)];
        bparam.data = [bparam.data; data];
        bparam.x = [bparam.x; x];
        bparam.y = [bparam.y; y];
    end
    
end


% Save the output file if we want
%  if exist('dirname','var') & ~isempty(dirname)
%    file_name = [dirname, 'params/param_', schname, sprintf('_%02i.mat',nrows)];
%	disp(['Saving file: ', file_name])
%    save(file_name,'bparam');
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


