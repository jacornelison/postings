function [ i_param ] = rps_fit_beam(sch, nsch, nrows,p,rpsopt,dirname,chans,DEBUG)
% function [ i_param ] = rps_fit_beam(sch, nsch, nrows,p,rpsopt,dirname,chans,DEBUG)
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
%   [azel]  Optional: Fits for beam centers in az/el/dk first and then
%           converts centers to x/y/phi. For prototyping. Default On.
%
%
% [ Output ]
%	i_param	Data structure containing best-fit paramters and fit diagnostics.
%
%	Fields for i_param data structure:
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
%       phi:        Concatenated phi angle data
%		guess:		First-guess modeled timestream data that's fed to matmin.
%
% Last update: 2018-Jul-02 JAC
%
%
% function [ i_param ] = rps_fit_beam(sch, nsch, nrows,p,rpsopt,dirname,chans,DEBUG)

if ~exist('DEBUG','var')
    DEBUG = 0;
end

if ~isfield(rpsopt,'fitopt')
    rpsopt.fitopt.type == 'norm';
end

if ~isfield(rpsopt,'grid')
    rpsopt.grid.horizontal = 0;
elseif ~isfield(rpsopt.grid,'horizontal')
    rpsopt.grid.horizontal = 0;
end


sch = sch(nsch);
i_param = {};
prx = 2 * sind(p.r / 2) .* cosd(p.theta) * 180 / pi;
pry = 2 * sind(p.r / 2) .* sind(p.theta) * 180 / pi;

% Loop over Elevation offsets.
for ii = nrows
    
    % Load in all timestreams for a given El offset in a schedule.
    schname = sch{1}.name(end-15:end-4);
    tods = {};
    for j = 1:sch{1}.nrps
        fname = [dirname 'tods/tod_', schname, sprintf('_%03i.mat', j+(ii-1)*sch{1}.nrps)];
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
    if length(tods)==0
        fprintf('No data in TODs!\n')
        continue
    end
    
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
        
        % Find any channel among TOD's.
        for i = 2:length(tods)
            channel = union(channel,tods{i}.ch);
        end
        ch_count = zeros(size(channel));
        
        % Cut out channels that have less than 75% of rasters
        for i = 1:length(tods)
            ch_count = ch_count + ismember(channel,tods{i}.ch);
        end
        channel = channel(ch_count>=length(tods)*0.75);
        
    end
    
    % Add NaN timestreams to the channels that pass, but don't have
    % everything.
    %keyboard()
    for i = 1:length(tods)
        for j = 1:length(channel)
            if ~any(tods{i}.ch==channel(j))
                dummydata = NaN(size(tods{i}.az));
                tods{i}.ch(end+1) = channel(j);
                tods{i}.todcos(:,end+1) = NaN(size(tods{i}.todcos(:,1)));
                tods{i}.todsin(:,end+1) = NaN(size(tods{i}.todsin(:,1)));
                tods{i}.todquad(:,end+1) = NaN(size(tods{i}.todquad(:,1)));
            end
        end
    end
    
    % Check if TOD's have any channels at all.
    if isempty(channel)
        fprintf('Row %02i doesn''t have any channels. Skipping...\n',ii)
        continue
    end
    fprintf('Fitting a total of %i channels. This could take awhile.\n', length(channel))
    
    % rpstod.x and rpstod.y are actually det-centered coords (x' and y')
    % Change this to boresight-centered coordinates (x and y)
    
    rot = 1:length(tods);
    for j = 1:length(tods)
        [r, theta, psi] = ...
            keck_beam_map_pointing(tods{j}.az, tods{j}.el, tods{j}.dk, ...
            rpsopt.mount,rpsopt.mirror, rpsopt.source, bs);
        
        tods{j}.x = 2 * sind(r / 2) .* cosd(theta) * 180 / pi;
        tods{j}.y = 2 * sind(r / 2) .* sind(theta) * 180 / pi;
        tods{j}.phi = psi;
        %Grab source angles as well.
        rot(j)  = tods{j}.scan.abs_rot-rpsopt.grid.horizontal;
    end
    
    % Loop over channels
    for i = 1:length(channel)
        fprintf('Fitting channel index %i of %i channels...\n',i,length(channel))
        ch = channel(i);
        
        i_param{end+1}  = [];
        i_param{end}.ch = ch;
        
        % Build best-guess parameter array.
        % Func form: F_gauss(x,y)*(cos(2*(th+psi))-(eps+1)/(eps-1))*(N!*cos(th+N2)+1)
        % [x0 y0 sigx sigy rho B0 psi xpol N1 N2 A]
        %
        %
        % Guess psi based off Colin's keck_beam_map_pointing model
        % Concat timestream data together
        [data,x,y,phi,az,el,dk,quad] = rps_concat_data(tods,ch);
        
        %keyboard()
        
        % Need to be able to fit any number of rastersets.
        [A_guess, A_guess_i] = deal(ones(1,length(tods)));
        
        for j=1:length(tods)
            % Griddata takes longer, but is more reliable.
            %[A_guess(j), A_guess_i(j)] = max(tods{j}.todcos(:,tods{j}.ch == ch));
            A_guess(j) = griddata(tods{j}.x,tods{j}.y,tods{j}.todcos(:,tods{j}.ch==ch),prx(ch),pry(ch));
        end
        
        % Guess for beams
        guess = [prx(ch) pry(ch) p.fwhm_maj(ch)^2 p.fwhm_maj(ch)^2 0 A_guess];
        
        
        % Trying matmin at Victor's suggestion
        % Found that it's infinitely better than fminsearch.
        % Build free parameter structure
        % [x0 y0 sigx sigy rho A1 ... AN]
        freepar.free = [1 1 1 1 1 ~isnan(A_guess)];
        freepar.lb = [-30 -30 0 0 -1 repmat(-1e4,size(A_guess))];
        freepar.ub = [30 30 2 2 1 repmat(1e4,size(A_guess))];
        
        % Very rough estimate of data rms noise by masking all data within 10sigma
        % of the beam, then calculate std.
        dist = sqrt((x-guess(1)).^2+(y-guess(2)).^2);
        %ind = find(dist>mean(guess(3:4))*10);
        sd = nanstd(data((dist>mean(guess(3:4))*10)));
        
        % Estimate parameters
        [bparam, berr, gof, stat, cov] = matmin('chisq',...
            guess, freepar,	'rps_get_model',data, sd, tods, ch, rpsopt);
        
        
        % Interpolate to find phi of RPS at estimated beam center
        phi_s = griddata(x,y,phi,bparam(1),bparam(2));
        i_param{end}.phi_s = phi_s;
        
        i_param{end}.bparam = bparam;
        i_param{end}.berr = berr;
        i_param{end}.bchi2 = gof;
        i_param{end}.bpte = 1-chi2cdf(gof,length(data)-length(bparam));
        i_param{end}.bstat = stat;
        i_param{end}.bcov = cov;
        
        % Extract best fit model and calculate residuals
        model = rps_get_model(bparam,tods,ch,rpsopt);
        res = model-data;
        
        i_param{end}.var = nanvar(res);
        i_param{end}.dk = tods{1}.scan.dk_ofs;
        i_param{end}.data_rms = sd;
        
        % debugging
        if DEBUG
            i_param{end}.model = model;
            i_param{end}.data = data;
            i_param{end}.x = x;
            i_param{end}.y = y;
            i_param{end}.phi = phi;
            i_param{end}.guess = rps_get_model(guess,tods,ch,rpsopt);
        end
        
        % Fit modulation curve.
        
        if strcmp(rpsopt.fitopt.type,'legacy')
            % Recycle parameter bounds
            % [psi eps N1 N2 A]
            freepar.free = [1 1 1 1 1];
            freepar.lb = [-360 -1 -1e4 -1e4 -1e4];
            freepar.ub = [360 1 1e4 1e4 1e4];
            
            % First guess.
            ang_guess = atand(tand(p.chi(ch)+p.chi_thetaref(ch)-phi_s));
            
            guess = [ang_guess 0.003 0 0 max(bparam(6:end))/2];
            
            
            % Estimate parameters
            [aparam, aerr, agof, astat, acov] = matmin('chisq',...
                guess, freepar,	'rps_get_mod_model',bparam(6:end),berr(6:end),rot);
            
            % Phi of detector co-polar axis is pol response angle - phi_s:
            % Why is this minus? Because of the mirror?
            amodel = rps_get_mod_model(aparam,rot);
            phi_d = aparam(1)-phi_s;
        elseif strcmp(rpsopt.fitopt.type,'norm')
            % Grab the pointing and orientation vectors based on best fit
            % beam centers.
            % Create Boresight pointing, orientation, and aperture position
            
            [A0,B3_0,B1_0] = kbmp_mount(az,el,dk,rpsopt.mount,0);
            
            % Reflect with mirror
            [mpos mnorm] = kbmp_mirror(az,el,rpsopt.mount,rpsopt.mirror);
            [A0,B3_0,B1_0] = kbmp_reflect(A0, B3_0, B1_0, mpos, mnorm);
            
            [A B3 B1] = deal(NaN(1,3));
            for k = 1:3
                A(1,k) = griddata(x,y,A0(:,k),bparam(1),bparam(2));
                B3(1,k) = griddata(x,y,B3_0(:,k),bparam(1),bparam(2));
                B1(1,k) = griddata(x,y,B1_0(:,k),bparam(1),bparam(2));
            end
                       
            % Use ideal pol angle as initial guess
            ang_guess = atand(tand(p.chi(ch)+p.chi_thetaref(ch)));
            ag = 1;
            
            xpol_guess = 0;
            xp = 1;
            
            % RPS mispointing is in local RPS az/el.
            if isfield(rpsopt.fitopt,'align')
                aln_guess = rpsopt.fitopt.align;
                al = [0 0];
            else
                aln_guess = [0 0];
                al = [1 1];
            end
            if isfield(rpsopt.fitopt,'nut')
                nut_guess = rpsopt.fitopt.nut;
                nt = [0 0];
            end
            
            gain_guess = max(bparam(6:end))/2;
            gn = 1;
            
            guess = [ang_guess xpol_guess nut_guess aln_guess gain_guess];
            
            freepar.free = [ag xp nt al gn];
            freepar.lb = [-10+ang_guess -1 -15 -15 -15 -15 0];
            freepar.ub = [10+ang_guess 1 15 15 15 15 1e6];
            
            % Estimate parameters
            [aparam, aerr, agof, astat, acov] = matmin('chisq',...
                guess, freepar,	'rps_get_mod_model_vectors',bparam(6:end),berr(6:end),rot,A,B3,B1,rpsopt.source);
            
            amodel = rps_get_mod_model_vectors(aparam,rot,A,B3,B1,rpsopt.source);
            
            phi_d = aparam(1);
        end
        

        i_param{end}.aparam = aparam;
        i_param{end}.phi_d = phi_d;
        i_param{end}.aerr = aerr;
        i_param{end}.agof = agof;
        i_param{end}.apte = 1-chi2cdf(agof,length(rot)-length(guess));
        i_param{end}.astat = astat;
        i_param{end}.acov = acov;
        i_param{end}.amodel = amodel;
        i_param{end}.rot = rot;
        i_param{end}.quad = quad;
        if DEBUG
            i_param{end}.aguess = guess;
        end
    end
    % Save the output file if we want
    if exist('dirname','var') & ~isempty(dirname) & ~DEBUG
        file_name = [dirname, 'params/param_', schname, sprintf('_%02i.mat',ii)];
        disp(['Saving file: ', file_name])
        save(file_name,'i_param');
    end
end

disp('Complete!')


% Subfunction: rps_concat_data
% ---
% Loops through source angles and
% concats data for this channel.
function [data,x,y,phi,az,el,dk,quad] = rps_concat_data(tods,ch)

[data,x,y,phi,az,el,dk,quad] = deal([]);

for i=1:length(tods)
    ch_ind = find(tods{i}.ch == ch);
    data = [data; tods{i}.todcos(:,ch_ind)];
    quad = [quad; tods{i}.todsin(:,ch_ind)];
    x = [x; tods{i}.x];
    y = [y; tods{i}.y];
    phi = [phi; tods{i}.phi];
    az = [az; tods{i}.az];
    el = [el; tods{i}.el];
    dk = [dk; tods{i}.dk];
    
end


