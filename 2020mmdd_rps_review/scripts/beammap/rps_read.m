function rpstod = rps_read(scan, ch, p, rpsopt)
% Updated by JAC
% rpstod = rps_read(scan, ch, p, rpsopt)
%
% Reads time-ordered data for one RPS scan from arc files and demodulates.
%
% [Arguments]
%   scan    RPS scan data structure, returned by rps_log_read.m
%   ch      List of channels, indexed from 1 and running across all MCEs,
%           to read. For example, ch = [529, 530] would select the first
%           two detectors from rx1.
%   p       Detector parameter structure from get_array_info.m. This
%           should be the full p structure, even though we will downselect
%           based on the specified channels.
%   rpsopt  Data structure controlling options for RPS analysis.
%
%   Fields for rpsopt data structure:
%     rpsopt.mount   Mount model data structure. Default value supplied by
%                    keck_beam_map_pointing.m
%     rpsopt.mirror  Mirror model data structure. Default value supplied
%                    by keck_beam_map_pointing.m
%     rpsopt.source  Source position data structure. Default value supplied
%                    by keck_beam_map_pointing.m
%     rpsopt.flip_dk Flips the sign of dk encoder offset for BICEP2 RPS
%                    schedules that use mirror online pointing model.
%     rpsopt.fbit    Specify the value of the feature bit that indicates
%                    scanning motion. Defaults to 3 (i.e. f0+f1).
%     rpsopt.do_deconv       Set to true if you want deconvolution applied
%                            to bolometer timestreams (default), false
%                            otherwise. Deconvolution code doesn't seem to
%                            work if rpsopt.rx is set to something other
%                            than all rx.
%     rpsopt.chop_channel    Specify which pmac fast_aux_input is used to
%                            record the chop reference. Defaults to 1.
%     rpsopt.chop_threshold  Specify the analog level used to differentiate
%                            between chop reference high/low states. Default
%                            value is 1250, which works well for Keck.
%     rpsopt.chop_shift      Number of samples to shift the chop reference.
%                            Defaults to zero, but you will want to set this
%                            to a value that maximizes signal in the cosine
%                            demodulation channel.
%     rpsopt.chop_deglitch   Set to true to remove single sample
%                            transitions in the chop reference, which lead
%                            to glitchy looking spikes in the demodulated
%                            timestreams. False by default.
%
% [Returns]
%   rpstod  Data structure containing demodulated TOD for selected channels.
%
%   Fields for rpstod data structure:
%     rpstod.scan    Record of the scan data structure argument.
%     rpstod.ch      Record of the channel list argument.
%     rpstod.rpsopt  Record of the rpsopt argument.
%     rpstod.p       Detector parameter data structure, cut down to include
%                    only selected channels.
%     rpstod.todcos  Cos demodulated time-ordered data for all detectors.
%     rpstod.todsin  Sin demodulated time-ordered data for all detectors.
%                    Note that all of the pointing arrays have been
%                    downsampled to the locations of the cos demodulated
%                    timestream.
%     rpstod.az      Azimuth coordinates for telescope boresight.
%     rpstod.el      Elevation coordinates for telescope boresight.
%     rpstod.pa      Dk angle coordinate for telescope.
%     rpstod.r       r coordinate of the source for all detectors.
%     rpstod.theta   theta coordinate of the source for all detectors.
%     rpstod.psi     psi coordinate of the detector co-polar axes projected
%                    onto the source pointing.

% Last update: 2018-02-01 JAC
disp('Initializing variables')
% Fill out rpsopt data structure with default parameters, if necessary.
if ~exist('rpsopt', 'var')
    rpsopt = [];
end
% Use default mount model from keck_beam_map_pointing.m
if ~isfield(rpsopt, 'mount')
    rpsopt.mount = [];
end
% Use default mirror model from keck_beam_map_pointing.m
if ~isfield(rpsopt, 'mirror')
    rpsopt.mirror = [];
end
% Use default source location from keck_beam_map_pointing.m
if ~isfield(rpsopt, 'source')
    rpsopt.source = [];
end
% For BICEP2, set the flip_dk keyword to 1 so that it correctly
% handles the online pointing model that was used for RPS observations.
if ~isfield(rpsopt, 'flip_dk')
    rpsopt.flip_dk = 0;
end
% By default, scan feature bit is 3 (=f0+f1).
if ~isfield(rpsopt, 'fbit')
    rpsopt.fbit = 3;
end
% Deconvolve timestreams.
if ~isfield(rpsopt, 'do_deconv')
    rpsopt.do_deconv = 0;
end
% Which fast_aux_input channel is used for the chop reference.
if ~isfield(rpsopt, 'chop_channel')
    % Usual value for Keck.
    rpsopt.chop_channel = 1;
    % Other useful values:
    % BICEP2: chop_channel = 4
end
% Threshold somewhere in between chop reference high and low state.
if ~isfield(rpsopt, 'chop_threshold')
    % Chose value for Nov 2012 Keck RPS data.
    rpsopt.chop_threshold = 1250;
end
% Can shift the chop reference by integer values.
if ~isfield(rpsopt, 'chop_shift')
    rpsopt.chop_shift = 0;
end
% By default, don't deglitch chop reference.
if ~isfield(rpsopt, 'chop_deglitch')
    rpsopt.chop_deglitch = false;
end
% By default, Choose 'hires' square demod.
if ~isfield(rpsopt, 'demodtype')
    rpsopt.demodtype = 'hires';
end

% Store input arguments in output data structure.
rpstod.scan = scan;
rpstod.ch = ch;
rpstod.rpsopt = rpsopt;

% Exit if no channels were hit.
if length(ch) == 0
    return
end

%% Setup register list
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

%% Load ARC
disp('Loading arc files. This might take a minute.')
% Read arc files, selected channels only.
list_arc_files('arc', mjd2datestr(scan(1).t1), mjd2datestr(scan(1).t2));
data = load_arc('arc', mjd2datestr(scan(1).t1), mjd2datestr(scan(1).t2), reg);
% If flip_dk keyword is set, flip the sign for the dk axis parts of
% encoder_off.
if rpsopt.flip_dk
    data.antenna0.tracker.encoder_off(:,3) = ...
        -1 * data.antenna0.tracker.encoder_off(:,3);
end

% Downselect p structure based on channel list.
if isfield(p, 'expt')
    p = rmfield(p, 'expt');
end
p = structcut(p, ch);
% Add to results data structure.
rpstod.p = p;

% Combine all mce into mce0 field of the data structure.
if ~isfield(data, 'mce0')
    data.mce0 = [];
end
for m = {'mce1', 'mce2', 'mce3', 'mce4'};
    if isfield(data, m)
        data.mce0 = structcat(2, [data.mce0, getfield(data, m{1})]);
        data = rmfield(data, m);
    end
end

%% Apply Inverse Pointing Model
% Schedule encoder_zeros -52.621, -7.0, 2.626
% Normal encoder_zeros 127.37899, -82.41876, -2.626 # 2016 Feb 2 - BICEP3 from 20160122
% Schedule encoder_cals: 2304000, -2304000, -574400, ant=all # az,el,dk counts per turn, reverse el&dk direction for ffflat
% Normal encoder_cals: 2304000, 2304000, 574400

data.antenna0.tracker.encoder_mul = repmat([2304000, 2304000, 574400],...
    size(data.antenna0.tracker.encoder_mul,1),1);
data.antenna0.tracker.encoder_off = repmat(([127.37899, -82.41876, -2.626])*3.6e6,...
    size(data.antenna0.tracker.encoder_off,1),1);
data.antenna0.tracker.encoder_sign = abs(data.antenna0.tracker.encoder_sign);

%data.antenna0.tracker.encoder_mul(1,:) = [2304000, -2304000, -574400];
%data.antenna0.tracker.encoder_off = repmat(([-52.621,-7.0, 2.626])*3.6e6,...
%    size(data.antenna0.tracker.encoder_off,1),1);


% Apply offline pointing model.
if strcmp(expt,'bicep3')
    pm = get_pointing_model(scan.t1,0, data);
    pm.az_tilt_ha = 0;
    pm.az_tilt_lat = 0;
    pm.el_tilt = 0;
    data = invpointing_model(data, pm);
else
    pm = get_pointing_model(scan.t1, 0, data);
    data = invpointing_model(data, pm);
end

% Select only data with the correct feature bit set.
[mask_slow, mask_fast] = select_feature_bit(data, rpsopt.fbit);

% Also get rid of skipped frames. They should be the same for all channels.

ind = find(sum(data.mce0.data.fb == 0,2)==size(data.mce0.data.fb,2));
mask_fast(ind) = false;

if ~any(mask_fast == 1)
    fprintf('No data in selected channels!\n')
    rpstod.flag = 'No data in selected channels';
    return
end

% Calculate downsampled encoder positions.
az = data.pointing.hor.az(mask_fast);
el = data.pointing.hor.el(mask_fast);
dk = data.pointing.hor.dk(mask_fast);



%% Deconvolve Channels
if rpsopt.do_deconv
    % Deconvolve data.
    % Just one continuous block for deconvolution.
    dc.s = min(find(mask_slow));
    dc.e = max(find(mask_slow));
    dc.sf = min(find(mask_fast));
    dc.ef = max(find(mask_fast));
    
    if ~isfield(rpsopt,'xferopts')
        fprintf('No xferopts detected, using legacy mode...\n');
        rpsopt.xferopts.legacy = true;
    elseif rpsopt.xferopts.legacy==false & ~isfield(rpsopt.xferopts,'delay')
        fprintf('No convolution delay detected. Defaulting to 0 ms delay...\n');
        rpsopt.xferopts.delay = 0;
    elseif rpsopt.xferopts.legacy==false
        fprintf('Using %i ms delay in deconvolution...\n',rpsopt.xferopts.delay);
    end
    % Calculate deconvolution kernel.
    data = get_xferfunc(data, dc, 0, 0,rpsopt.xferopts);
    %  data = get_xferfunc(data, dc, 0, 0);
    % Work-around for a bug when running get_xferfunc with low-pass
    % filter off.
    for ii=1:length(data.tf)
        data.tf(ii).lpf = 1;
    end
    % Apply deconvolution.
    data = deconv_scans(data, p);
end

% Get chop reference.
chopref = data.antenna0.pmac.fast_aux_input(mask_fast, rpsopt.chop_channel);
chopref = chopref > mean(chopref);
% Shift chop reference.
%chopref = circshift(chopref, rpsopt.chop_shift);

% Deglitch chop reference.
if rpsopt.chop_deglitch
    chopref = chop_deglitch(chopref);
end
%keyboard()
% Demodulate data
disp('Demodulating data...')
switch rpsopt.demodtype
    case 'hires'
        disp('Demod type: square+hires')
        chopref = shift_sqw(chopref,rpsopt.chop_shift);
        chopout = cleanchopref(chopref);
        
        [todcos, todsin, indcos, indsin] = ...
            square_demod(data.mce0.data.fb(mask_fast,:), chopout.refclean_hires, 'highres');
    case 'norm'
        disp('Demod type: square')
        chopref = shift_sqw(chopref,rpsopt.chop_shift);
        [todcos, todsin, indcos, indsin] = ...
            square_demod(data.mce0.data.fb(mask_fast,:), chopref);
        
    case 'lockin'
        disp('Demod type: Lock-In')
        d = demod_lockin(data.mce0.data.fb(mask_fast,:),chopref,...
            data.antenna0.time.utcfast(mask_fast,2),...
            1,60,[19,21],'setLowPass',9);
        todcos = d.X;
        todsin = d.Y;
        indcos = d.ind;
        indsin = d.ind;
end

% Store results.
rpstod.todcos = todcos;
rpstod.todsin = todsin;
scos = size(todcos);
ssin = size(todsin);
if size(todcos) == size(todsin)
    rpstod.todquad = sqrt(todcos.^2+todsin.^2);
elseif scos(1) < ssin(1)
    rpstod.todquad = sqrt(todcos.^2+todsin(1:end-1,:).^2);
elseif scos(1) > ssin(1);
    rpstod.todquad = sqrt(todsin.^2+todcos(1:end-1,:).^2);
    rpstod.todcos = rpstod.todcos(1:end-1,:);
    indcos = indcos(1:end-1,:);
end


rpstod.az = az(indcos);

rpstod.el = el(indcos);

rpstod.dk = dk(indcos);



% Calculate pointing for detectors relative to source.
% This calculates the source position in x'/y' coords. We really don't need
% it so I'm disabling it for now.
if 0
    if strcmp(expt, 'bicep3')
        [r, theta, psi] = ...
            keck_beam_map_pointing(rpstod.az, rpstod.el, rpstod.dk, rpsopt.mount, ...
            rpsopt.mirror, rpsopt.source, rpstod.p);
    else
        [r, theta, psi] = ...
            keck_beam_map_pointing(rpstod.az, rpstod.el, rpstod.dk, rpsopt.mount, ...
            rpsopt.mirror, rpsopt.source, rpstod.p);
        
    end
    
    
    rpstod.r = r;
    rpstod.theta = theta;
    rpstod.psi = psi;
    rpstod.x = 2. * sind(rpstod.r/2).*cosd(rpstod.theta)*180.0/pi;
    rpstod.y = 2. * sind(rpstod.r/2).*sind(rpstod.theta)*180.0/pi;
end
rpstod.p.mce = p.mce;
rpstod.flag = 'Complete';


