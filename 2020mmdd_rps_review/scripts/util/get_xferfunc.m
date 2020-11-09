function d=get_xferfunc(d,dc,combine,uselpf,xferopts)
% d=get_xferfunc(d,dc,combine,uselpf,xferopts)
%
% input:
%   dc = continuous deconv period
%   combine = explicitly state combine on/off (1/0) [optional]
%           = default is on (1)
%   uselpf = apply a low pass filter when deconvolving
%          = default is 1, which applies 5Hz lpf for reduc_initial
%          = can be a user-defined freq-domain lpf
%          = if 0, then no lpf is applied
%          Note that LPF design assumes archived data at 20 samples/s.
%          For combine off data or other archive rate, pass band will
%          change while filter coefficients stay the same.
%
%   xferopts = a structure of options modifying the behavior; see below.
%
% output:
%   d.tf(nrx,ndc) where ndc = num deconv blocks, nrx = num MCEs
%   tf.deconvkern = deconvolution kernel to be conv with data
%                   inverse of filter functions convolved with
%                   low pass filter to keep deconv from ringing
%                   designed to be FIR
%   tf.kernel = impulse response of convolved mce/gcp filters
%               an array or matrix of kernel values
%               if a matrix includes measured tf
%   tf.seconds = array of length(tf.kernel) in seconds
%   tf.dc = deconvolution blocks to pass to deconv_scans.m
%   tf.lpf = impulse response function of low pass filter
%   tf.info = filter info
%   tf.gcp_taps = antialiasing filter applied by GCP
%   tf.gcpds = downsample factor applied in GCP
%
% xferopts: Any of the following options may be set.
%
%   .compact
%     Defaults to true. If false, the full ~10s deconvolution kernel is
%     used.
%
%     NB!! Setting to false will cause a single NaN in a timestream to grow to
%     cover ~10 sec!
%
%   .delay
%     Defaults to []. If empty, no constant time (linear phase) delay will
%     be applied. If true, aux_data will be searched for a delay specification.
%     Otherwise, may be a vector of time delays in milliseconds to apply for
%     each MCE.
%
%   .legacy
%     Defaults to true. If true, the original behavior used for B2, and
%     most of Keck's era will be used. See
%       http://bicep.rc.fas.harvard.edu/bkcmb/analysis_logbook/analysis/20180507_deconv/
%     for full details on what changed.
%
% calculate bolometer impulse response
% determine appropriate mce/gcp filter kernel from d structure
% empirical transfer functions (with mce/gcp filters deconvolved) can be added in optionally
% kernel is automatically truncated, starting with a length of 10sec
% convolve inverse of filter kernel with low pass filter to create time domain deconv kernel

% for season 1 bicep2 data you can't tell whether combine is on/off
% some data may be taken with combine on, n=1 which would apply the 
% fir filter, but not downsample.  For season 2 bicep2/keck data the
% fir filter taps are archived, so you can determine combine on/off
% assume that data was acquired with combine on if nsnap>1, and
% combine off if nsnap=1, unless you explicitly specify combine=0
if ~exist('combine','var')
  combine=[];
end

% default should be to apply lpf, as used in reduc_initial.m
if ~exist('uselpf','var') || isempty(uselpf)
  uselpf=1;
end

if ~exist('xferopts','var') || isempty(xferopts)
  xferopts = struct();
end
if ~isfield(xferopts, 'compact')
  xferopts.compact = true;
end
if ~isfield(xferopts, 'delay')
  xferopts.delay = [];
end
if ~isfield(xferopts, 'legacy')
  % Legacy implementation used up until 2018. See
  %   http://bicep.rc.fas.harvard.edu/bkcmb/analysis_logbook/analysis/20180507_deconv/
  % for description of original behavior and why it was changed.
  xferopts.legacy = true;
end

nmce=size(d.mce0.frame.status,2);

ndelay = length(xferopts.delay);
if isnumeric(xferopts.delay) && length(xferopts.delay) == 1
  xferopts.delay = repmat(xferopts.delay, nmce, 1);
elseif islogical(xferopts.delay)
  if xferopts.delay
    yyyymmdd = str2double(mjd2datestr(get_date(d), 'yyyymmdd'));
    auxdelay = find_file_by_date(yyyymmdd, 'aux_data/timeconst/delay');
    if isempty(auxdelay)
      xferopts.delay = [];
    else
      xferopts.delay = auxdelay.delay;
    end
  else
    xferopts.delay = [];
  end
end
if ~isempty(xferopts.delay) && length(xferopts.delay) ~= nmce
  error('Mismatch in %d delays for %d MCEs', length(xferopts.delay), nmce);
end

tf=[];
for ii=1:length(dc.sf)

  % Loop over MCEs
  for jj=1:nmce

    % store the deconv block for deconv_scans
    tf(jj,ii).dc.sf=dc.sf(ii);
    tf(jj,ii).dc.ef=dc.ef(ii);
    tf(jj,ii).dc.s=dc.s(ii);
    tf(jj,ii).dc.e=dc.e(ii);

    % determine the mce settings
    data_rate=double(nanmedian(d.mce0.cc.data_rate(dc.s(ii):dc.e(ii),jj)));
    num_rows=double(nanmedian(d.mce0.cc.num_rows(dc.s(ii):dc.e(ii),jj)));
    row_len=double(nanmedian(d.mce0.cc.row_len(dc.s(ii):dc.e(ii),jj)));

    % determine the sampling rate before gcp downsampling is applied
    fsamp=50e6/(num_rows*row_len*data_rate);

    % determine gcp downsample factor
    gcpds=double(nanmedian(d.array.frame.nsnap(dc.s(ii):dc.e(ii))));

    % Check if MCE is off
    if (data_rate==0 && num_rows==0 && row_len==0)
      tf(jj,ii).info{1}=['mce' num2str(jj-1) ' is off.'];
      continue
    end

    % Make npts some multiple of the total downsample factor to keep impulse
    % functions uniformly sampled.
    % Start with kernel of length ~10s.
    kernlen=10;
    dstot=data_rate*gcpds;
    % Points in MCE impulse response function (at full internal rate).
    nptmce = ceil(fsamp * data_rate * kernlen);
    nptmce = nptmce - mod(nptmce, dstot);
    % Number of points in gcp impulse response (at non-combined rate).
    nptgcp = nptmce / data_rate;

    % calculate mce filter transfer function
    ds = ifelse(xferopts.legacy, dstot, data_rate);
    [hmce,imce,fltr_coeff]=get_mce_xfer(d,dc,ii,jj,ds,nptmce);
    if false
      figure(); set(gcf(), 'Name', 'MCE');
      plot_zfunc(hmce, fsamp / gcpds * dstot/ds, dstot/ds);
    end

    % calculate gcp filter transfer function
    ds = ifelse(xferopts.legacy, gcpds, 1);
    [hgcp,igcp,gcp_taps]=get_gcp_xfer(d,dc,ii,fsamp,ds,nptgcp,combine);
    if false
      figure(); set(gcf(), 'Name', 'GCP');
      plot_zfunc(hgcp, fsamp / ds, gcpds / ds);
    end

    % include measured bolometer transfer functions
    % the mce/gcp transfer functions must already be deconv
    % add this in later, for now just assume it's flat response
    hmeas=ones(length(hgcp),1);

    % convolve transfer functions
    htot=hgcp.*hmce.*hmeas;
    if false
      figure(); set(gcf(), 'Name', 'MCE * GCP');
      plot_zfunc(htot, fsamp / ds, gcpds / ds);
    end

    % back to the time domain
    imptot=ifft(htot,'symmetric');

    if ~xferopts.legacy
      % decimate by gcp downsample factor in the time domain
      imptot = downsample(imptot, gcpds);
    end

    % make filter kernel as compact as possible; keep elements greater than
    % approx rounding error (~4.5 * eps(1.0), to be specific).
    if xferopts.compact
      imptot = imptot(abs(imptot) > 1e-15);
    end

    % normalize kernel to one at dc, do it after compactifying
    % this is the combined impulse response function of the 
    % mce & gcp filters
    kern=imptot./abs(nansum(imptot));

    % Build the low-pass filter (if appropriate)
    [lpf,ilpf,lpfinfo] = get_lpf(uselpf, length(kern));
    N=length(lpf)-1;
    if false
      figure(); set(gcf(), 'Name', 'Low-pass filter');
      plot_zfunc(lpf, fsamp/gcpds);
    end

    % take normalized, truncated kernel back to f domain
    % zero pad it to the same length as the lpf
    % but shift right, pad, shift left to preserve phase
    % shift point closest to zero to end for padding
    [a,shftind]=min(abs(kern));
    shft=length(kern)-shftind;
    impp=circshift([circshift(kern,shft);zeros(N+1-length(kern),1)],-shft);
    impf=fft(impp);

    % The arbitrary delay is inherently infinite impulse response, so we
    % construct it here to match the size of the full kernel.
    if ~isempty(xferopts.delay)
      axf = cvec(2*pi * fftfreq(N+1, 1 / (fsamp / gcpds)));
      hdelay = exp(-1i * axf * xferopts.delay(jj) / 1000);
      impf = impf .* hdelay;
      if false
        figure(); set(gcf(), 'Name', 'Delay');
        plot_zfunc(hdelay, fsamp / gcpds);
        figure(); set(gcf(), 'Name', 'MCE * GCP * Delay')
        plot_zfunc(impf, fsamp / gcpds);
      end
    end

    % convolve lpf with inverse of mce/gcp filter transfer functions
    % and normalize to 1 at DC
    % this is the final time domain deconv kernel to be convolved
    % with data in deconv_scans
    deconvkern=fftshift(ifft(lpf./impf,'symmetric'));
    deconvkern=deconvkern./abs(nansum(deconvkern,1));
    if false
      figure(); set(gcf(), 'Name', 'Deconv');
      plot_zfunc(lpf./impf, fsamp/gcpds);
    end

    % replace tf.gcp structure with tf.kernel matrix
    % if no empirical transfer function is included,
    % tf.kernel will be 1 channel wide
    tf(jj,ii).kernel=kern;
    tf(jj,ii).lpf=ilpf;
    tf(jj,ii).deconvkern=deconvkern;
    tf(jj,ii).sec_kern=(0:length(kern)-1)'./(fsamp/gcpds);
    tf(jj,ii).sec_deconvkern=(0:length(deconvkern)-1)'./(fsamp/gcpds);
    tf(jj,ii).gcp_taps=gcp_taps;
    tf(jj,ii).gcpds=gcpds;
    tf(jj,ii).legacy=xferopts.legacy;

    tf(jj,ii).info{1}=['mce filter coefficients = ' int2str(fltr_coeff)];
    tf(jj,ii).info{2}=[int2str(length(gcp_taps)) ' gcp filter taps = ' num2str(gcp_taps,3)];
    tf(jj,ii).info{3}=['final data rate = ' num2str(fsamp/gcpds,6) ' Hz'];
    tf(jj,ii).info{4}=lpfinfo;
    if ~isempty(xferopts.delay)
      tf(jj,ii).info{5} = sprintf(['with phase term for time delay of ' ...
            '%0.3f milliseconds'], xferopts.delay(jj));
    end
  end

end

% tack tf structure on
d.tf=tf;

end

function val=ifelse(cond,trueval,falseval)
% val=ifelse(cond, trueval, falseval)
%
% Simulate C's ternary operator,
%   val = cond ? trueval : falseval

  if cond
    val = trueval;
  else
    val = falseval;
  end
end

%%%%%%%%%%%%%%%%%%%%%%

function [zfunc,ifunc,fltr_coeff]=get_mce_xfer(d,dc,ii,jj,dstot,nptmce)

% default (discrete) delta response - i.e. no filter
ifunc = [1; zeros(nptmce/dstot - 1, 1)];
zfunc = ones(nptmce/dstot, 1);
fltr_coeff = 0;

% determine mce settings
% apply mce filter only if data_mode is filtered feedback
data_mode=nanmedian(d.mce0.rc1.data_mode(dc.s(ii):dc.e(ii),jj));

% if the data_mode is not for filtered data just
% return a flat mce filter transfer function
if ~any(data_mode == [2,6:10])
  return
end

% check to see if mce type2 filter parameters are being used
% for bicep2 type2 filter is the default, but for keck this can be programmed.
% for type2 mce filter, the filter cutoff scales as frame_rate/400;
if isfield(d.mce0.rc1,'fltr_type')
  fltr_type=int32(nanmedian(d.mce0.rc1.fltr_type(dc.s(ii):dc.e(ii),jj)));

  % Check whether this MCE is on.
  if (fltr_type==0)
    warning('quad_analysis:get_xferfunc:mce_off', 'mce%d is off.', jj-1);
    return
  end

  fltr_coeff=double(nanmedian(d.mce0.rc1.fltr_coeff(dc.s(ii):dc.e(ii),:),1));
  % Get coefficients for this MCE only
  nmces=size(d.mce0.rc1.data_mode,2);
  nparams=size(fltr_coeff,2)/nmces;
  fltr_coeff=reshape(fltr_coeff, nparams, nmces)';
  fltr_coeff=fltr_coeff(jj,:);

else % no fltr_type field
  % mce type2 filter coeff, valid for all sampling rates
  fltr_coeff=double([32295, 15915, 32568, 16188, 3, 14]);
end

% first pass at correct mce filter transfer function
% try to include filter coefficient quantization effects
% filters quantized to 15bit resolution
% mce filter is realized as 2 cascaded biquads, with
% coeff specified for the 2 second order sections (SOS)
sosmatrix=[1 2 1 1 -fltr_coeff(1)/2^14 fltr_coeff(2)/2^14;...
           1 2 1 1 -fltr_coeff(3)/2^14 fltr_coeff(4)/2^14];

% generate the impulse response functions at the appropriate sampling rates
% [impmce,timemce]=impz(bmce,amce,nptmce,fsamp.*data_rate);
impmce=sosfilt(sosmatrix,[1; zeros(nptmce-1,1)]);

% downsample mce impulse response function to the gcp sample rate
% right now, do it dumbly and don't introduce phase offset
ifunc = downsample(impmce, dstot);
zfunc = fft(ifunc);
% ensure normalization
ifunc = ifunc ./ zfunc(1);
zfunc = zfunc ./ zfunc(1);

end

%%%%%%%%%%%%%%%%%%%%%%

function [zfunc,ifunc,taps]=get_gcp_xfer(d,dc,ii,fsamp,gcpds,nptgcp,combine)

% determine if combine is on/off if possible (need d.array.filter)
% if gcpds>1 you know combine is on and a filter is applied
% if gcpds=1 combine can be on or off, if you can't determine it 
% from d.array.filter and no combine input is specified, assume on
if gcpds==1 
  if isfield(d.array,'filter')
    tapsum=nansum(d.array.filter.taps(dc.sf(ii):dc.ef(ii)));
    combine = tapsum ~= 0;
  else
    % Is this case actually consistent with the comment above? Shouldn't
    % this not exist and be allowed to fall through to the next elseif?
    combine = false;
  end
elseif isempty(combine)
  combine = true;
end

% if combine is off, return a flat gcp filter transfer function
if ~combine
  zfunc = ones(nptgcp/gcpds, 1);
  ifunc = [1; zeros(nptgcp/gcpds - 1, 1)];
  taps = 0;
  return
end

% if d.array.filter exists, read in filter taps
% bicep2 first season didn't record taps so we have
% to determine the applied filter from the date here
% what precision is actually used in gcp filter calcs?
if isfield(d.array,'filter')
  filtstart=find(d.array.filter.idx(dc.sf(ii):dc.ef(ii))==0);
  ntaps=nanmedian(diff(filtstart));
  filtindx=ceil(length(filtstart)/2)+1;
  taps=[];
  while length(taps)<ntaps
    disp(['Have filtindx=' num2str(filtindx) ', dropping by one to ' num2str(filtindx-1) '.']);
    filtindx=filtindx-1;
    taps=transpose(d.array.filter.taps((dc.sf(ii)+filtstart(filtindx))-1:(dc.sf(ii)+filtstart(filtindx+1))-2));
  end
else
  % determine bicep2 first year configuration from the date
  % d.t is mjd, convert to gregorian
  dat=str2num(datestr(datenum('nov-17-1858','mmm-dd-yyyy')+median(d.antenna0.time.utcfast),'yyyymmddHH'));
  if(dat>2010032519 & dat<2011031015)
    % starting tag: 20100325B01_dk265  ending tag:  20110309C10_dk068
    % starting date:  25-Mar-2010:20:15:53  ending date:  10-Mar-2011-14:35:46
    % fir filter fir-n-60-100-x5
    % 5x gcp downsample
    N=120;
    fp=0.12;
    fs=0.2;
    taps=firpm(N,[0 fp fs 1], [1 1 0 0]);
    taps=taps/sum(taps);
  end
  if(dat>2010101 & dat<2010032514)
    % ending tag 20100322I10_dk130
    % ending date:  25-Mar-2010:13:36:41  
    % fir filter fir-80-90-x5
    % 5x gcp downsample
    N=50;
    fp=0.16;
    fs=0.18;
    taps=fir2(N,[0 fp fs 1], [1 1 0 0]);
  end
end

% generate the impulse response functions at the appropriate sampling rates
impgcp=impz(taps,1,nptgcp,fsamp);

% the gcp impulse response function is linear-phase, to make it
% zero-phase shift it left by (length(bgcp)-1)/2
shft=(length(taps)-1)/2;
impgcpzp=circshift(impgcp,-shft);

% downsample gcp impulse response function to 
% the final downsampled gcp sample rate.
% as long as the downsampling of the impulse response function is
% even, this is the correct thing to do.
ifunc = downsample(impgcpzp, gcpds);
zfunc = fft(ifunc);
% ensure normalization
ifunc = ifunc ./ zfunc(1);
zfunc = zfunc ./ zfunc(1);

end

%%%%%%%%%%%%%%%%%%%%%%

function [zfunc,ifunc,info]=get_lpf(uselpf,kernlen)
  if uselpf == 0
    info = 'No low-pass filter';
    ifunc = [1; zeros(kernlen-1, 1)];
    zfunc = ones(kernlen, 1);
    return
  end

  if length(uselpf) > 1
    % Keep user-supplied low-pass filter.
    info = 'User-supplied low-pass filter';
    zfunc = uselpf ./ uselpf(1);
    ifunc = ifft(zfunc, 'symmetric');
    return
  end

  % design low pass filter (3-5Hz stop band) to be convolved with
  % inverse of filter kernel, this is what is used in deconv_scans
  % low pass filter should be FIR to make invalid regions finite
  % N=32 gives the minimum for all filter combinations
  % coefficients found from chisq min of ripple < 3Hz
  N=32;
  if N < kernlen
    N=2*ceil(kernlen/2);
  end
  fp=0.286069192737340;
  fs=0.5; % 0.928499208688734;
  hlpf=firpm(N,[0 fp fs 1],[1 1 0 0])';
  % make it zero-phase
  ifunc=circshift(hlpf,-N/2);
  % go to the f domain and normalize
  zfunc=fft(ifunc);
  ifunc=ifunc./zfunc(1);
  zfunc=zfunc./zfunc(1);

  info = sprintf('%0.3f, ', hlpf');
  info = info(1:end-2); % remove trailing ', '
  info = sprintf(['low-pass filter firpm, fp = %0.4f, fs = %0.4f, with %d ' ...
                  'taps = [ %s ]'], fp, fs, length(hlpf), info);
end

function yyyymmdd=get_date(d)
  tval = d.antenna0.time.utcslow;
  yyyymmdd = tval(find(tval(:,1) ~= 0), :);
  if size(yyyymmdd, 2) > 1
    yyyymmdd = yyyymmdd(1) + yyyymmdd(2) / 86400;
  end
end

%%%%%%%%%%%%%%%%%%%%%%

% useful for interactive debug plotting

function plot_zfunc(z, fsamp, downsamp, debounce)
  if ~exist('downsamp','var') || isempty(downsamp)
    downsamp = 1;
  end
  if ~exist('debounce','var') || isempty(debounce)
    debounce = true;
  end

  N = length(z);
  dsamp = fsamp / downsamp;

  f = fsamp * [0:(N-1)]' ./ N;
  % truncate at the Nyquist
  z = z(1:ceil(end/2));
  f = f(1:ceil(end/2));

  xlim_full = @() xlim([f(2), f(end)]);
  xlim_zoom = @() xlim([f(2), dsamp/2]);

  clf(); setwinsize(gcf(), 1000, 800);
  set(gcf(), 'Renderer', 'painters')
  set(gcf(), 'DefaultAxesXMinorTick', 'on')
  set(gcf(), 'DefaultAxesYMinorTick', 'on')

  subplot(3,2,1); cla(); box on
    xlim_full(); ylim([-0.05, 1.05])
    hold on; set(gca(), 'xscale', 'log')
    plot(xlim(), [0,0], 'k:')
    plot(xlim(), [1,1], 'k:')
    plot(dsamp/2 .* [1,1], ylim(), 'k--')
    semilogx(f, abs(z));
    title('Amplitude response')
    ylabel('power fraction')
    %xlabel('frequency [Hz]')
    set(gca(), 'XMinorTick','on', 'XMinorGrid','on')
  if downsamp > 1
    subplot(3,2,2); cla(); box on
      xlim_zoom(); ylim([0.995, 1.005])
      hold on; set(gca(), 'xscale', 'log')
      plot(xlim(), [1,1], 'k:')
      semilogx(f, abs(z));
      title('Amplitude response (zoomed)')
      ylabel('power fraction')
      xlabel('frequency [Hz]')
      set(gca(), 'XMinorTick','on', 'XMinorGrid','on')
  end

  % Clean up the low-amplitude portion of the response function a bit;
  % mainly helps with visual impression of phase and delays
  ztmp = z;
  if debounce
    blk = find_blk(abs(real(ztmp))<1e-3 & abs(imag(ztmp))<1e-3);
    if length(blk.e) == 1 && blk.e == length(ztmp)
      tmp = angle(z(blk.s-1));
      %ztmp(blk.s:blk.e) = eps(abs(z(blk.s-1)))/sqrt(2) .* exp(1i*tmp);
      if abs(tmp) < 1e-6
        ztmp(blk.s:blk.e) = 0;
      else
        ztmp(blk.s:blk.e) = NaN;
      end
    end
  end

  ph = round(unwrap(angle(ztmp)) ./ (2*pi) * 1e10) ./ 1e10;
  if downsamp > 1
    ee = 1 + floor(length(f) / downsamp); % where f ~ dsamp
  else
    ee = find(abs(ztmp) > 1e-2, 1, 'last');
  end
  subplot(3,2,3); cla(); box on
    xlim_full();
    if ph(end) < 0
      ylim([ph(end), -0.1 * ph(end)])
    elseif min(ph(1:ee)) < max(ph(1:ee))
      ylim([min(ph(1:ee)), max(ph(1:ee))])
    else
      ylim(nanmedian(ph(1:ee)) + [-1,1])
    end
    hold on; set(gca(), 'xscale', 'log')
    plot(dsamp/2 .* [1,1], ylim(), 'k--')
    semilogx(f, ph);
    title('Phase response')
    ylabel('Angle [2\pi\times{}rad]')
    %xlabel('frequency [Hz]')
    set(gca(), 'XMinorTick','on', 'XMinorGrid','on')
  if downsamp > 1
    subplot(3,2,4); cla(); box on
      xlim_zoom();
      if ph(ee) < 0
        ylim([ph(ee), -0.1 * ph(ee)])
      elseif min(ph(1:ee)) < max(ph(1:ee))
        ylim([min(ph(1:ee)), max(ph(1:ee))])
      else
        ylim(nanmedian(ph(1:ee)) + [-1,1])
      end
      hold on; set(gca(), 'xscale', 'log')
      semilogx(f, ph);
      title('Phase response (zoomed)')
      ylabel('Angle [2\pi\times{}rad]')
      xlabel('frequency [Hz]')
      set(gca(), 'XMinorTick','on', 'XMinorGrid','on')
  end

  gd_msec = 1e3 .* groupdelay(ztmp, fsamp, N);
  min_msec = min(0, min(gd_msec(1:ee)));
  max_msec = max(gd_msec(1:ee)); % min/max delay of response
  max_samp = 1e3 / dsamp; % decimation sample spacing
  pl_msec = 1.05 * max(max_msec, max_samp); % max delay w/ buffer for plots
  subplot(3,2,5); cla(); box off
    xlim_full(); ylim([min_msec, pl_msec])
    hold on; set(gca(), 'xscale', 'log')
    plot(dsamp/2 .* [1,1], ylim(), 'k--')
    plot(xlim(), [0,0], 'k:')
    plot(xlim(), max_samp .* [1,1], 'k:')
    semilogx(f, gd_msec);
    xlabel('frequency [Hz]')
    title('Group delay')
    set(gca(), 'XMinorTick','on', 'XMinorGrid','on')
    y1 = gca();
    set(get(y1, 'YLabel'), 'String', 'Delay [msec]')
    y2 = axes('YAxisLocation', 'right', 'XAxisLocation','top', 'Color','none');
    set(y2, 'Position', get(y1, 'Position'));
    set(y2, 'XTick', get(y1, 'XTick'), 'XTickLabel',{});
    set(y2, 'YLim', dsamp/1e3 .* [min_msec, pl_msec]);
    set(get(y2, 'YLabel'), 'String', '[decimated samples]')

  pd_msec = 1e3 .* phasedelay(ztmp, fsamp, N);
  min_msec = min(0, min(pd_msec(1:ee)));
  max_msec = max(pd_msec(1:ee)); % max delay of response
  max_samp = 1e3 / dsamp; % decimation sample spacing
  pl_msec = 1.05 * max(max_msec, max_samp); % max delay w/ buffer for plots
  subplot(3,2,6); cla(); box off
    xlim_full(); ylim([min_msec, pl_msec])
    hold on; set(gca(), 'xscale', 'log')
    plot(dsamp/2 .* [1,1], ylim(), 'k--')
    plot(xlim(), [0,0], 'k:')
    plot(xlim(), max_samp .* [1,1], 'k:')
    semilogx(f, pd_msec);
    xlabel('frequency [Hz]')
    title('Phase delay')
    set(gca(), 'XMinorTick','on', 'XMinorGrid','on')
    y1 = gca();
    set(get(y1, 'YLabel'), 'String', 'Delay [msec]')
    y2 = axes('YAxisLocation', 'right', 'XAxisLocation','top', 'Color','none');
    set(y2, 'Position', get(y1, 'Position'));
    set(y2, 'XTick', get(y1, 'XTick'), 'XTickLabel',{});
    set(y2, 'YLim', dsamp/1e3 .* [min_msec, pl_msec]);
    set(get(y2, 'YLabel'), 'String', '[decimated samples]')
end

function dly=phasedelay(H, samprate, N)
  % N = full-length Fourier vector; should be 2*length(H) + [0,1] if H is
  %     transfer function truncated to the even/odd Nyquist frequency
  if ~exist('N','var') || isempty(N)
    N = length(H);
  end
  ph = unwrap(angle(H(:)));
  fq = 2*pi*samprate * [0:(length(H)-1)]' ./ N;
  dly = -ph ./ fq;
end

function dly=groupdelay(H, samprate, N)
  % N = full-length Fourier vector; should be 2*length(H) + [0,1] if H is
  %     transfer function truncated to the even/odd Nyquist frequency
  if ~exist('N','var') || isempty(N)
    N = length(H);
  end
  ph = unwrap(angle(H(:)));
  dfq = 2*pi*samprate * 1 / N;
  dly = -diff(ph) ./ dfq;
  % assuming the DC phase delta is the same as the next
  dly = [dly(1); dly];
end
