function d = demod_lockin(data,chop,time,Fc,slope,chopFreqLimits,varargin)
% function d = demod_lockin(data,chop,time,Fc,slope,chopFreqLimits,varargin)
%
% Demodulates data using a basic lock-in amplification routine. Still in
% development. Speak to James Cornelison for details.
%
% [Arguments]
%   d       Timestream data. d = load_arc() configuration expected.
%
%   Fc     Roll-off frequency. Related to lock-in timeconstant by
%          tau=1/(2 pi Fc) (default: one-quarter of the chop frequency)
%
%   slope   Rolloff of butterworth filter in db/Octave. Typically multiples
%           of 6. (default: 24 db/Octave)
%
%   chopFreqLimits  Frequency range in which the reference signal is
%           expected. If scalar, overrides chop freq search.
%           (default [5,15] Hz)
%
%   [Optional Args]
%   'chopShift', double
%           Shifts the reference signal by a (possibly non-integer)
%           value in samples. Useful if you want to get all the signal in X
%           so you're not rectifying the noise by using R. (default: 0)
%
%   'keepChopWaveForm'
%           The reference waveform from
%           d.antenna0.pmac.fast_aux_input is not filtered when *true* is
%           chosen. This is useful when the modulated detector time
%           streams don't look sine wave like such as in SkyNet test data.
%           Note that in this case the user must take care that provided
%           reference corresponds to a wave with amplitude 1.
%
%   'DEBUG'   Outputs diagnostic information to the console.
%
%   'expand'    Don't downsample after low-pass filtering.
%   [Example]
%       d = demod_lockin(d,chopref,t,5,24,[19,21],'chopShift',1.4,'DEBUG','expand')
%
% function d = demod_lockin(data,chop,time,Fc,slope,chopFreqLimits,varargin)

if (~exist('slope','var') | isempty(slope)), slope=24; end
if ~exist('chopFreqLimits','var'), chopFreqLimits=[5,15]; end

% initialize optional args
keepChopWaveform=0;
DEBUG = false;
chopShift = 0;
doDownSamp = 1;
PLOT = false;

for i = 1:length(varargin)
    if strcmp(varargin{i},'keepChopWaveForm')
        keepChopWaveform=1;
    end
    
    if strcmp(varargin{i},'DEBUG')
        DEBUG = 1;
    end
    
    if strcmp(varargin{i},'chopShift')
        chopShift = varargin{i+1};
    end
    
    if strcmp(varargin{i},'expand')
        doDownSamp = 0;;
    end
    
    if strcmp(varargin{i},'setLowPass')
        lpf = varargin{i+1};
    end
    
    if strcmp(varargin{i},'PLOT')
        DEBUG = true;
        PLOT = varargin{i+1};
    end
    
end
if DEBUG
    tic


end



% If the scan crosses midnight, we'll get a weird number for out sample
% rate.
dt = diff(time);
sampFreq = 1/mean(dt(dt>0));



% If chop limit is scalar, override chop-freq search
if length(chopFreqLimits)==1
    chopFreq = chopFreqLimits;
elseif length(chopFreqLimits)==2
    if 1%DEBUG
        fprintf('\nLooking for reference frequency...\n')
    end
    [chopSpec,cFreqs]=pwelch(chop,[],[],[],sampFreq);
    chopSpec(cFreqs<chopFreqLimits(1)|cFreqs>chopFreqLimits(2))=0;
    [c,ind]=max(abs(chopSpec));
    chopFreq=cFreqs(ind);
    if chopFreq > chopFreqLimits(1) & chopFreq < chopFreqLimits(2)
        fprintf('Found reference frequency at: %2.2f Hz\n',chopFreq);
    else
        fprintf('Found reference frequency at: %2.2f Hz\n',chopFreq);
        fprintf('Which is outside limits of %2.2f and %2.2f Hz.\n',chopFreqLimits(1),chopFreqLimits(2));
        fprintf('Setting reference frequency to %2.2f Hz.\n',mean(chopFreqLimits));
        chopFreq = mean(chopFreqLimits);
    end
else
    error('Chop limits should either be size 2 vector or scalar')
end

% If no roll-off freq was provided, make it 50% of the chop frequency
if (~exist('Fc','var') | isempty(Fc))
    Fc = chopFreq/2.0;
end
if (~exist('lpf','var') | isempty(lpf))
    lpf = Fc;
end

    if PLOT

        fig = figure(1);
clf(fig)
set(fig,'Position',[500 -200 950 450])
set(fig,'PaperPosition',[-2 2.5 [12.6 6]*0.5])
hold off
subplot(2,1,1)
periodogram(chop,[],[],sampFreq)
hold on
plot([0 0]+(chopFreq+Fc/2),[-1 1]*200,'r')
plot([0 0]+(chopFreq-Fc/2),[-1 1]*200,'r')
ylim([-100 30])

    end

%keyboard()
if ~keepChopWaveform
    if DEBUG
        fprintf('\nFiltering reference signal...\n')
    end
    % Apply filter around chopFreq with the bandwidth of the lowpass filter bandwidth, eg 1Hz
    % The bandwidth must be sufficient to catch changes of the chop freq.
    % during datataking
    %[b,a]=butter(3,[chopFreq-bandwidth/2,chopFreq+bandwidth/2]/(sampFreq/2));
    [b,a]=butter(4,[chopFreq-Fc/2,chopFreq+Fc/2]/(sampFreq/2));
    % The filtered chop signal is close to a sinewave but following the frequency
    % changes of the chop during datataking?
    chop = filtfilt(b,a,chop);
    
    %normalize
    %chop = chop/sinefit(time,chop,chopFreq);
end

    if PLOT

subplot(2,1,2)
hold off
periodogram(chop,[],[],sampFreq)
hold on
plot([0 0]+(chopFreq+Fc/2),[-1 1]*200,'r')
plot([0 0]+(chopFreq-Fc/2),[-1 1]*200,'r')
%ylim([-100 30])

    end


if chopShift ~= 0
    chop = fraccircshift(chop,chopShift);
end

% Rotate chop by 90deg.
chop2 = imag(hilbert(chop));


% Get number of channels
nchan = size(data,2);

% Mix reference with input signal
mixSig = data.*repmat(chop,1,nchan);
mixSig90 = data.*repmat(chop2,1,nchan);

[mixSigFilt, mixSigFilt90] = deal(zeros(size(mixSig)));

if slope==0

    mixSigFilt = mixSig;
    mixSigFilt90 = mixSig90;

else
    nfilt = floor(slope/6);
    [b,a]=butter(nfilt,lpf/(sampFreq/2),'low');
    for i = 1:size(mixSig,2)
    mixSigFilt(:,i) = filtfilt(b,a,mixSig(:,i));
    mixSigFilt90(:,i) = filtfilt(b,a,mixSig90(:,i));
    end
end

% Try downsampling to the new Nyquist Freq value
if doDownSamp
    ds = round(sampFreq/2/lpf);
else
    ds = 1;
end


d.X = downsample(mixSigFilt,ds);
d.Y = downsample(mixSigFilt90,ds);
d.R = 2*sqrt(d.X.^2+d.Y.^2);
d.theta = atan2(d.Y,d.X)*180/pi;
d.t = downsample(time,ds);
d.ind = (downsample(1:length(mixSigFilt),ds))';

% We get ringing on the order of 20 samples (at 20 Hz chop) at the ends of the
% timestream. Easiest thing to do is to truncate.
ind = false(size(d.ind));
ind(20:end-20) = true;

d.X = d.X(ind,:);
d.Y = d.Y(ind,:);
d.R = d.R(ind,:);
d.theta = d.theta(ind,:);
d.t = d.t(ind);
d.ind = d.ind(ind);


%keyboard()
if DEBUG
    toc
end


