function d = demod_choplet(d,nWave,lowPass,bandwidth,chopFreqLimits,keepChopWaveform,alternatePhase)
% function d = demod_choplet(d,nWave,lowPass,bandwidth,chopFreqLimits,keepChopWaveform,alternatePhase)
% 
% Demodulate using a filtered and windowed chop signal
% d is the time stream for instance from d=load_arc('arc',...);
% 
% Optional arguments:
%
% nWave : number of chopsignal wavelengths to window the chopsignal with (default 10)
%
% lowPass : lowPass frequency in Hz for final filtering (default: half the chop frequency)   
%
% bandwidth: bandwidth in Hz of the filter to render the step function of the chop to a sine-like wave (default: 2 Hz)
%
% chopFreqLimits: frequency range in which the reference signal is expected (default [5,15] Hz)
%
% keepChopWaveForm: default false, the reference waveform handed in over d.antenna0.pmac.fast_aux_input is not filtered when *true* is
% chosen. This is usefull when the modulated detector time streams don't look sine wave like such as in SkyNet test data. 
% Note that in this case the user must take care that provided reference corresponds to a wave with amplitude 1.
%
% alternatePhase: true (default), if true rotates the phase of the choplet by 90 deg. When the reference and the choped signal
% are fully in phase, all the signal will reside in the demodulated time series returned by this algorithm. If the phasing
% between the two is off the demodulated signal will be deminished (vanishes for 90 deg off phase). The signal portion that
% was lost can be accessed by demodulating with a reference/choplet whose phase is rotated by 90 deg (sin vs cos part of the
% other demodulators). Setting alternatePhase=true uses this alternative phasing to do the signal demodulation.
%
% See http://bicep.caltech.edu/~spuder/keck_analysis_logbook/analysis/20120310_ChopletDemodulator/ChopletDemodulator.html

if ~exist('nWave','var'), nWave=10; end
if ~exist('lowPass','var'), lowPass=0; end
if ~exist('bandwidth','var'), bandwidth=2; end
if ~exist('chopFreqLimits','var'), chopFreqLimits=[5,15]; end
if ~exist('keepChopWaveform','var'), keepChopWaveform=false; end
if ~exist('alternatePhase','var'), alternatePhase=1; end
  
% Trick to get the proper sampling frequency
%samprate = length(d.antenna0.time.utcfast)/((d.antenna0.time.utcfast(end)-d.antenna0.time.utcfast(1))*3600*24);
samprate = length(d.t)/((d.t(end)-d.t(1))*3600*24);
% Some preparation, semi-automatic search for the frequency of the chop:
chop = double(d.antenna0.pmac.fast_aux_input(:,1));
[chopSpec,cFreqs]=pwelch(chop,[],[],[],samprate);
% Expect the chop to be between 5 and 15 Hz
chopSpec(cFreqs<chopFreqLimits(1)|cFreqs>chopFreqLimits(2))=0;
[c,ind]=max(abs(chopSpec));
chopFreq=cFreqs(ind);
if ~keepChopWaveform
  % Apply filter around chopFreq with the bandwidth of xxx, eg 1Hz
  % The bandwidth must be sufficient to catch changes of the chop freq.
  % during datataking
  [b,a]=butter(3,[chopFreq-bandwidth/2,chopFreq+bandwidth/2]/(samprate/2));
  % The filtered chop signal is close to a sinewave but following the frequency
  % changes of the chop during datataking
  chop = filtfilt(b,a,chop);
end
  
if alternatePhase
  % Use the alternative phasing, this shifts the phase of the chop by 90deg:
  chop2 = imag(hilbert(chop));
end
  
% Number of points to catch nWave wavelengths of the chop:
% nWave should be a compromise:
% If too long the source will be too wide in the final beammap
% If too short the comparison between windowed chop and the signal becomes insensitive
windowPoints = (1./chopFreq * nWave)*samprate;
%  Round to the closest uneven number, to have a point symmetric window
windowPoints = 2.*round((windowPoints+1)/2)-1;
% The index of the middle of the window: window(windowOffset) != 1
windowOffset = round(windowPoints/2);
%  Make the actual window
window = hann(windowPoints);
% The window is normalized such that, when the chop waveform matches
% the data waveform, the amplitude of the data waveform is fully
% recovered: (effective bandwidth)  
window = window/(sum(window)/2);
% The wavelet will be the chop that is cut by the window 
% We'll use the envelope of the chop at the window center as normalization
if ~keepChopWaveform
  envelope = abs(hilbert(chop));
else
  envelope = ones(size(chop));
end
signal = d.mce0.data.fb;
% Here the correlation coefficients will be stored
nchan = size(signal,2);
crossCorrelation= zeros(length(chop),nchan);
crossCorrelation2 = zeros(length(chop2),nchan);
% Loop over the data, at the beginning and the end, where the window 
% cannot reach properly, zeros will remain in the crossCorrelation array
disp('demodulating...')
tic
for shift = 0:length(chop)-windowPoints
  % Cut the filtered chop signal with the window to do the wavetrain
  % and renormalize the wavelet amplitudes (this should be a small relative correction only):
  wavelet = chop(1+shift:shift+windowPoints).*window./envelope(shift+windowOffset);
  wavelet2 = chop2(1+shift:shift+windowPoints).*window./envelope(shift+windowOffset);
  % Do the cross correlation without the shifting (the shifting is done implicitly by the adaption
  % of the wavelet to the current chopFilt:
  % here the human readable version:
  % cc = sum(wavelet.*signal(1+shift:shift+windowPoints,:)); 
  % This is the same but slicing all channels at once:
  cc = sum(bsxfun(@times,wavelet,signal(1+shift:shift+windowPoints,:)));
  cc2 = sum(bsxfun(@times,wavelet2,signal(1+shift:shift+windowPoints,:)));
  crossCorrelation(shift+windowOffset,:) = cc;
  crossCorrelation2(shift+windowOffset,:) = cc2;
end
disp('...took:')
toc
 
% Under perfect conditions all frequency components of the filtered chop signal should now be
% removed in the crossCorrelation. However, if there was a phase shift between the chop and the
% signal, the chop frequency and higher orders of it will appear.
% Postulate that the signal we are interested in, should not be contained in the higher frequencies
% and low-pass filter the crossCorrelation:
if lowPass==0
  lowPass = chopFreq/2;
end
[b,a]=butter(3,lowPass/(samprate/2),'low');
crossCorrelationFilt = filtfilt(b,a,crossCorrelation);
crossCorrelationFilt2 = filtfilt(b,a,crossCorrelation2);
%d.mce0.data.fb=crossCorrelationFilt; 
%d.cos = crossCorrelationFilt;
%d.sin = crossCorrelationFilt2;
%d.ind = [1:length(d.cos)]';

% Downsample -- use just the indices at the peaks of the chop signal
[pks,loc] = findpeaks(chop);
d.cos = crossCorrelationFilt(loc,:);
d.sin = crossCorrelationFilt2(loc,:); % probably wrong...
d.ind = loc';

d.ref = d.ref(loc);
d.pos = d.pos(loc,:);
d.com = d.com(loc,:);
d.skips = d.skips(loc,:);
d.fb = d.fb(loc,:);
      
return

