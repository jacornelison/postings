function sn = bm_addsn(signal,noise,snr)
% sn = bm_addsn(signal,noise,snr)
%
% Add a signal and noise beam map with the A/B/S/D struct for residbeam
% work.  By itself can simply add two maps together, or can scale the
% noise map by some factor snr
%
% INPUTS
%        signal: standard beam struct (fields A/B/S/D with 'i','f')
%        noise:  standard beam struct
%        snr:    factor by which to scale 'noise' map (can be empty)
%                If present, scales relative to peak of 'signal' map
%                i.e. noise_new = noise*max(signal)/snr
%
% OUTPUTS
%        sn:     signal+noise beam map

if ~isempty(snr)
  maxsig = (nanmax(real(signal.Ai(:))) + nanmax(real(signal.Bi(:))))/2;
  scalefac = maxsig/snr;
  noise.Ai = noise.Ai.*scalefac;
  noise.Bi = noise.Bi.*scalefac;
end

sn.Ai = signal.Ai + noise.Ai;
sn.Bi = signal.Bi + noise.Bi;
sn.Si = (sn.Ai + sn.Bi)/2;
sn.Di = (sn.Ai - sn.Bi)/2;            

return