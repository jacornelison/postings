function freqs=fftfreq(N, s)
% freqs = fftfreq(N, s)
%
% Returns the Discrete Fourier Transform sample frequencies, in cycles per
% unit. (Duplicated from numpy's fftfreq function.)
%
% INPUTS
%   N       Number of samples in the DFT.
%   s       Sample spacing in the original domain. Defaults to 1.0.
%
% OUTPUTS
%   freqs   Fourier axis frequencies of the N-point DFT with sample spacing
%           s. The elements have FFTW-compatible ordering; i.e. the positive
%           frequencies appear before the negative frequencies.
%
% EXAMPLE
%
%   t = 0:0.1:10;
%   x = sin(pi*t/2);
%   p = fft(x);
%   f = fftfreq(numel(t), 0.1);
%   hold on
%   stem(fftshift(f), fftshift(real(p)), 'b');
%   stem(fftshift(f), fftshift(imag(p)), 'r');
%

if ~exist('s', 'var') || isempty(s)
  s = 1.0;
end
del = 1.0 / (N * s);
if mod(N, 2) == 0
  f1 = floor(N / 2) - 1;
  f2 = -f1 - 1;
else
  f1 = (N - 1) / 2;
  f2 = -f1;
end
freqs = [0:f1, f2:-1] .* del;

end
