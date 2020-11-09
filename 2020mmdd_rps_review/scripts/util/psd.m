function [freq,spec,phasearr]=psd(input, fsample)
%[freq,spec,phase]=psd(input, fsample)
%function to take the psd of a timestream and return also the freq axis and the phase. 
% input should be 1D
% fsample in Hz
 
input = double(input);
if floor(length(input) / 2) * 2 ~= length(input)
  input=input(2:end);
end
n = length(input) / 2;
transform = (fft(input));
transform=transform(1:n+1);
freq = (fsample / 2.0) * (linspace(0,n,n+1)) / n;
factor = cvec([1.0, repmat([2.0],1,n-1), 1.0]);
spec = sqrt(factor / freq(2)) .* transform/ (2*n);
%phasearr = phase(spec);
spec = abs(spec);

	
return
end


