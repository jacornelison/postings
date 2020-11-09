function r=rms(x,dim)
% r=rms(x)
%
% Calculate the rms r=sqrt(mean(x.*conj(x)))

if(~exist('dim','var'))
  r=sqrt(mean(x.*conj(x)));
else
  r=sqrt(mean(x.*conj(x),dim));
end
