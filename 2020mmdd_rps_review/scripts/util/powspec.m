function [ps,rb]=powspec(ad,ft1,ft2)
% [ps,rb]=powspec(ad,ft1,ft2)
%
% Make a 1D power spectrum of a 2D array
% If two input arrays makes cross power spectrum
%
% Try to make the bins the same width as ft ones and
% the centers run from just below 0 to just above maxr
%
% Twin function is realize to go the other way.

if(~exist('ft2','var'))
  ft2=[];
end

if(isempty(ft2))
  ft2=ft1;
end

m=max(ad.u_r(:));
top=m+ad.del_u(1)/2;
bot=-ad.del_u(1)/2;
span=top-bot;
nbin=round((span/ad.del_u(1)));

% Allow for taking the powspec of a stack of images
u_r=repmat(ad.u_r,[1,1,size(ft1,3)]);

% Remember abs(x).^2=x.*conj(x)=real(x).^2+imag(x).^2
% Second is fastest
[rb,ps]=hprof(u_r,ft1.*conj(ft2),nbin,bot,top);

% Kludge to make the first and last bin centers span
% 0 and maxr
rb(1)=-ad.del_u(1)/1e3;
rb(end)=m+ad.del_u(1)/1e3;

return
