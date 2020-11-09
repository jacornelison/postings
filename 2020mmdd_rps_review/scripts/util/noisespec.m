function ns=noisespec(p,f,ukarcmin)
% ns=noisespec(p,f)
%
% Generic noise power spectrum function
% Analog of powlaw.m, gauss.m etc
%
% p(1) = white noise level, in uK^2 (or similar unit) unless switch used
% p(2) = 1/f knee
% p(3) = exponent of power law component (-1 for true 1/f)
%
% f = frequency (or ells)
%
% ukarcmin (optional) = p(1) is in (u)K-arcmin - convert to power
%
% e.g.:
% l=2:1000; ns=noisespec([1,10,-2],l); loglog(l,ns);
% 
% l=2:4096;
% cls(:,1)=noisespec([1e-6,200,-4],l);
% cls(:,2)=0;
% cls(:,3)=noisespec([1e-6,10,-1],l);
% cls(:,4)=cls(:,3);
% write_fits_cls('cmb-s4/noise_cls.fits',l,cls);

if(~exist('ukarcmin','var'))
  ukarcmin=false;
end

if(ukarcmin)
  omega_arcmin=4*pi/(41253*3600); % solid angle of arcmin^2 pixel
  p(1)=omega_arcmin*p(1)^2;
end

ns=p(1)*(1+(f/p(2)).^p(3));

return
