function [fwhm_maj,fwhm_min,theta]=egauss2_scp2mmt(s,c,p)
%  [fwhm_maj,fwhm_min,theta]=egauss2_mmt2scp(s,c,p)
%
% There are several ways to parameterize an elliptical
% gaussian. This function converts from  sigma, c, p formulation  to the traditional fwhm_maj,
% fwhm_min, axis-angle formulation 
% as described in this posting:
% http://bicep.caltech.edu/~spuder/analysis_logbook/analysis/20120712_EllipseParameterization/EllipsesParameterization.html
%
% note that fwhm_maj,fwhm_min will be in whatever units s is in
% (c and p are dimensionless), theta in deg

e = sqrt(c.^2+p.^2);
sig_maj=s.*sqrt(1+e);
sig_min=s.*sqrt(1-e);
theta = rad2deg(asin(c./(e))/2);
theta(p<0) = 90-theta(p<0);

x=2*sqrt(2*log(2));
fwhm_maj = sig_maj.*x;
fwhm_min = sig_min.*x;

%in the case of circular beams theta is undefined, set it to 0:
theta(fwhm_maj==fwhm_min)=0;

return

