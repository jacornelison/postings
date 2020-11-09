function [s,c,p]=egauss2_mmt2scp(fwhm_maj,fwhm_min,theta)
% [s,c,p]=egauss2_mmt2scp(maj,min,theta)
%
% There are several ways to parameterize an elliptical
% gaussian. This function converts from the traditional fwhm_maj,
% fwhm_min, axis-angle formulation to the sigma, c, p formulation
% as described in this posting:
% http://bicep.caltech.edu/~spuder/analysis_logbook/analysis/20120712_EllipseParameterization/EllipsesParameterization.html
%
% note that s will be in whatever units fwhm_maj,fwhm_min is in
% (c and p are dimensionless), theta in deg

% convert from fwhm to sigma
x=2*sqrt(2*log(2));
sig_maj=fwhm_maj./x;
sig_min=fwhm_min./x;

% convert angle to rad
theta=theta*pi/180;

% convert the parameters
s=sqrt((sig_maj.^2+sig_min.^2)/2);
c=sin(2*theta).*(sig_maj.^2-sig_min.^2)./(sig_maj.^2+sig_min.^2);
p=cos(2*theta).*(sig_maj.^2-sig_min.^2)./(sig_maj.^2+sig_min.^2);

return
