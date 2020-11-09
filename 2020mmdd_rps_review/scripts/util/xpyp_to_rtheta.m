% [r,theta]=xpyp_to_rtheta(rcen,thetacen,xp,yp)
%
% Calculates the coordinates (r,theta) of a point in focal-plane coordinates
% given the x_P and y_P coordinates relative to a detector pair centered at
% (rcen, thetacen).  All (r,theta) coordinates are in focal plane
% coordinates as in fp_data, and the (x_P,y_P) coordinates are as defined
% in the beam map coordinate system described here:
% http://bicep0.caltech.edu/~spuder/analysis_logbook/analysis/20131005_sidelobe_coordinates/sidelobe_coordinates.pdf
%
% xp, yp can be vectors.
%
% rcen, thetacen can be scalars, or else vectors with the same length as
% xp and yp.
%
% The inverse of this function is rtheta_to_xpyp.
%
% Example: Add pointing offset to A detector beams
%   [p.r(ind.a),p.theta(ind.a)]=xpyp_to_rtheta(rcen,thetacen,0.1,0);

function [r,th]=xpyp_to_rtheta(rcen,thcen,xp,yp);

% Azimuth back toward boresight
% Note that we're taking boresight as the north pole, so direction back toward
% boresight is always the direction of north, az=0
% az0=azimuth(90-rcen,thcen,0,0);
az0=0;

dpr=sqrt(xp.^2+yp.^2);
rpr=2*asind(pi/180 * dpr/2);
thpr=180/pi * atan2(yp,xp);
az=thcen-thpr-(180-az0);
[r,th]=reckon(90-rcen,thcen,rpr,az);
r=90-r;

return

