% [xp,yp]=rtheta_to_xpyp(rcen,thetacen,r,theta)
%
% Calculates the x_P and y_P coordinates of a point (r,theta) in detector
% pair-centered coordinates (x_P,y_P) relative to a detector pair centered
% at (rcen,thetacen).  All (r,theta) coordinates are in focal plane
% coordinates as in fp_data, and the (x_P,y_P) coordinates are as defined
% in the beam map coordinate system described here:
% http://bicep0.caltech.edu/~spuder/analysis_logbook/analysis/20131005_sidelobe_coordinates/sidelobe_coordinates.pdf
%
% r, theta can be vectors.
%
% rcen, thetacen can be scalars, or else vectors with the same length as
% r and theta.
%
% The inverse of this function is xpyp_to_rtheta.
%
% Example: Find location of A detector beams in (x_P,y_P)
%   [rcen,thetacen]=paircenter(p.r(ind.a),p.theta(ind.a),p.r(ind.b),p.theta(ind.b));
%   [xp,yp]=rtheta_to_xpyp(rcen,thetacen,p.r(ind.a),p.theta(ind.a));

function [xp,yp]=rtheta_to_xpyp(rcen,thcen,r,th)

% Azimuth back toward boresight
% Note that we're taking boresight as the north pole, so direction back toward
% boresight is always the direction of north, az=0
% az0=azimuth(90-rcen,thcen,0,0);
az0=0;

[rpr,az]=distance(90-rcen,thcen,90-r,th);
dpr=180/pi * 2*sind(rpr/2);
thpr=thcen-az-(180-az0);
xp=dpr.*cosd(thpr);
yp=dpr.*sind(thpr);

return

