% [rcen,thcen]=paircenter(r1,theta1,r2,theta2)
%
% Find the pair center in focal plane coordinates given the location of the
% A and B detector beams.  This has been done in various ways in different
% bits of code; the current function chooses the point that's halfway
% between along a great circle.  I think we should standardize on this
% definition.
%
% [rcen,thcen,xprime,yprime]=paircenter(...) also returns the x' and y'
% coordinates of the A and B beams: A is at (xprime{1},yprime{1}) and
% B is at (xprime{2},yprime{2}).
%
% r1, theta1, etc. can be vectors.
%
% Example: [rc,thc,xp,yp]=paircenter(p.r(ind.a),p.theta(ind.a),p.r(ind.b),p.theta(ind.b));

function [rcen,thcen,xprime,yprime]=paircenter(r1,th1,r2,th2)

% First, find the pair center.
% Recipe: find great circle between the two points, and move halfway
% along it from A to B.  Use boresight (array center) as north pole.
[rng,az]=distance(90-r1,th1,90-r2,th2);
[rcen,thcen]=reckon(90-r1,th1,rng/2,az);
rcen=90-rcen;

% Second, find the per-detector beam positions in xprime, yprime coordinates

[xprime{1},yprime{1}]=rtheta_to_xpyp(rcen,thcen,r1,th1);
[xprime{2},yprime{2}]=rtheta_to_xpyp(rcen,thcen,r2,th2);

return

