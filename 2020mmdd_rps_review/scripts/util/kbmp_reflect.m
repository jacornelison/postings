function [rpos, rpnt, rort] = kbmp_reflect(pos, pnt, ort, mpos, mnorm)
% [rpos, rpnt, rort] = kbmp_reflect(pos, pnt, ort, mpos, mnorm)
% 
% Calculates reflected aperture position, boresight pointing, and boresight 
% orientation for specified mirror position.
%
% Reflected boresight pointing and boresight orientation are obtained by 
% reflecting those pointing vectors across the mirror normal axis. The 
% reflected aperture position is calculated by reflecting the actual 
% aperture position across the plane of the mirror (i.e. it is like the 
% reflected image that you see in a mirror).
%
% This function is usually not called directly. Use the wrapper function 
% keck_beam_map_pointing.m
%
% [Arguments]
%   All arguments should be arrays with shape [N,3], where the first
%   dimension is the number of samples and the second dimension are Cartesian
%   horizon coordinates.
%
%   First three arguments should be obtained from keck_mount.
%   pos     Telescope aperture position vector.
%   pnt     Telescope boresight pointing vector.
%   ort     Telescope boresight orientation vector.
%
%   Last two arguments should be obtained from keck_mirror.m
%   mpos    Mirror position vector.
%   mnorm   Mirror normal vector.
%
% [Returns]
%   All returned values have the same size/shape as the inputs, [N,3].
%
%   rpos    Reflected aperture position.
%   rpnt    Reflected boresight pointing vector.
%   rort    Reflected boresight orientation vector.

% Last update: 2014-02-08 CAB

% Check that input arguments all have matching size.
s = size(pos);
nsamp = s(1);
if (s(2) ~= 3)
  disp('[kbmp_reflect] ERROR: argument pos must have size [N,3]');
  return
end
s = size(pnt);
if (s(1) ~= nsamp) | (s(2) ~= 3)
  disp('[kbmp_reflect] ERROR: argument pnt must have size [N,3]');
  return
end
s = size(ort);
if (s(1) ~= nsamp) | (s(2) ~= 3)
  disp('[kbmp_reflect] ERROR: argument ort must have size [N,3]');
  return
end
s = size(mpos);
if (s(1) ~= nsamp) | (s(2) ~= 3)
  disp('[kbmp_reflect] ERROR: argument mpos must have size [N,3]');
  return
end
s = size(mnorm);
if (s(1) ~= nsamp) | (s(2) ~= 3)
  disp('[kbmp_reflect] ERROR: argument mnorm must have size [N,3]');
  return
end

% Calculate reflected aperture position.
%   rpos = pos + 2 * ((mpos - pos) DOT mnorm) * mnorm
dotnorm = sum((mpos - pos) .* mnorm, 2);
for i=1:3
  sclnorm(:,i) = 2 * dotnorm .* mnorm(:,i);
end
rpos = pos + sclnorm;

% Calculate reflected pointing vector.
%   rpnt = pnt - 2 * (pnt DOT mnorm) * mnorm
pdotn = sum(pnt .* mnorm, 2);
for i=1:3
  pproj(:,i) = pdotn .* mnorm(:,i);
end
rpnt = pnt - 2 * pproj;

% Calculate reflected orientation vector.
%  rort = ort - 2 * (ort DOT mnorm) * mnorm
odotn = sum(ort .* mnorm, 2);
for i=1:3
  oproj(:,i) = odotn .* mnorm(:,i);
end
rort = ort - 2 * oproj;
