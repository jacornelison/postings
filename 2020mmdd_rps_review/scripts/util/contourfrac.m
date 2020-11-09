function [c,h]=contourfrac(x,y,z,fracs,varargin)
% [c,h]=contourfrac(x,y,z,fracs,varargin)
%
% Standard contour draws lines at specified levels - this variant
% draws them to encompass fractions of the total

% make a high res grid
xp=linspace(x(1),x(end),1000);
yp=linspace(y(1),y(end),1000);

% interp to these grid vals
[xg,yg]=meshgrid(xp,yp);
zp=interp2(x,y,z,xg,yg);

% sort and take cumsum
zpp=sort(zp(:));
cs=cumsum(zpp); cs=cs./cs(end);

% find levels
l=interp1(cs,zpp,1-fracs);

% plot
[c,h]=contour(xp,yp,zp,l,varargin{:});

return
