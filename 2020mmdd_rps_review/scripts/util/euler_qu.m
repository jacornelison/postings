function [lon,lat,q,u]=euler_qu(loni,lati,qi,ui,s,iau)
% [lon,lat,q,u]=euler_qu(lon,lat,q,u,select,iau)
% select = 1  (ra,dec) -> (l,b)
% select = 2  (l,b) -> (ra,dec)
%
% iau = true (default) - Q/U in IAU convention
%       false - Q/U in Healpix convention
%
% Inputs in degrees.

if ~exist('iau','var') || isempty (iau)
  iau=true;
end

if ~iau
  ui=-ui;
end

% Coordinates of galactic pole in celestial coords or celestial pole in galactic coords
s2=abs(s-3);
[lon_pole,lat_pole] = euler(0,90,s2,0);

% In the input coordinate system, the azimuth toward north is always 0. What is the
% azimuth toward the pole of the other coordinate system?
az=azimuth(lati,loni,lat_pole,lon_pole);
az(az>180)=az(az>180)-360;

% Polarization angle in input coordinate system
th=0.5*atan2(ui,qi)*180/pi;
th=wrappolang(th);

% Polarization angle in new coordinate system
if iau
  th2=th+az;
else
  th2=th-az;
end

% Polarized intensity
p=sqrt(qi.^2+ui.^2);

% Assign to q and u 
q=p.*cosd(2*th2);
u=p.*sind(2*th2);

% Put back in original convention
u=-u;

% Return new lat/lon just for kicks
[lon,lat]=euler(loni,lati,s,0);
  

return


%%%%%%%%%%%%%%%%%%%%
function x=wrappolang(x)

x(x>180)=x(x>180)-180;
x(x<-180)=x(x<-180)+180;
x(x<0)=x(x<0)+180;

return


