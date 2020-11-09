function [ra,dec]=tanproj_to_radec(x,y,racen,deccen)
%
% [ra,dec]=tanproj_to_radec(x,y,racen,deccen)
%

a0=racen*pi/180; d0=deccen*pi/180;
y=y*pi/180;
x=x*pi/180;

a=a0+atan2(x,cos(d0)-y.*sin(d0));
dela=a-a0;
d=atan2(cos(dela).*(y.*cos(d0)+sin(d0)),cos(d0)-y.*sin(d0));

ra=a*180/pi;
dec=d*180/pi;

return

