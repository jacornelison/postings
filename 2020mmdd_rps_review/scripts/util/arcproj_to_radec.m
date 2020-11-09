function [ra,dec]=arcproj_to_radec(x,y,racen,deccen)
%
% [ra,dec]=arcproj_to_radec(x,y,racen,deccen)
%

a0=racen*pi/180; d0=deccen*pi/180;
y=y*pi/180;
x=x*pi/180;

theta=sqrt(x.^2+y.^2);
d=asin(y.*cos(d0).*sin(theta)./theta + sin(d0).*cos(theta));
a=a0+asin(sin(theta).*x./(theta.*cos(d)));

% protect against identical zero
ind=theta==0; d(ind)=d0; a(ind)=a0;

ra=a*180/pi;
dec=d*180/pi;

return







