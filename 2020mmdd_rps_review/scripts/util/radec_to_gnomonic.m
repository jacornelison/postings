function [x,y]=radec_to_gnomonic(ra,dec,racen,deccen)
%
% [x,y]=radec_to_gnomic(ra,dec,racen,deccen)
%

ra=ra*pi/180;
dec=dec*pi/180;
racen=racen*pi/180;
deccen=deccen*pi/180;

c = acos(sin(deccen).*sin(dec) + cos(deccen).*cos(dec).*cos(ra-racen));

y=(cos(deccen).*sin(dec)-sin(deccen).*cos(dec).*cos(ra-racen))./cos(c);
x=cos(dec).*sin(ra-racen)./cos(c);

x=x*180/pi;
y=y*180/pi;

return
