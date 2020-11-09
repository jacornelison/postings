function [x,y]=radec_to_tanproj(ra,dec,racen,deccen)
%
% [x,y]=radec_to_tanproj(ra,dec,racen,deccen)
%

a=ra*pi/180; d=dec*pi/180;
a0=racen*pi/180; d0=deccen*pi/180;

dela=a-a0;
deld=d-d0;

x=cos(d).*sin(dela)./(sin(d)*sin(d0)+cos(d).*cos(d0).*cos(dela));
y=(sin(d).*cos(d0)-cos(d).*sin(d0).*cos(dela))./...
    (sin(d).*sin(d0)+cos(d).*cos(d0).*cos(dela));

x=x*180/pi;
y=y*180/pi;

return















