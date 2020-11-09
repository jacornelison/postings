function [ra,dec]=gnomonic_to_radec(x,y,racen,deccen)
%
% [ra,dec]=gnomonic_to_radec(x,y,racen,deccen)
%

x=x*pi/180;
y=y*pi/180;
deccen=deccen*pi/180;
racen=racen*pi/180;

rho=sqrt(x.^2+y.^2);
c=atan(rho);

ind=find(c~=0);

dec=zeros(size(x));
dec(ind)=asind( cos(c(ind)).*sin(deccen) + y(ind).*sin(c(ind)).*cos(deccen)./rho(ind) );
dec(c==0)=deccen;

ra=racen + atan2(x.*sin(c),rho.*cos(deccen).*cos(c)-y.*sin(deccen).*sin(c));
ra=ra*180/pi;

return
