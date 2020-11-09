function [ra,dec]=azeq_to_radec(x,y,racen,deccen)
%
% [ra,dec]=azeq_to_radec(x,y,racen,deccen)
%

dec=y+deccen;
ra=x./cos(dec*pi/180)+racen;

return


