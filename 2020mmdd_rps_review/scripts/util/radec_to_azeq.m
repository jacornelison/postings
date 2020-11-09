function [x,y]=radec_to_azeq(ra,dec,racen,deccen)
%
% [x,y]=radec_to_azeq(ra,dec,racen,deccen)
%

y=dec-deccen;
x=(ra-racen).*cos(dec*pi/180);

return

