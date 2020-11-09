function [r,theta]=average_pairs(r,theta)
% average a pair's r and theta in a gonomonic projection centered at
% (ra,dec)=(0,0) that will recover the pair's observed average beam
% centroid 
%  [r,theta]=average_pairs(p.r(3:4),p.theta(3:4))
%

[dec,ra]=reckon(0,0,r,theta+90);
[x,y]=radec_to_gnomonic(ra,dec,0,0);
x=mean(x);
y=mean(y);

[ra,dec]=gnomonic_to_radec(x,y,0,0);
r=distance(0,0,dec,ra);
theta=azimuth(0,0,dec,ra)-90;

theta(theta<0)=theta(theta<0)+360;

return  
