function p=displace_beamcen(p,dx,dy,ind,sym)
% p=displace_beamcen(p,dx,dy,ind,symmetric)
%
% Adds an offset dx and dy (specified in degrees) to the beam centers
% specified by p in the indices ind.
%
% If the symmetric flag = 1, ind is the regular ind structure containing ind.a and
% ind.b fields and a symmetric offset +/- dx/2 and +/- dy/2 are added.
%
% ex. pd=displace_beamcen(p,1/2,0,ind.a);
%     pd=displace_beamcen(pd,-1/2,0,ind.b);
%
%     yields the same as
%    
%     pd=displace_beamcen(p,1,0,ind,1);

if ~exist('sym','var')
  sym=0;
end

% dither beam centers 
[decd,rad]=reckon(0,0,p.r,p.theta+90);
  
% Use -dec since x and y are defined for az/el
[xd,yd]=radec_to_gnomonic(rad,-decd,0,0);

if sym
  % Add symmetric offsets:
  xd(ind.a)=xd(ind.a)+dx/2;
  yd(ind.a)=yd(ind.a)+dy/2;
  xd(ind.b)=xd(ind.b)-dx/2;
  yd(ind.b)=yd(ind.b)-dy/2;
else
  xd(ind)=xd(ind)+dx;
  yd(ind)=yd(ind)+dy;
end

% Convert offsets to ra/dec
[rad,decd]=gnomonic_to_radec(xd,-yd,0,0);

% Invert the action of reckon:
rdith=distance(0,0,decd,rad);
thetadith=azimuth(0,0,decd,rad)-90;

dtheta=thetadith-p.theta;
dtheta(dtheta<=-180)=dtheta(dtheta<=-180)+360;
p.theta=p.theta+dtheta;

p.r=rdith;

return
