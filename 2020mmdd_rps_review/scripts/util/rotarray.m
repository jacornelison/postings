function p=rotarray(p,dk)
% p=rotarray(p,dk)
%
% Rotate array info to given deck angle

% rotate array to given deck angle
p.theta=p.theta-dk;

% rotate pix_phys
[t,r]=cart2pol(p.pix_phys_x,p.pix_phys_y);
[p.pix_phys_x,p.pix_phys_y]=pol2cart(t-dk*pi/180,r);
p.pix_phys_x=round(p.pix_phys_x*1e6)/1e6;
p.pix_phys_y=round(p.pix_phys_y*1e6)/1e6;

% Don't rotate beam ellipse angle, because now the axis is defined w.r.t. the
% r,theta vector. Since theta gets rotated with dk, we do not need to rotate alpha. 
%if(isfield(p,'alpha'))
%  p.alpha=p.alpha-dk;
%end

% rotate theta_ref
if(isfield(p,'chi_thetaref'))
  p.chi_thetaref=p.chi_thetaref-dk;
end

if(isfield(p,'polofs_x'))
  tmp=p.polofs_x*cosd(dk)+p.polofs_y*sind(dk);
  p.polofs_y=p.polofs_y*cosd(dk)-p.polofs_x*sind(dk);
  p.polofs_x=tmp;
end

return
