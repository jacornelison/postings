function [Q,U]=eb2qu(u,v,E,B,convention)
% [Q,U]=eb2qu(u,v,E,B,convention)
%
% Convert E,B visibilities to Q,U
%
% convention='iau' or 'healpix'

if(~exist('convention','var'))
  convention='iau';
end

switch convention
  case 'iau'
    % pol angle measured from N towards E
    chi=-atan2(v,u)+pi/2;
  case 'healpix'
    % pol angle measured from N towards W
    chi=atan2(v,u)-pi/2;
end

%changed sign of B (dB jun 3 2008) to agree with spicepol and healpix sign of TB and EB.
%equivalent change in qu2eb

Q=+E.*cos(2*chi)+B.*sin(2*chi);
U=+E.*sin(2*chi)-B.*cos(2*chi);

return
