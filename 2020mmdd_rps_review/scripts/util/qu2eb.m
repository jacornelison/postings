function [E,B,Evar,Bvar]=qu2eb(u,v,Q,U,Qvar,Uvar,convention)
% [E,B,Evar,Bvar]=qu2eb(u,v,Q,U,Qvar,Uvar,convention)
%
% Convert Q,U visibilities to E,B
% If variances are input transform these also
%
% convention='iau' or 'healpix'
%
% Inverse of eb2qu

if(~exist('Qvar','var'))
  Qvar=[];
end
if(~exist('Uvar','var'))
  Uvar=[];
end


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

c=cos(2*chi); s=sin(2*chi);
E=+Q.*c+U.*s;
%changed from this original one (B=-Q.*s+U.*c) to have the sign of TB and EB agree with
%spicepol and healpix from Eric.dB Jun 3 2008
B=+Q.*s-U.*c;  

if(~isempty(Qvar))
  Evar=+Qvar.*c.^2+Uvar.*s.^2;
  Bvar=+Qvar.*s.^2+Uvar.*c.^2;
end

return
