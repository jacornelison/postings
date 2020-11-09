function sa=spaceangle(az1,el1,az2,el2,units)
% sa=spaceangle(az1,el1,az2,el2,units)
%
% Calculate the space angle between vectors az1,el1 and az2,el2
%
% If optional parameter units is 'deg' inputs and outputs in
% degrees. Otherwise radians

if(~exist('units','var'))
  units='rad';
end

if(strcmp(units,'deg'))
  d2r=pi/180;
  az1=az1*d2r; az2=az2*d2r;
  el1=el1*d2r; el2=el2*d2r;
end

[x1,y1,z1]=sph2cart(az1,el1,1);
[x2,y2,z2]=sph2cart(az2,el2,1);

sa=acos(x1.*x2+y1.*y2+z1.*z2);

if(strcmp(units,'deg'))
  sa=sa/d2r;
end

return
