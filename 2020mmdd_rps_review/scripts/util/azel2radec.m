function [ra,dec]=azel2radec(az,el,mjd,lat,lon)
% [ra,dec]=azel2radec(az,el,mjd,lat,lon)
%
% Convert az,el to ra,dec for given time and location on the Earth's
% surface
%
% All in/out in degrees apart from mjd in days
%
% Inputs can be any scalar or array in any combination
%
% e.g.
% Find the ra/dec when looking north at el=45 at mjd 0 in Chile on
% Jan 1 2016
% mjd=date2mjd(2016,1,1);
% [ra,dec]=azel2radec(0,45,mjd,-24,0)

% expand any arg which needs it
j{1}=size(az); j{2}=size(el); j{3}=size(mjd); j{4}=size(lat); j{5}=size(lon); 
[dum,i]=max([numel(az),numel(el),numel(mjd),numel(lat),numel(lon)]);

% expand all to same size/shape
if(isscalar(az))
  az=az*ones(j{i});
end
if(isscalar(el))
  el=el*ones(j{i});
end
if(isscalar(mjd))
  mjd=mjd*ones(j{i});
end
if(isscalar(lat))
  lat=lat*ones(j{i});
end
if(isscalar(lon))
  lon=lon*ones(j{i});
end

% check all now same size
if(any(size(az)~=size(el))|any(size(az)~=size(mjd))|any(size(az)~=size(lat))|any(size(az)~=size(lon)))
  error('problems with size of inputs')
end

% convert from mjd to jd
jd=mjd+2400000.5;

% call the mex function
[ra,dec]=azel2radecc(az,el,jd,lat,lon);

return
