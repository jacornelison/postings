function [az,el]=radec2azel(ra,dec,mjd,lat,lon)
% [az,el]=radec2azel(ra,dec,mjd,lat,lon)
%
% Convert ra,dec to az,el for given time and location on the Earth's
% surface
%
% All in/out in degrees apart from mjd in days
%
% Inputs can be any scalar or array in any combination 
%
% e.g.
% The path across the sky of the BK field center as seen from Pole
% and from Chile on Jan 1 2016:
% mjd=date2mjd(2016,1,1);
% [az,el]=radec2azel(0,-57.5,mjd+linspace(0,1,1000),-90,0);
% plot(az,el,'.-b'); grid
% [az,el]=radec2azel(0,-57.5,mjd+linspace(0,1,1000),-24,0);
% hold on
% plot(az,el,'.-r'); grid

% expand any arg which needs it
j{1}=size(ra); j{2}=size(dec); j{3}=size(mjd); j{4}=size(lat); j{5}=size(lon); 
[dum,i]=max([numel(ra),numel(dec),numel(mjd),numel(lat),numel(lon)]);

% expand all to same size/shape
if(isscalar(ra))
  ra=ra*ones(j{i});
end
if(isscalar(dec))
  dec=dec*ones(j{i});
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
if(any(size(ra)~=size(dec))|any(size(ra)~=size(mjd))|any(size(ra)~=size(lat))|any(size(ra)~=size(lon)))
  error('problems with size of inputs')
end

% convert from mjd to jd
jd=mjd+2400000.5;

% call the mex function
[az,el]=radec2azelc(ra,dec,jd,lat,lon);

% ensure output same shape as input
az=reshape(az,j{i});
el=reshape(el,j{i});

return
