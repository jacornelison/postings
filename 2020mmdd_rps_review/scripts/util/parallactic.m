function p=parallactic(az,el,mjd,lat,lon)
% p=parallactic(az,el,mjd,lat,lon)
%
% calculate parallactic angle for J2000
%
% algorithm made up by Walt
% See analysis posting:
% http://bicep.rc.fas.harvard.edu/bkcmb/analysis_logbook/analysis/20130424_parallactic/parallactic.html
% (this is "Fix #2")

DELTA = 1e-6;
[ra1,dec1]=azel2radec(az,el,mjd,lat,lon);
[ra2,dec2]=azel2radec(az,el+DELTA,mjd,lat,lon);
p=azimuth(dec1,ra1,dec2,ra2);

return
