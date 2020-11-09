function rotang=curveskyrotbeam(boresight_dec,cel_dk,r,theta)
% rogang=curveskyrotbeam(boresight_dec,cel_dk,r,theta)
%
% Compute beam map rotation angle in degrees.

[dec_bol,ra_bol]=reckon(boresight_dec,0,r,theta-cel_dk-90);
az=azimuth(dec_bol,ra_bol,boresight_dec,0);
rotang = az+90-theta;

return

