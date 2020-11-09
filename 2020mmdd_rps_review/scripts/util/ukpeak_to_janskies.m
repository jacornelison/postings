function c=ukpeak_to_janskies(nu,beamwid)
% c=ukpeak_to_janskies(nu,beamwid)
%
% Return conversion factor required to scale point source peak heights
% taken from a CMB intensity map in units of thermodynamic uK to
% source flux in Janskies assuming a Gaussian beam of width beamwid
% (and that the source is unresolved).
%
% NB: beamwid is the "sigma" of the Gaussian beam - NOT the FWHM - if
% you have FWHM divide it by 2.35 before feeding into this function
%
% Beam width in radians!

% Calc delta intensity in W m^-2 Hz^-1  Sr^-1 K^-1
dIdT=planck_dIdT(nu,2.73);

% change to Jy Sr^-1 uK^-1
dIdT=dIdT*1e26/1e6;

% the peak height at which a source of flux 1Jy appears in a
% brightness map depends on the "concentration" of the beam -
% i.e. it's peak to integral ratio. For a 2D gaussian shape
% z=exp(-((x^2+y^2)/s^2)) the peak height is 1 and the area is
% s^2*2*pi so we simply multiply by the later
c=dIdT.*2*pi.*beamwid.^2;

return
