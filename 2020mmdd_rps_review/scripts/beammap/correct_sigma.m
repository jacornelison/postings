function sigma_true = correct_sigma(sigma_meas,expt,year)
% sigma_true = correct_sigma(sigma_meas,expt,year)
%
% One-point correction for finite chopper aperture - see postings
% 20141010_beamprofile
% 20160414_keck2016_beamparam 
% 
% INPUT: MEASURED best-fit Gaussian sigma with a particular chopper given
%        by expt, year - can be array
% OUTPUT: CORRECTED best-fit Gaussian sigma - should have noticeable
%         correction for e.g. Keck 220s and BICEP3, negligible for Keck95
%
% Note that this is not the most accurate way to do this correction.
% Ideally you measure the B_l, then divide by the chopper B_l, and then
% fit sigma to that.  Here we just ask for the chopper correction at a
% single ell and make the assumption that the measured B_l at that ell is
% that of a Gaussian, which is not true!

% Determine the B_l correction based on the chopper used that year

% Distances vary by experiment
switch expt
  case {'bicep2','bicep3'}
    dist2mast = 195; % meters to MAPO mast
  case 'keck'
    dist2mast = 211; % meters to DSL mast
end

% Choose the chopper diameter
switch year
  case {2016,2017}
    diam = 24 * 0.0254; % Carbon fiber chopper diameter in meters
  otherwise
    diam = 18 * 0.0254; % Uberchopper diameter in meters
end

% Chopper diameter in radians
diam = atan( diam / dist2mast );

% Measured sigma in radians
sigma_meas = sigma_meas * pi/180;

% Analytical chopper correction
% Choose l=300 as a reasonable pivot point for all frequencies
l = 300;
jincarg = l*diam/2;
bl_chopper = 2*besselj(1,jincarg)./jincarg;

% Measured Gaussian B_ls - Dodelson 11.45
b_meas = exp(-0.5 * (sigma_meas).^2 * l.^2);

% Now find the "true" B_l at this ell after correcting for the chopper
b_true = b_meas/bl_chopper;

% This corresponds to a sigma of...
sig_true = sqrt(-2 * log(b_true)/l^2);
sigma_true = sig_true * 180/pi;

return