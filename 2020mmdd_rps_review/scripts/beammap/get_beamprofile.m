function [r,profile] = get_beamprofile(ad,beam,method,norm)
% [r,profile] = get_beamprofile(ad,beam)
%
% Bins and averages a 2D beam map into a 1D radial profile in degrees
% 
% INPUTS
%   ad:      standard ad struct describing beam
%   beam:    2D beam map of dimensions ad.N_pix
%   method:  'mean' (default) - take average in radial bin
%            'median' - take median of each radial bin
%   norm:    'none' (default) - no normalization
%            'peak' - divide by max
%            'dBi' - compared to isotropic radiator
%
% OUTPUTS
%   r:       vector of radial bins, ad-dependent
%   profile: profile value of beam map in that r bin
%            If in dBi, plot 10*log10(profile)

if isempty(method)
  method = 'mean';
end
if isempty(norm)
  norm = 'none';
end

switch norm
  case 'dBi'
    % Total map power is 1
    beam = beam / nansum(beam(:));
    isomap = ones(size(beam));
    isomap = isomap / nansum(isomap(:));
    % Spread out the isotropic radiator over the whole sky
    area = ad.Field_size_deg(1)*ad.Field_size_deg(2);
    isomap = isomap * (area / (4*pi*(180/pi)^2));
    beam = beam./isomap;
end

% Radius at each point in degrees
r = ad.t_r*180/pi;
dr = ad.del_t(1)*180/pi;

% Bin out to this radius
rmax = min(abs([ad.t_val_deg{1}(1),ad.t_val_deg{1}(end),...
                ad.t_val_deg{2}(1),ad.t_val_deg{2}(end)]));

% Bin and average
%nbin = round(ad.N_pix(1)/2);
rbin = linspace(0,rmax,round(ad.N_pix(1)/2) + 1);

for ii = 1:length(rbin)
  % This kills the need for reinterpolating later
  z = beam(r >= (rbin(ii) - dr/2) & r < (rbin(ii) + dr/2));
  %z = beam(r >= rbin(ii) & r < rbin(ii+1));
  switch method
    case 'mean'
      bavg(ii) = nanmean(z(:));
    case 'median'
      bavg(ii) = nanmedian(z(:));
  end
end

% Take the average of each bin
%ravg = (rbin(1:nbin) + rbin(2:nbin+1))/2;
r = rbin;
profile = bavg;

switch norm
  case 'peak'
    profile = profile./nanmax(profile);
end

return
