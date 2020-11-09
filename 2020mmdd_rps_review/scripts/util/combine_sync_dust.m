function out = combine_sync_dust(dust150, correlation, ...
                                 beta_sync, beta_dust)
% out = combine_sync_dust(dust150, correlation, beta_sync, beta_dust)
%
% Test simple models that combine synchrotron and dust to explain
% the BICEP2 150 GHz BB signal.
%
% Parameters:
%   dust150      Fraction of 150 GHz r=0.2 signal that is due to dust. The 
%                remaining part of the signal will be filled in with 
%                synchrotron (the model includes *no* CMB B-modes). Can be 
%                an array.
%   correlation  Correlation coefficient for synchrotron and dust. Set to 
%                0 (default) for uncorrelated foregrounds; set to 1 for 
%                completely correlated foregrounds. 
%   beta_sync    Power law spectral index for synchrotron. Defaults to -3.
%   beta_dust    Power law spectral index for dust. Defaults to +1.6.

% 2014-05-14 CAB

% Default to uncorrelated.
if nargin < 2
  correlation = 0;
end
% Default to beta_sync = -3
if nargin < 3
  beta_sync = -3.0;
end
% Default to beta_dust = 1.6
if nargin < 4
  beta_dust = 1.6;
end

% Adjust the actual frequencies used for "150 GHz" and "100 GHz".
nu150 = 149.8;
nu100 = 96;
nu23 = 23;

% By definition, synchrotron and dust will add up to one at 150 GHz.
auto150 = ones(size(dust150));

% Calculate what fraction (in power) of the 150 GHz signal must be
% due to synchrotron (in order to explain the entire 150 GHz
% signal).
sync150 = (-correlation * sqrt(dust150) + ...
          sqrt(1 - (1 - correlation^2) * dust150)).^2;

% Extrapolate synchrotron and dust amplitudes to 100 GHz.
sync100 = sync150 * (nu100 / nu150)^(2 * beta_sync);
dust100 = dust150 * (nu100 / nu150)^(2 * beta_dust);

% Calculate 100 GHz auto spectrum.
auto100 = sync100 + dust100 + ...
          2 * correlation * sqrt(sync100 .* dust100);

% Calculate 100x150 GHz cross spectrum.
cross = sqrt(sync100 .* sync150) + sqrt(dust100 .* dust150) + ...
        correlation * sqrt(sync100 .* dust150) + ...
        correlation * sqrt(dust100 .* sync150);

% Estimate beta.
beta = log(cross ./ auto150) / log(nu100 / nu150);

% Calculate sigma for WMAP K-band x BICEP2.
% BB signal detected at 150 GHz.
r150 = 0.2;
% One sigma uncertainty, expressed in units of r, for WMAP-K x B2
% cross-correlation, extrapolating WMAP-K to 150 GHz using beta_sync.
sigmar_wmapK_b2 = 0.0071 * (nu150 / nu23)^(beta_sync + 3);
% Number of sigma by which we can rule out this model using 
% WMAP-K x B2.
sigma = (sync150 + correlation * sqrt(sync150 .* dust150)) * ...
        (r150 / sigmar_wmapK_b2);

% Calculate sigma for WMAP K-band x BICEP1.

% One sigma uncertainty, in units of r, for WMAP-K x B1
% cross-correlation, extrapolating WMAP-K to 96 GHz using beta_sync.
% The uncertainty can be reproduced by paper_plots_b2_respap14
sigmar_wmapK_b1 = 0.058 * (nu100 / nu23) ^ (beta_sync + 3);
% sync100 and dust100 are in antenna temperature, but sigmar_wmapK_b1 is in CMB temperature


% Return results.
out.sync100 = sync100;
out.dust100 = dust100;
out.sync150 = sync150;
out.dust150 = dust150;
out.auto150 = auto150;
out.auto100 = auto100;
out.cross = cross;
out.beta = beta;
out.sigma = sigma;
