function [coeff T] = ffbm_deprojmap(map,ad,sigma,modes)
% function [coeff T] = ffbm_deprojmap(map,ad,sigma,modes)
%
% Function to return the best-fit map-space Gaussian templates corresponding
% to our standard deprojection templates
%
% INPUTS
% 
%   map:     map of the A-B pair difference beam
%   ad:      normal struct describing pixelization
%   sigma:   nominal beamwidth to generate templates to fit (degrees)
%   modes:   which modes to regress
%
% OUTPUTS
% 
%   coeff:   best fit coefficients for deproj modes
%   T:       Maps of the modes, such that coeff(i).*T{i} is what you want to
%            subtract from map to get the deprojected version.
%            T{1} = relgain
%            T{2} = diff x pointing
%            T{3} = diff y pointing
%            T{4} = diff beamwidth
%            T{5} = diff ellip plus
%            T{6} = diff ellip cross
%
% USAGE
%
%   To generate a map-space deprojection of a pair in our fiducial 2014
%   composite beam maps:
%   >> load('beammaps/maps_composite/ffbm_2014_all_onlycommonab');
%   >> diffmap = map(22).T - map(21).T; % First RGL pair diff
%   >> [coeff T] = ffbm_deprojmap(diffmap,ad,0.3); % nominal beamwidth
%   >> imagesc(diffmap-coeff(1).*T{1}-coeff(2).*T{2}-coeff(3).*T{3});
%   which would give you the map-space equivalent of dp1100

if(isempty(modes))
  modes = 1:6;
end

% Generate map space templates
T = ffbm_makelintemp(sigma,ad);

% Reshape into 1d vectors
for ii = 1:numel(modes)
  A(:,ii) = T{modes(ii)}(:);
end
map = map(:);

% Regress templates
coeff = A\map;

return
