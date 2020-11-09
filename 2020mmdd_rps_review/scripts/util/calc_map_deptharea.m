function [s,a]=calc_map_deptharea(m,map,wmask,varmap)
% [s,a]=calc_map_deptharea(m,map,wmask,varmap)
%
% Calculate the depth of the deepest part of the map along with
% the effective area at this depth.  The varmap should be the 
% real variance map -- no smoothing, etc.  You can smooth the
% wmask if you want.
%
% See posting
% http://bicep0.caltech.edu/~spuder/analysis_logbook/analysis/20140309_map_depth/map_depth.pdf

if ~exist('varmap','var') || isempty(varmap)
  invvarmap=wmask;
else
  invvarmap=1./varmap;
end

% vectorize map and wmask
map=map(:);
wmask=wmask(:);
invvarmap=invvarmap(:);
invvarmap(~isfinite(wmask))=NaN;

% normalize to unit maximum
wmask=wmask./max(wmask);

% take area as sum of weights times pixel area
a=nansum(wmask)*m.pixsize.^2;

% assuming the weights are proportional to true inverse
% variance, we can construct a statistical estimate of the
% max depth in the deepest pixel
s=sqrt(nansum(map.^2.*invvarmap)/nansum(invvarmap)*nanmean(wmask));

% convert to uK per 1 deg^2 pixel
s=s*m.pixsize;

return
