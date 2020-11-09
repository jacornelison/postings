function n = hfillb(vals,binedge,weights)
% n = hfillb(vals,binedge,weights)
%
% Histogram data into bins specified by a vector of edges.
%
% vals are the events to be binned
% binedge is a vector of bin edges
% weight of the events default to 1
% 
% If no output specified plots the resulting histogram.
%
% eg: hfillb(randn(100),[0,0.3,0.4,0.8,0.9,1]);

if(nargin<2)
  error('Must provide data and bin edges');
end

if(~exist('weights','var'))
  weights=[];
end

if(isempty(weights))
  weights=ones(size(vals));
end

% Make sure data is vector
vals   =vals(:)';
weights=weights(:)';

nbin=length(binedge)-1;
n = zeros(1,nbin);
for i=1:nbin
  ind=find(vals>binedge(i)&vals<=binedge(i+1));
  n(i)=sum(weights(ind));
end

if(nargout==0)
  hplotb(binedge,n);
end
