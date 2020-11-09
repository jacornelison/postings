function [wm,wv,neff]=wmean(x,w,dim)
% [wm,wv,neff]=wmean(x,w,dim)
%
% Return the weighted mean, weigted variance and
% effective number of measurements

if(~exist('dim'))
  dim=1;
end

if(isscalar(w))
  w=ones(size(x));
end

% Follow math at https://en.wikipedia.org/wiki/Weighted_arithmetic_mean under
% "Reliability weights"

% normalize the weights to sum of one
% (i.e. set V1=1 in wikipedia math)
n=ones(1,ndims(w));
n(dim)=size(w,dim);
w=w./repmat(nansum(w,dim),n);

% weighted mean
wm=nansum(w.*x,dim);

% take the sum of weights squared 
v2=nansum(w.^2,dim);

% the weighted variance with bias correction - see link above
wv=(1./(1-v2)).*nansum(w.*(x-repmat(wm,n)).^2,dim);

% the effective number of obs for the given weights -
% gets less surprisingly slowly for non-uniform weights
neff=1./v2;

return
