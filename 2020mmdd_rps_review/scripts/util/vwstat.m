function [xbar,xvar]=vwstat(x,v,vpow,dim)
% [xbar,xvar]=vwstat(x,v,vpow,dim)
%
% Return the mean and variance of a group of observations of value x
% and variance v combining using the vpow power of the input
% variances combining along dimension dim

% Default is standard minimum variance combine
if(~exist('vpow'))
  vpow=[];
end

if(~exist('dim'))
  dim=[];
end

if(isempty(vpow))
  vpow=-1;
end

if(isempty(dim))
  dim=1; % default to act over col like mean etc
end

% convoluted code so actually works for n-dim input
z=v.^vpow;
n=ones(1,ndims(z));
n(dim)=size(z,dim);
w=z./repmat(nansum(z,dim),n);

xbar=nansum(w.*x,dim);
xvar=nansum((w.^2).*v,dim);

return
