function [bincenter,mu,sig,mu_e,sig_e,n] = hprof(x,y,nbin,low,high)
% [bincen,mu,sig,mu_e,sig_e,n] = hprof(x,y,nbin,low,high)
%
% Generate a "profile histogram" of mean and sigmas
% for data y in bins of x specified by number and range.
%
% x,y are the coords of events to be profiled
% nbin defaults to 10
% low edge of the bins defaults to min(x)
% high edge of the bins defaults to max(x)
% 
% bincenter = bin center
% mu = average of points in bin
% sig = standard deviation of points in bin
% mu_e = error on mean
% sig_e = error on sigma
% n = points in bin
%
% If no output specified plots the resulting profile histogram
% using errorbar(bincenter,mu,mu_e,'*')
%
% eg: hprof(rand(1,1000),randn(1,1000),30,0,1);

if(~exist('nbin','var'))
  nbin=[];
end
if(~exist('low','var'))
  low=[];
end
if(~exist('high','var'))
  high=[];
end

% Vectorize data
x=x(:)'; y=y(:)';

% If there is non-finite data remove it
if(sum(~isfinite(y))>0)
  ind=isfinite(y);
  x=x(ind);
  y=y(ind);
 % warning('Non-finite input data (removed)');
end

if(isempty(nbin))
  nbin=10;
end
if(isempty(low))
  low=min(x(:));
end
if(isempty(high))
  high=max(x(:));
  % Need this or loose one item of data
  high=high+(high-low)*1e-9;
end

% C code for speed
[bincenter,mu,sig,n] = hprofc(double(x),double(y),double(nbin),double(low),double(high));

% calculate the errors avoiding div by zero
np=n; np(np<1)=1;
mu_e=sig./sqrt(np);
np=n-1; np(np<1)=1;
sig_e=sig./sqrt(2*np);

if(nargout==0 & nbin>1)
  errorbar(bincenter,mu,mu_e,'*');
  bw=bincenter(2)-bincenter(1);
  xlim([bincenter(1)-bw/2,bincenter(end)+bw/2]);
end
