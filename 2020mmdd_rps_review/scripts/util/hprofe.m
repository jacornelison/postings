function [bincenter,mu,sig,mu_e,sig_e,n] = hprofe(x,y,binedge)
% [bincen,mu,sig,mu_e,sig_e,n] = hprofe(x,y,binedge)
%
% Alt version of hprof which bins between bin edges specified 
% 
% eg: hprofe(1000*rand(1,1000),randn(1,1000),[10:10:100,1000]);

% Vectorize data
x=x(:)'; y=y(:)';

if(~exist('binedge','var'))
  error('Must provide bin edges')
end

if(~exist('remindptr','var'))
  remindptr=0;
end

% NaNs will be ignored below
if(sum(~isfinite(y))>0)
  warning('Non-finite input data (removed)');
end

nbin=length(binedge)-1;

% find indecies falling in each bin
for i=1:nbin
  hprofe_ind{i}=find(x>binedge(i)&x<binedge(i+1));
end

% loop over bins
for i=1:nbin
  mu(i)=nanmean(y(hprofe_ind{i}));
  sig(i)=nanstd(y(hprofe_ind{i}));
  n(i)=length(hprofe_ind{i});
  bincenter(i)=mean(binedge(i:(i+1)));
end

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

return
