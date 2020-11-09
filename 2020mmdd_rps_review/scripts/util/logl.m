function logl = logl(pars,func,n,varargin)
% logl = logl(pars,func,n,x...)
%
% Calculate -2 times the log of the joint probability of a 
% set of observations n at values x if they are Poisson deviates from
% mean values described by function func with parameters
% pars.
%
% eg: [x,n]=hfill(randn(1,1000),100,-3,3);
%     p=fminsearch('logl',[30,0,1],[],'gauss',n,x)
%     [p,pe]=matmin('logl',[30,0,1],[],'gauss',n,x)
%
% See also chisq

% table of log factorials
persistent LFACTS

% Ensure table is big enough for this data set
% A bit wasteful to test every time and don't need to recalc
% whole table just because max value changes...
maxn=max(n);
if(length(LFACTS)<maxn+1)
  LFACTS=cumsum([0,log(1:maxn)]); 
  % this is the same as LFACTS=log(factorial(0:maxn))
end

mu=feval(func,pars,varargin{:});

% when the model prediction is identically zero it causes problems
ind=mu~=0;
n=n(ind); mu=mu(ind);
%  sum(ind)

% calculate the total likelihood
logl=-2*sum(n.*log(mu)-mu-LFACTS(n+1));
