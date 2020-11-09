function [par,err,gof,stat,cov] = matmin(gof_func,inpar,freepar,func,y,varargin)
% [par,err,gof,stat,cov] = matmin(gof_func,inpar,freepar,func,y,...)
%
% Use Minuit to fit an arbitrary function and return
% fitted parameters WITH THEIR ERRORS!
%
% gof_func is the goodness-of-fit rule to use - chisq, logl etc.
% inpar are the initial parameter values to start minimization from.
% freepar is a logical array indicating if the corresponding parameter
% par should be varied in the fit. If freepar=[] then all par are
% taken to be free.
% func is the name of an N parameter user supplied MATLAB
% function to fit which must have the form [y]=func(par,x).
% y,... are the data required to calc gof. For logl
% minimum would be n,x - for opt chisq min would be y,e,x -
% but y,e,x1,x2 with [y]=func(par,x1,x2) etc is needed
% for some problems.
%
% If freepar is a structure it will contain elements
% free = indicates which par are free
% lb = lower parameter bounds
% ub = upper parameter bounds
%
% par are the output parameter values and err the errors
% on those parameters.
% gof is the goodness-of-fit at the par values given (chi
% square or -2 times the log joint probability).
% stat is the status of the fit (see MINUIT manual - 3 is AOK).
%
% eg: [x,n]=hfill(randn(1,1000),100,-3,3);
%     [p,pe]=matmin('chisq',[30,0,1],[],'gauss',n,sqrt(n),x)
%     [p,pe]=matmin('logl',[30,0,1],[],'gauss',n,x)

if(isempty(freepar))
  freepar=ones(size(inpar));
end

if(isstruct(freepar))
  lb=freepar.lb;
  ub=freepar.ub;
  freepar=freepar.free;
else
  lb=zeros(size(inpar));
  ub=zeros(size(inpar));
end

if(numel(inpar)~=numel(freepar))
  error('freepar must be same size as inpar (or empty)')
end

[par,err,gof,stat,cov] = matminc(gof_func,inpar,freepar,lb,ub,func,y,varargin{:});

return
