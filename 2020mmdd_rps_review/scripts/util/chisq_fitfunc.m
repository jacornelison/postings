function y = chisq_fitfunc(p,x)
% y = chisq_fitfunc(p,x)
%
% Compute values at x on a scaled chisq curve
% with parameters:
%
% p(1) = height scaling parameter
% p(2) = degrees of freedom
% p(3) = x-axis scaling parameter
%
% p(4) = optional zero offset 
%
% eg: x=[0:0.1:20]; plot(x,chisq_fitfunc([10,3,2],x));
%
% This function is analogous to gauss.m - suitable for fitting to a
% histogram

if(length(p)<4)
  p(4)=0;
end
  
y=p(1)*chi2pdf((x-p(4))*p(3),p(2));
