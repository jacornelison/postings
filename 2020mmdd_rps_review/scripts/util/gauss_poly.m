function y = gauss_poly(p,x)
% y = gauss_poly(p,x)
%
% Compute values at x on a Gaussian shaped curve
% sitting on a polynomial baseline
%
% p(1) = peak height
% p(2) = mean
% p(3) = sigma
%
% p(4:end) = polynomial pars

y=p(1)*exp(-(x-p(2)).^2/(2*p(3)^2))+polyval2(p(4:end),x);
