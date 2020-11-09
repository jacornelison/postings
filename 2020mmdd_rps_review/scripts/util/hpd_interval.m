function [xint, xmax] = hpd_interval(x, P, cl)
% [xint, xmax] = hpd_interval(x, P, cl)
%
% Calculate highest posterior density interval for probability density 
% function P(x). Can also calculate the highest probability point based 
% on a quadratic fit to the pdf near its peak.
%
% NOTES:
%   1. This function assumes that x is specified as a regular grid and will 
%      give incorrect results otherwise.
%   2. Assumes that P(x) is unimodal, so the HPD interval is simply
%      connected.
%
% [Input]
%   x   Independent variable with pdf given by P.
%   P   Probability density function for x.
%   cl  Desired confidence level for the interval. 
%       Default value is 0.68 (i.e. one sigma interval).
%
% [Output]
%   xint  Two element array containing the lower and upper bounds of the 
%         HPD interval.

% Default confidence level is 68%.
if nargin < 3
  cl = [];
end
if isempty(cl)
  cl = 0.68;
end

x=cvec(x); P=cvec(P);

% Normalize pdf.
P = P / sum(P);
% Sort from highest to lowest probability density.
[Psort, idx] = sort(P, 1, 'descend');
% Find density threshold that contains desired level of total probability.
Pc = interp1(cumsum(Psort), Psort, cl);
% Find interval lower bound.
imax = find(P == max(P), 1, 'first');
if P(1) > Pc
  xint(1) = x(1);
else
  xint(1) = interp1(P(1:imax), x(1:imax), Pc);
end
% Find interval upper bound.
imax = find(P == max(P), 1, 'last');
if P(end) > Pc
  xint(2) = x(end);
else
  xint(2) = interp1(P(imax:end), x(imax:end), Pc);
end

if (idx(1) == 1) || (idx(1) == numel(Psort))
  % Highest probability density is at lower or upper boundary.
  xmax = x(idx(1));
else
  % Fit a quadratic function to the three highest probability points.
  fit = polyfit(x(idx(1:3)), P(idx(1:3)), 2);
  % Find x position with maximum probability density.
  xmax = -0.5 * fit(2) / fit(1);
end
