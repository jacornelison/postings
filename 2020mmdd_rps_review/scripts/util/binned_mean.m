% BINNED_MEAN  Average a dependent variable within bins
%   in an independent variable.
%
% M=BINNED_MEAN(X,Y) takes the mean of Y within bins in X.
% M=BINNED_MEAN(X,Y,N) uses N bins in X.
% M=BINNED_MEAN(X,Y,XC) uses bins centered at XC.
% M=BINNED_MEAN(X,Y,XC,C) applies a cut C to X before averaging.
% [M S]=BINNED_MEAN(...) finds the standard deviation as well.
% [M S B]=BINNED_MEAN(...) also returns the bin centers.
% [M S B I]=BINNED_MEAN(...) also returns the counts in each bin.

% RWO 090519

function [m s b_out i_out] = binned_mean (x, y, b, cc)

if (nargin < 3) || isempty (b)
	b = 10;
end;
if isscalar (b)
	nb = b;
	b = min (x) + (max(x) - min(x))/(nb+1) * (1:nb);
end;
b = sort (b);

if (nargin >= 4) && ~isempty (cc)
	y (~cc) = NaN;
end;

if (size(x,2) == 1)
	x = x';
end;
x = double (x);
if (size(y,2) == 1)
	y = y';
end;
y = double (y);

ib = interp1 (b, 1:length(b), x, 'nearest');
ib (x < b(1)) = 1;
ib (x > b(end)) = length(b);

tab = tab_by_bin (y, ib, length(b));

m = tab(2,:) ./ tab(1,:);
if (nargout >= 2)
	s = sqrt (tab(3,:) ./ tab(1,:) - m.*m);
end;

if (nargout >= 3)
	b_out = b;
end;

if (nargout >= 4)
	i_out = tab(1,:);
end;
