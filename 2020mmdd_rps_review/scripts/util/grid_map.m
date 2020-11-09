% GRID_MAP  Accumulates data from scans to make maps.
%
% ZI = GRID_MAP(X,Y,Z,XI,YI) makes a map from the data in
%   X, Y, and Z, with grid points given by XI and YI.  The
%   data points are binned according to the grid, and averaged
%   within each bin to give ZI.
%
% ZI = GRID_MAP(X,Y,Z, NX, NY) uses NX bins in the X direction
%   and NY bins in the Y direction.
%
% ZI = GRID_MAP(X,Y,Z) uses NX=NY=10.
%
% ZI = GRID_MAP(X,Y,Z,XI,YI,CC) applies the logical mask array
%   CC for data selection before making the map.
%
% [ZI SI] = GRID_MAP(...) also returns the standard deviation
%   within each grid bin.
%
% [ZI SI XI YI] = GRID_MAP(...) also returns the grid vectors
%   XI and YI.
%
% [ZI SI XI YI I] = GRID_MAP(...) also returns the bin index
%   for each input data point.

% RWO 090729

function [ZI SI XI_out YI_out I_out] = grid_map (X, Y, Z, XI, YI, CC)

% If given # bins, or not specified, construct
% XI and YI vectors
if (nargin < 4) || isempty (XI)
	XI = 10;
end;
if isscalar (XI)
	nb = XI;
	XI = min (X) + (max(X) - min(X))/(nb+1) * (1:nb);
end;
XI = sort (XI);

if (nargin < 5) || isempty (YI)
        YI = 10;
end;
if isscalar (YI)
        nb = YI;
        YI = min (Y) + (max(Y) - min(Y))/(nb+1) * (1:nb);
end;
YI = sort (YI);

if (nargin >= 6) && ~isempty (CC)
	Z(~CC) = NaN;
end;

% Check that we have row vectors, and class double
if (size(X,2) == 1)
	X = X';
end;
X = double (X);
if (size(Y,2) == 1)
	Y = Y';
end;
Y = double (Y);
if (size(Z,2) == 1)
        Z = Z';
end;
Z = double (Z);

% Convert X and Y to a single integer bin index
ii = interp1 (XI, 1:length(XI), X, 'nearest');
ii (X < XI(1)) = 1;
ii (X > XI(end)) = length(XI);
jj = interp1 (YI, 1:length(YI), Y, 'nearest');
jj (Y < YI(1)) = 1;
jj (Y > YI(end)) = length(YI);
ibin = double(jj + (ii-1) * length(YI));

% Handle NaN values of bin index -- these break the mex
cnan = ~isfinite(ibin);
ibin(cnan) = 0;
Z(cnan) = NaN;

% Do the tabulation in mex for speed
tab = tab_by_bin (Z, ibin, length(XI) * length(YI));
m = tab(2,:) ./ tab(1,:);
ZI = reshape (m, length(YI), length(XI));
if (nargout >= 2)
	s = sqrt (tab(3,:) ./ tab(1,:) - m.*m);
	SI = reshape (s, length(YI), length(XI));
end;

if (nargout >= 3)
	XI_out = XI;
end;
if (nargout >= 4)
	YI_out = YI;
end;

if (nargout >= 5)
	I_out = tab(1,:);
end;
