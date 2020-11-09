% NORMFIT2D Fit a 2-D Gaussian to data.
%   A = NORMFIT2D(X,Y,Z) performs a chi-squared fit of
%   a 2-D Gaussian to the measurements Z taken at
%   coordinates X and Y.  X and Y may be matrices or
%   vectors.
%
%   The Gaussian parameters are output in A:
%     A(1) = normalization
%     A(2) = mean in X
%     A(3) = mean in Y
%     A(4) = standard deviation in X
%     A(5) = standard deviation in Y
%     A(6) = X-Y correlation coefficient
%     A(7) = baseline offset
%
%   [A, Zout] = NORMFIT2D(...) also returns the fitted
%   Gaussian approximation of Z, Zout.
%
%   NORMFIT2D(X,Y) assumes that X and Y are data points
%   drawn from a 2-D normal distribution.  They are
%   binned, and then a 2-D Gaussian is fitted to the
%   2-D histogram.
%
%   NORMFIT2D(...,'window',true) constructs a window around
%   the region where Z is largest, and zeroes out all
%   data outside this window.  This is useful when a
%   large part of the total power is not in the main
%   beam.
%
%   NORMFIT2d(...,'A0',A0) uses A0 as initial guess in
%   minimizing chi-squared.  The definition is the same
%   as for output parameter list A.
%
%   NORMFIT2d(...,'F',F) specifies a mask of parameters
%   to fit.  F(i)=1 means parameter i will be fit;
%   F(i)=0 means that parameter i will be taken from A0.
%
%   UPDATE 20180619 TSG
%   NORMFIT2d(...,'OPTS',OPTS) uses OPTS as "OPTIONS"
%   input to FMINSEARCH.  OPTS should be defined
%   using OPTIMSET.  Use this to change tolerance and
%   max number of iterations during fitting.

function [A, Z_out, A0_out, fval] = normfit2d (X, Y, Z, varargin)

V = varargin;
if (nargin >= 3) & ischar(Z) || isstruct(Z)
	V = [{Z},V];
	Z = [];
end;
if (nargin < 3)
	Z = [];
end

if length(V)==1 && isstruct(V{1})
	OPT = V;
	clear V
else
	OPT.window = false;
	OPT.a0 = [];
	OPT.f = true(7,1);
        OPT.opts = [];
	OPT = parse_opts (OPT, V);
end

% Inputs might be X-Y pairs from normal distribution
if isempty(Z) && all(size(X)==size(Y))
	GRIDSIZE = 100;
	xx = (max(X) - min(X)) * (0:(GRIDSIZE+1))/GRIDSIZE + min(X);
	yy = (max(Y) - min(Y)) * (0:(GRIDSIZE+1))/GRIDSIZE + min(Y);
	zz = hist3 ([X(:), Y(:)], {xx yy});

	X = xx(2:(end-1));
	Y = yy(2:(end-1));
	Z = zz(2:(end-1),2:(end-1));
end;

if ((min(size(X)) == 1) && (min(size(Y)) == 1) && length(X)*length(Y)==length(Z(:)))
	[X, Y] = meshgrid (X, Y);
end;

% In window mode, pick out window, apply, and call self
if OPT.window
	% Bin in X and Y
	[mx sx bx] = binned_mean (X(:), Z(:), 100);
	[my sy by] = binned_mean (Y(:), Z(:), 100);

	% Take min of neighbors to get rid of glitchy spikes
	mxmin = min ([mx(1), mx(1:(end-1))], [mx(2:end), mx(end)]);
	mymin = min ([my(1), my(1:(end-1))], [my(2:end), my(end)]);
	mx = min (mx, mxmin);
	my = min (my, mymin);

	THRESH = 1/2;
	Xw1 = 0;
	Yw1 = 0;
	while (Xw1 < median(diff(bx))) || (Yw1 < mean(diff(by)))
		Xw0 = mean (bx (mx >= max(mx(:)) * THRESH));
		Xw1 = 5*std (bx (mx >= max(mx(:)) * THRESH));
                Yw0 = mean (by (my >= max(my(:)) * THRESH));
                Yw1 = 5*std (by (my >= max(my(:)) * THRESH));
		THRESH = THRESH * 0.75;
	end;
        cw = (X >= Xw0-Xw1) & (X <= Xw0+Xw1) & (Y >= Yw0-Yw1) & (Y <= Yw0+Yw1);
	% Correct for the fact that the input grid
	% coordinates may not be rectangular
	[nx inx] = max (sum (cw, 2));
	[ny iny] = max (sum (cw, 1));
	cw2 = cw;
	cw2 (cw (:, iny), cw (inx, :)) = 1;
	cw = cw2;

	xx = reshape (X(cw), ny, nx);
	yy = reshape (Y(cw), ny, nx);
	zz = reshape (Z(cw), ny, nx);
	OPT.window = false;
	[A tmpZout tmpA0out] = normfit2d(xx, yy, zz, OPT);
	if (nargout > 1)
		Z_out = zeros (size (Z));
		Z_out (cw) = tmpZout;
        end;
        if (nargout >= 3)
                A0_out = tmpA0out;
        end;
        return;
end;

% Condition inputs
X_ofs = mean (X(:));
Y_ofs = mean (Y(:));
X_scale = max(X(:)) - min(X(:));
Y_scale = max(Y(:)) - min(Y(:));
Z_scale = mean (abs (Z(:))) / 100;

% Initial guess
if isempty(OPT.a0)
	A0_tmp = initial_guess(X,Y,Z);
else
	A0_tmp = OPT.a0;
	% We fit in variance, not standard deviation
	A0_tmp(4:5) = A0_tmp(4:5).^2;
end
% Transform initial guess into conditioned coordinates
A0_tmp = (A0_tmp ...
        - [      0,   X_ofs,   Y_ofs,          0,       0,    0,       0]) ...
       ./ [Z_scale, X_scale, Y_scale, X_scale.^2, Y_scale.^2, 1, Z_scale];

% Run the fit on conditioned quantities
[A_tmp, Z_tmp, fval] = do_the_fit ( ...
	(X - X_ofs) / X_scale, ...
	(Y - Y_ofs) / Y_scale, ...
	Z / Z_scale, A0_tmp, OPT.f, OPT.opts);

% We fit in variance, convert now to standard deviation
A_tmp(4:5) = sqrt (abs (A_tmp(4:5)));
A0_tmp(4:5) = sqrt (abs (A0_tmp(4:5)));

% Transform output back into real coordinates
A_tmp = A_tmp ...
     .* [Z_scale, X_scale, Y_scale, X_scale, Y_scale,     1, Z_scale] ...
      + [      0,   X_ofs,   Y_ofs,       0,       0,     0,       0];

A0_tmp = A0_tmp ...
     .* [Z_scale, X_scale, Y_scale, X_scale, Y_scale,     1, Z_scale] ...
      + [      0,   X_ofs,   Y_ofs,       0,       0,     0,       0];
  
Z_tmp = Z_tmp * Z_scale;

% Select output variables
A = A_tmp;
if (nargout >= 2)
	Z_out = Z_tmp;
end;
if (nargout >= 3)
	A0_out = A0_tmp;
end;

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Reasonable guess at initial parameters for search
%

function [A0] = initial_guess (X, Y, Z)

% Compute means, correlation matrix
% for initial guess at fit parameters
Zp=Z;
Zp(Z<0)=0;
Zpsum = sum (sum (Zp));
Zsum  = sum (sum (Z));

A0(2) = sum (sum (X .* Zp)) / Zpsum;
A0(3) = sum (sum (Y .* Zp)) / Zpsum;
A0(4) = sum (sum (Zp .* (X - A0(2)).^2)) / Zpsum;
A0(5) = sum (sum (Zp .* (Y - A0(3)).^2)) / Zpsum;
A0(6) = sum (sum (Zp .* (X - A0(2)) .* (Y - A0(3)))) / Zpsum;
A0(6) = A0(6) / sqrt (A0(4) * A0(5));
A0(1) = 1;
A0(7) = 0;
A0(1) = Zsum / sum (sum (gauss2d (A0, X, Y)));
A0(7) = min(min(Z));

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% The real fitting routine
%

function [A, Z_out, fval] = do_the_fit (X, Y, Z, A0, F, opts)

% Select which parameters to fit and which to hold fixed
F = logical(F);
A0F = A0(F);

% Minimize chi2, assuming uniform uncertainties on Z
[Atmp fval] = fminsearch (@(A1) gauss2d_chi2_uniform (A1, X, Y, Z, A0, F), ...
    A0F, opts);
A = A0;
A(F) = Atmp(F);
Z_out = gauss2d (A, X, Y);

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Parse parameter-value arguments
%

function opt = parse_opts(opt,v)

for i=1:2:length(v)
  param = v{i};
  if length(v)>=i+1
    val = v{i+1};
  else
    val = true;
  end
  if ~ischar(param)
    error(['Parameter name is not a string.']);
  end
  if ~isfield(opt,lower(param))
    error(['Unrecognized parameter ' param]);
  end
  opt.(lower(param)) = val;
end

return

%%%%%%%%%%%%%%%%%%%%%%%
% Helper functions
%

function C = gauss2d_chi2_uniform (A, X, Y, Z, A0, F)
	Atmp = A0;
	Atmp(F) = A;
	C = mean (mean ((Z - gauss2d(Atmp,X,Y)).^2));
	return;

function Z = gauss2d (A, X, Y)
	xx = X - A(2);
	yy = Y - A(3);
	corrxy = A(6) * sqrt (abs (A(4) * A(5)));
	if (abs(corrxy) > 1)
		corrxy = sign (corrxy);
	end;
	C = [A(4), corrxy; corrxy, A(5)];
	Cinv = C ^ -1;
	Z = A(1) ...% / (2*pi) / sqrt(det(C)) ...
          * exp (-1/2 * (xx.^2 * Cinv(1,1) + yy.^2 * Cinv(2,2) + xx.*yy * 2*Cinv(1,2)))+A(7);

	return;

