function Z = gauss2d (A, X, Y, stdvar)
% Get values of a 2D elliptical Gaussian
% 
% Inputs:
%   A:         A(1) = normalization
%              A(2) = mean in X
%              A(3) = mean in Y
%              A(4) = variance in X 
%              A(5) = variance in Y
%              A(6) = X-Y correlation coefficient
%              A(7) = baseline offset
%   X,Y:       Axis 1/2 axes - can either be vector or grid
%              If vector, we make a meshgrid
%   stdvar:    Argument to tell if A(4) and A(5) are in standard
%              deviation or variance 
%              Default output of normfit2d is std
%              'std', 'var' (default)
% Output:
%   Z:         2D grid of dimensionality X by Y with the Gaussian
%              evaluated at each point

if ~exist('stdvar','var')
  stdvar = 'var';
end

switch stdvar
  case 'std'
    A(4) = A(4).^2;
    A(5) = A(5).^2;
end

if length(A) < 7
  A(7) = 0;
end

if ((min(size(X)) == 1) && (min(size(Y)) == 1))
  [X, Y] = meshgrid (X, Y);
end;

xx = X - A(2);
yy = Y - A(3);

corrxy = A(6) * sqrt (abs (A(4) * A(5)));
if (abs(corrxy) > 1)
  corrxy = sign (corrxy);
end;

C = [A(4), corrxy; ...
     corrxy, A(5)];
Cinv = C ^ -1;
Z = A(1) ...% / (2*pi) / sqrt(det(C)) ...
    * exp(-0.5*(xx.^2 * Cinv(1,1) + yy.^2 * Cinv(2,2) + xx.*yy * 2*Cinv(1,2)))...
    + A(7);

return