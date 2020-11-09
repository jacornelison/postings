% POLYFITND  Fit a polynomial to data with multiple
%   independent variables. This is a generalized
%   equivalent to POLYFIT.
%
%   P = POLYFITND(X,Y,N) fits an Nth-order polynomial
%   to the dependent variable Y given the independent
%   variables X. Y is a simple vector. X may be given
%   as a cell array, or as a matrix. The fit is per-
%   formed using a matrix inversion (as with POLYFIT)
%   rather than a chi-squared minimization.
%
%   The terms are listed in order of increasing rank.
%   For example, in a case with two independent vari-
%   ables x1 and x2, and a second-order polynomial,
%   the polynomial is:
%
%      y = p(1)      + p(2)*x1    + p(3)*x2
%        + p(4)*x1^2 + p(5)*x1*x2 + p(6)*x2^2
%
%   [P,T] = POLYFITND(...) also returns a matrix with
%   exponent of each variable for each term.  In the
%   example above, T has values corresponding to the
%   powers of x1 and x2 in each term:
%
%      T = [0  0
%           1  0
%           0  1
%           2  0
%           1  1
%           0  2]
%
%   [P,T,YOUT] = POLYFITND(...) also returns the values
%   of the polynomial evaluated at the inputs X1,X2,...,
%   that is, the model values corresponding to Y.
%
%   If Y is a matrix of values on a grid, X1, X2, ...
%   may be specified as vectors of grid values; they
%   may also be omitted, in which case integer grid
%   indices will be filled in appropriately.
%
%   Example: [P,T] = polyfitnd({X1,X2},Y,2)
%   will fit a 2-d parabola in (X1,X2) to Y.
%
% RWO 110822

function [p,t,yout,tstr] = polyfitnd (x, y, n, varargin)

% We are expecting multiple independent variables
% and one dependent variable.  The independent
% variables are packed into the input "x" in
% one way or another.

% Number of independent variables to be determined.
k = NaN;
% Number of data points.

% Maybe we're not given any x and just have to
% assume y is on a grid.
if isempty(x)
  tmp = size(y);
  for i=1:length(tmp)
    x{i} = 1:size(y,i);
  end
end

% Maybe we have a cell array.
if iscell(x)
  x = {x{:}};
  k = length(x);
  for i=1:k
    % Maybe this is a gridded map and we only
    % have the axis ticks
    if (length(x{i})<length(y(:))) && (length(x{i})==size(y,i))
      tmp = size(y);
      tmp(i) = 1;
      x{i} = repmat(x{i}, tmp);
    end

    % Make sure we end up with a column vector
    x{i} = x{i}(:);
  end
else
  [s1, s2] = size(x);
  % each variable in its own column
  if (s1 == length(y(:)))
    k = s2;
    tmp = {};
    for i=1:k
      tmp{i} = x(:,i);
    end
    x = tmp;
  % each variable in its own row
  elseif (s2 == length(y(:)))
    k = s1;
    tmp = {};
    for i=1:k
      tmp{i} = x(i,:);
    end
    x = tmp;
  else
    error('Dependent and independent variable lengths do not match.');
  end
end

% Terms in the polynomial are given by combinations of the
% k independent variables and unity, with some manipulation
trm = combnk(1:(k+n), n);
for i=1:n
  trm(:,i) = trm(end:-1:1,i) - i;
end

% Massage y
szin = size(y);
y = y(:);

% Throw out NaNs?
cc = isfinite(y);
y = y(cc);
for i=1:k
  x{i} = x{i}(cc);
end
m = length(y);

% Construct Vandermonde matrix
V = ones(m,size(trm,1));
t = zeros(size(trm,1),k);
for i=1:size(trm,1)
  for j=1:n
    if trm(i,j)>0
      V(:,i) = V(:,i) .* x{trm(i,j)};
      t(i,trm(i,j)) = t(i,trm(i,j)) + 1;
    end 
  end
end

% Now solve least-squares problem
[Q,R] = qr(V,0);
p = R\(Q'*y);

% Construct model values of y, if desired
if (nargout>=3)
  yout = NaN * ones(szin);
  yout(cc) = V*p;
end

% Construct terms as strings, if desired
if (nargout>=4)
  tstr = {};
  for i=1:size(t,1)
    tstr{i} = '';
    for j=1:size(t,2)
      if t(i,j)>0
        if ~isempty(tstr{i})
          tstr{i} = [tstr{i} ' * '];
        end
        tstr{i} = [tstr{i} 'x' num2str(j) '^' num2str(t(i,j))];
      end
    end
    if isempty(tstr{i})
      tstr{i} = '1';
    end
  end
end

return
