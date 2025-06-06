function c = chebpoly(n, x)
%CHEBPOLY Chebyshev polynomial of the first kind.
%
%   CHEBPOLY(N) returns the coefficients of the polynomial of degree N.
%
%   CHEBPOLY(N, X) returns the polynomial of degree N evaluated in X.
%
%   N must be a non-negative scalar integer of class double.  X may be any
%   double array, real or complex.
%
%   These polynomials are orthogonal on the interval [-1,1], with respect to
%   the weight function w(x) = (1-x^2)^(-1/2).
%
%   See also CHEBROOT, CHEBEXTR.
%
%   Author:      Peter J. Acklam
%   Time-stamp:  2004-09-22 18:40:57 +0200
%   E-mail:      pjacklam@online.no
%   URL:         http://home.online.no/~pjacklam
%
% Taken from: http://home.online.no/~pjacklam/matlab/software/util/polyutil/chebpoly.m

   % Check number of input arguments.
   nargsin = nargin;
   error(nargchk(1, 2, nargsin));

   % Check array class.
   if ~isa(n, 'double')
      error('N must be double.');
   end

   % Check array size.
   if any(size(n) ~= 1)
      error('N must be scalar.');
   end

   % Check array values.
   if ~isreal(n) || (n < 0) || (n ~= round(n))
      error('N must be a real non-negative integer');
   end

   if nargsin == 1                  % calculate coefficients

      if n == 0                     % use explicit formula when N = 0.
         c = 1;
      elseif n == 1                 % use explicit formula when N = 1.
         c = [ 1 0 ];
      else                          % use recursive formula when N > 1.
         a = 1;
         b = [ 1 0 ];
         for k = 2 : n
            c = 2 * [ b 0 ] - [ 0 0 a ];
            a = b;
            b = c;
         end
      end

      return

      %
      % This non-recursive algorithm is faster, but less accurate.
      %

      if n == 0
         c = 1;
      else
         c = pow2(n - 1) * poly(cos((pi / n) * (0.5 : n-0.5)));
         c = round(c);              % remove numerical noise
         c(2 : 2 : n+1) = 0;
      end

   else                         % evaluate polynomial

      % Check array class.
      if ~isa(n, 'double')
         error('X must be double.');
      end

      if n == 0                 % use explicit formula when N = 0
         c = ones(size(x));
      elseif n == 1             % use explicit formula when N = 1
         c = x;
      else                      % use recursive formula when N > 1
         a = 1;
         b = x;
         for k = 2 : n
            c = 2 * x .* b - a;
            a = b;
            b = c;
         end
      end

      return

      %
      % This non-recursive algorithm is faster, but less accurate.
      %

      k = ( abs(x) < 1 );
      c = zeros( size(x) );
      c(k)  = cos(  n * acos(  x(k)  ) );
      c(~k) = cosh( n * acosh( x(~k) ) );

      % Integer input should give integer output.
      k = x == round(x);
      c(k) = round(c(k));

      % Output is real when X is real.
      k = imag(x) == 0;
      c(k) = real(c(k));

   end
