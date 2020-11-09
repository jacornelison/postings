function leg=calc_legendre(z,lmax,verbose)
% function leg=calc_legendre(z,maxl,verbose)
%
% Calculates Legendre polynomials at 'z', up to degree 'lmax', using the
% forward-recursion formulae in astro-ph/0012120v3, Appendix A, "How to
% measure CMB polarization power spectra without losing information" (2001)
% M. Tegmark & A. Oliveira-Costa.
%
% INPUTS
%   z          Arguments to Legendre polynomial and weight functions
%   lmax       Maximum ell value to compute
%   verbose    If true or non-zero, extra logging information is printed
%              to standard out.
%
% OUTPUTS
%   leg    Each output array member in leg is length(z) by (lmax+1) where
%          (e.g.) leg.P(:,ell+1) gives the m=0 Legendre function of degree ell
%          evaluated at each z. The structure members are:
%
%            P      The m=0 Legendre polynomials
%            F10    Weighting function for T-Pol cross-correlation
%            F12    Polarization weight function
%            F22    Polarization weight function
%            z      The z coordinates, for reference
%
%          Note that when lmax < 4, output may contain some extra values.

  if ~exist('verbose', 'var') || isempty(verbose)
    verbose=true;
  end

  if (verbose)
    disp('Calculating Legendre polynomials...');
    fprintf(1, 'using %d z values up to lmax=%d\n', length(z), lmax);
  end

  ticall = tic();

  % Make sure z is a column vector.
  if (size(z,1)==1)
    z = z';
  end

  if issparse(z)
    if (verbose)
      fprintf(1, 'Converting to dense z (%0.4f sparse fraction)...', ...
        (prod(size(z))-nnz(z))/prod(size(z)));
      ticfull = tic();
    end
    z = full(z);
    if (verbose)
      toc(ticfull)
    end
  end
  [P,F10,F12,F22] = calc_legendre_c(double(z), int32(lmax), int32(verbose));

  if verbose
    fprintf(1, 'Legendre calculations complete. ');
    toc(ticall)
  end

  leg.P   = P;   clear P;
  leg.F10 = F10; clear F10;
  leg.F12 = F12; clear F12;
  leg.F22 = F22; clear F22;
  leg.z   = z;   clear z;

end
