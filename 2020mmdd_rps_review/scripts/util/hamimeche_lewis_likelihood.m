function lik = hamimeche_lewis_likelihood(prep, C_l_hat, C_l, loglike)
%Computes the likelihood
%Input
%  prep - fiducial model structure returned by hamimeche_lewis_prepare()
%  C_l_hat - [n_field, n_field, n_ell_band] array of data bandpowers. Including noise bias. Ordering should match hamimeche_lewis_prepare() arguments
%  C_l - array of theory+noise bandpowers
%  loglike - if true, return the log-likelihood. Otherwise, return the likelihood

% Fractional tolerance for imaginary likelihood values.
  itol = 1e-12;

  n_fields = size(C_l,1); %Gaussian fields in map
  n_spectra = n_fields*(n_fields+1)/2; %power spectra types
  n_ell_bins = size(C_l,3);
  X_g = zeros(n_spectra*n_ell_bins, 1);

  %From http://bmode.caltech.edu/~bicep/analysis_logbook_north/20121030_toy_model_b/#section6
  %Step 8 on that posting
  %Calculate X_gl
  %Loop bins
  for iBin = 1:n_ell_bins
    C_l_inv_12 = inv(sqrtm(squeeze(C_l(:,:,iBin))));
    [u, d] = eig(C_l_inv_12 * squeeze(C_l_hat(:,:, iBin)) * C_l_inv_12);
    %Applying g to zero offdiagonal elements might give strange result
    gd = sign(diag(d) - 1) .* sqrt(2 * (diag(d) - log(diag(d)) - 1));
    gd = diag(gd);
    xg_elements = prep.C_fl_12(:,:,iBin) * u*gd*(u') * prep.C_fl_12(:,:,iBin);
    X_g(n_spectra*(iBin-1)+1:n_spectra*iBin) = vecp(squeeze(xg_elements));
  end

  lik = -0.5*X_g' *prep.M_f_inv*X_g;
  
  % Commented this out. Function will now just go ahead and return
  % complex values when they occur.
  % 2015-02-04 CAB
  % 
  %Check for crazy model generating unphysical negative power beyond noise bias
  % if ~isreal(lik)
  %   % If the imaginary part is really small, just ignore it.
  %   if abs(imag(lik)) / abs(real(lik)) < itol
  %     lik = real(lik);
  %   else
  %     disp('Imaginary # detected -> Setting likelihood to 0');
  %     lik = -1e20;
  %   end
  % end

  if ~loglike
    lik = exp(lik);
  end
return

%function vec = vecp(M) %See vecp.m
%get unique elements of symmetric matrix in a vector  
%  dim = size(M, 1);
%  vec = zeros(dim*(dim+1)/2, 1);
%  
%  counter = 1;
%
%  for iDiag = 0:dim-1
%    vec(counter:counter+dim-iDiag-1) = diag(M, iDiag);
%    counter = counter + dim - iDiag;
%  end
%return
