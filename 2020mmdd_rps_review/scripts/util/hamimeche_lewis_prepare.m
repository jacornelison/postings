function prep = hamimeche_lewis_prepare(C_fl, M_f)
%Perform the precomputation steps of the Hamimeche-Lewis likelihood approximation.
%Inputs
%  C_fl - [n_fields, n_fields, n_ell_bands] array of fiducial bandpowers
%  M_f - covariance matrix of fiducial bandpowers. Ordering of fields must match C_fl
%
%Outputs the inverse covariance matrix (M_f_inv) and square roots of C_fl (C_fl_12)
  prep.M_f_inv = inv(M_f);
  C_fl_12 = C_fl;

  %Loop ell bins
  for iBin = 1:size(C_fl,3)
    C_fl_12(:,:, iBin) = sqrtm(C_fl(:,:, iBin));
  end

  prep.C_fl_12 = C_fl_12;
return
