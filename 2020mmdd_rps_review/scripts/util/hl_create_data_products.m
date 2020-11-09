function result = hl_create_data_products(r, ellbins, fields, n_offdiagonal)
%
%Inputs
%  r - results structure from reduc_final_comb
%  ellbins - ell bins to include in the likelihood
%  fields - 'T' 'TE' 'B' 'EB' 'TEB' etc. which CMB temperature/polarization fields to include
%  n_offdiagonal - if non-negative, only use this many offdiagonal elements of bandpower covariance matrix
%
%Outputs
%  result.N_l - [n_field, n_field, n_ell] noise bias
%  result.C_fl - fiducial model bandpowers + noise bias
%  result.C_l_hat - real data bandpowers + noise bias
%  result.M_f - bandpower covariance matrix

  n_ell = length(ellbins);
  n_sim = size(r.sim, 3);
  n_fields = length(fields);
  n_spectra = n_fields*(n_fields+1)/2; %power spectra types

  %Spectra ordering in pipeline: 'TT','TE','EE','BB','TB','EB','ET','BT','BE'
  switch fields
    case 'T'
      N_l = r.db(ellbins, 1);
      N_l = reshape(N_l, 1, 1, []);
      C_fl = mean(r.sim(ellbins, 1, :), 3);
      C_fl = reshape(C_fl, 1, 1, []);
      C_l_hat = r.real(ellbins, 1);
      C_l_hat = reshape(C_l_hat, 1, 1, []);
    
    case 'E'
      N_l = r.db(ellbins, 3);
      N_l = reshape(N_l, 1, 1, []);
      C_fl = mean(r.sim(ellbins, 3, :), 3);
      C_fl = reshape(C_fl, 1, 1, []);
      C_l_hat = r.real(ellbins, 3);
      C_l_hat = reshape(C_l_hat, 1, 1, []);

    case 'B'
      N_l = r.db(ellbins, 4);
      N_l = reshape(N_l, 1, 1, []);

      C_fl = mean(r.sim(ellbins, 4, :), 3);
      C_fl = reshape(C_fl, 1, 1, []);

      C_l_hat = r.real(ellbins, 4);
      C_l_hat = reshape(C_l_hat, 1, 1, []);
  
    case 'EB'
      N_l = zeros(2, 2, n_ell); 
      N_l(1, 1, :) = r.db(ellbins, 3); %EE
      N_l(1, 2, :) = r.db(ellbins, 6); %EB 
      N_l(2, 1, :) = N_l(1, 2, :);
      N_l(2, 2, :) = r.db(ellbins, 4); 

      C_fl = zeros(2, 2, n_ell); 
      C_fl(1, 1, :) = mean(r.sim(ellbins, 3, :), 3); %EE
      C_fl(1, 2, :) = mean(r.sim(ellbins, 6, :), 3); %EB
      C_fl(2, 1, :) = C_fl(1, 2, :);
      C_fl(2, 2, :) = mean(r.sim(ellbins, 4, :), 3); %BB

      C_l_hat = zeros(2, 2, n_ell); 
      C_l_hat(1, 1, :) = r.real(ellbins, 3); %EE
      C_l_hat(1, 2, :) = r.real(ellbins, 6); %EB
      C_l_hat(2, 1, :) = C_l_hat(1, 2, :);
      C_l_hat(2, 2, :) = r.real(ellbins, 4); %BB

    case 'TB'
      N_l = zeros(2, 2, n_ell); 
      N_l(1, 1, :) = r.db(ellbins, 1); %TT
      N_l(1, 2, :) = r.db(ellbins, 5); %TB 
      N_l(2, 1, :) = N_l(1, 2, :);
      N_l(2, 2, :) = r.db(ellbins, 4);
      C_fl = zeros(2, 2, n_ell); 
      C_fl(1, 1, :) = mean(r.sim(ellbins, 1, :), 3); %TT
      C_fl(1, 2, :) = mean(r.sim(ellbins, 5, :), 3); %TB
      C_fl(2, 1, :) = C_fl(1, 2, :);
      C_fl(2, 2, :) = mean(r.sim(ellbins, 4, :), 3); %BB
      C_l_hat = zeros(2, 2, n_ell); 
      C_l_hat(1, 1, :) = r.real(ellbins, 1); %TT
      C_l_hat(1, 2, :) = r.real(ellbins, 5); %TB
      C_l_hat(2, 1, :) = C_l_hat(1, 2, :);
      C_l_hat(2, 2, :) = r.real(ellbins, 4); %BB
    case 'TE'
      N_l = zeros(2, 2, n_ell); 
      N_l(1, 1, :) = r.db(ellbins, 1); %TT
      N_l(1, 2, :) = r.db(ellbins, 2); %TE 
      N_l(2, 1, :) = N_l(1, 2, :);
      N_l(2, 2, :) = r.db(ellbins, 3); %EE
      C_fl = zeros(2, 2, n_ell); 
      C_fl(1, 1, :) = mean(r.sim(ellbins, 1, :), 3); %TT
      C_fl(1, 2, :) = mean(r.sim(ellbins, 2, :), 3); %TE
      C_fl(2, 1, :) = C_fl(1, 2, :);
      C_fl(2, 2, :) = mean(r.sim(ellbins, 3, :), 3); %EE
      C_l_hat = zeros(2, 2, n_ell); 
      C_l_hat(1, 1, :) = r.real(ellbins, 1); %TT
      C_l_hat(1, 2, :) = r.real(ellbins, 2); %TE
      C_l_hat(2, 1, :) = C_l_hat(1, 2, :);
      C_l_hat(2, 2, :) = r.real(ellbins, 3); %EE

    case 'TEB'
      N_l = zeros(3, 3, n_ell); 
      N_l(1, 1, :) = r.db(ellbins, 1); %TT
      N_l(1, 2, :) = r.db(ellbins, 2); %TE 
      N_l(2, 1, :) = N_l(1, 2, :);
      N_l(2, 2, :) = r.db(ellbins, 3); %EE
      N_l(1, 3, :) = r.db(ellbins, 5); %TB 
      N_l(3, 1, :) = N_l(1, 3, :);
      N_l(2, 3, :) = r.db(ellbins, 6); %EB 
      N_l(3, 2, :) = N_l(2, 3, :);
      N_l(3, 3, :) = r.db(ellbins, 4);
      C_fl = zeros(3, 3, n_ell); 
      C_fl(1, 1, :) = mean(r.sim(ellbins, 1, :), 3); %TT
      C_fl(1, 2, :) = mean(r.sim(ellbins, 2, :), 3); %TE
      C_fl(2, 1, :) = C_fl(1, 2, :);
      C_fl(2, 2, :) = mean(r.sim(ellbins, 3, :), 3); %EE
      C_fl(1, 3, :) = mean(r.sim(ellbins, 5, :), 3); %TB
      C_fl(3, 1, :) = C_fl(1, 3, :);
      C_fl(2, 3, :) = mean(r.sim(ellbins, 6, :), 3); %EB
      C_fl(3, 2, :) = C_fl(2, 3, :);
      C_fl(3, 3, :) = mean(r.sim(ellbins, 4, :), 3); %BB
      C_l_hat = zeros(3, 3, n_ell); 
      C_l_hat(1, 1, :) = r.real(ellbins, 1); %TT
      C_l_hat(1, 2, :) = r.real(ellbins, 2); %TE
      C_l_hat(2, 1, :) = C_l_hat(1, 2, :);
      C_l_hat(2, 2, :) = r.real(ellbins, 3); %EE
      C_l_hat(1, 3, :) = r.real(ellbins, 5); %TB
      C_l_hat(3, 1, :) = C_l_hat(1, 3, :);
      C_l_hat(2, 3, :) = r.real(ellbins, 6); %EB
      C_l_hat(3, 2, :) = C_l_hat(2, 3, :);
      C_l_hat(3, 3, :) = r.real(ellbins, 4); %BB

    otherwise
      crash_not_implemented
  end

  %Signal+noise sims have had noise bias (and E->B) subtracted so add it back
  C_fl = C_fl + N_l;

  %Real data had noise bias subtracted so add it back
  C_l_hat = C_l_hat + N_l;

  %Covariance is unchanged by adding a constant so it doesn't matter whether we add noise bias before calculating the covariance matrix
  switch fields
    case 'T'
      sims = squeeze(r.sim(ellbins, 1, :));
    case 'E'
      sims = squeeze(r.sim(ellbins, 3, :));
    case 'B'
      sims = squeeze(r.sim(ellbins, 4, :));
    case 'EB'
      sims = zeros(n_ell .* 3, n_sim); 
      %Interleave the spectra in the order required for covariance matrix
      %Covariance matrix order is EE, BB, EB at each ell
      sims(1:3:(n_ell .* 3), :) = r.sim(ellbins, 3, :); %EE 
      sims(2:3:(n_ell .* 3), :) = r.sim(ellbins, 4, :); %BB
      sims(3:3:(n_ell .* 3), :) = r.sim(ellbins, 6, :); %EB
    case 'TB'
      sims = zeros(n_ell .* 3, n_sim); 
      %Interleave the spectra in the order required for covariance matrix
      %Covariance matrix order is TT, BB, TB at each ell
      sims(1:3:(n_ell .* 3), :) = r.sim(ellbins, 1, :); %TT 
      sims(2:3:(n_ell .* 3), :) = r.sim(ellbins, 4, :); %BB
      sims(3:3:(n_ell .* 3), :) = r.sim(ellbins, 5, :); %TB
    case 'TE'
      sims = zeros(n_ell .* 3, n_sim); 
      %Interleave the spectra in the order required for covariance matrix
      %Covariance matrix order is TT, EE, TE at each ell
      sims(1:3:(n_ell .* 3), :) = r.sim(ellbins, 1, :); %TT 
      sims(2:3:(n_ell .* 3), :) = r.sim(ellbins, 3, :); %EE
      sims(3:3:(n_ell .* 3), :) = r.sim(ellbins, 2, :); %TE

    case 'TEB'
      sims = zeros(n_ell .* 6, n_sim); 
      %Interleave the spectra in the order required for covariance matrix
      %Covariance matrix order is TT, EE, BB, TE, EB, TB at each ell
      sims(1:6:(n_ell .* 6), :) = r.sim(ellbins, 1, :); %TT 
      sims(2:6:(n_ell .* 6), :) = r.sim(ellbins, 3, :); %EE
      sims(3:6:(n_ell .* 6), :) = r.sim(ellbins, 4, :); %BB
      sims(4:6:(n_ell .* 6), :) = r.sim(ellbins, 2, :); %TE
      sims(5:6:(n_ell .* 6), :) = r.sim(ellbins, 6, :); %EB
      sims(6:6:(n_ell .* 6), :) = r.sim(ellbins, 5, :); %TB
    

    otherwise
      crash_not_implemented
  end
  M_f = cov(sims');
 
  %Trim offdiagonal elements of covariance matrix  
  if (n_offdiagonal >=0)
    M_f_trimmed = zeros(size(M_f));
    for iSpectrum = 1:n_spectra
      for jSpectrum = 1:n_spectra
        %Get each submatrix for a pair of spectrum types
        submatrix = M_f(iSpectrum:n_spectra:(n_ell.*n_spectra), ...
          jSpectrum:n_spectra:(n_ell.*n_spectra));
        %This is an n_ell * n_ell matrix so trim in ell
        submatrix = triu(tril(submatrix, n_offdiagonal), -n_offdiagonal);
        %Insert matrix back
        M_f_trimmed(iSpectrum:n_spectra:(n_ell.*n_spectra), ...
          jSpectrum:n_spectra:(n_ell.*n_spectra)) = submatrix;
      end
    end
    M_f = M_f_trimmed;
%    M_f = triu(tril(M_f, n_offdiagonal), -n_offdiagonal);
  end

  %Add systematic uncertainty after trimming covariance matrix.
  %The trimming is because we cannot estimate offdiagonal elements well with limited MC
  %However, the systematic uncertainty covariance is not limited by MC
  %so don't trim systematic uncertainty -- IDB 2013-04-08
  %Add systematic uncertainty
  if isfield(r, 'abscal_uncer') | isfield(r, 'beam_uncer')
    model = make_hl_vector(r.real, fields, ellbins, n_ell);
  end

  if isfield(r, 'abscal_uncer')
    gain_uncer = make_hl_vector(r.abscal_uncer, fields, ellbins, n_ell);
    M_f = M_f + (cvec(gain_uncer)*rvec(gain_uncer)).* (cvec(model)*rvec(model));
  end
  if isfield(r, 'beam_uncer')
    beam_uncer = make_hl_vector(r.beam_uncer, fields, ellbins, n_ell);
    M_f = M_f + (cvec(beam_uncer)*rvec(beam_uncer)).* (cvec(model)*rvec(model));
  end

  result.N_l = N_l;
  result.C_fl = C_fl;
  result.C_l_hat = C_l_hat;
  result.M_f = M_f;
return

function vec = make_hl_vector(source, fields, ellbins, n_ell)

    switch fields
      case 'T'
        vec = source(ellbins, 1);
      case 'E'
        vec = source(ellbins, 3);
      case 'B'
        vec = source(ellbins, 4);
      case 'EB'
        vec = zeros(n_ell .* 3, 1); 
      %Interleave the spectra in the order required for covariance matrix
      %Covariance matrix order is EE, BB, EB at each ell
        vec(1:3:(n_ell .* 3)) = source(ellbins, 3); %EE 
        vec(2:3:(n_ell .* 3)) = source(ellbins, 4); %BB
        vec(3:3:(n_ell .* 3)) = source(ellbins, 6); %EB
      case 'TB'
        vec = zeros(n_ell .* 3, 1); 
      %Interleave the spectra in the order required for covariance matrix
      %Covariance matrix order is TT, BB, TB at each ell
        vec(1:3:(n_ell .* 3)) = source(ellbins, 1); %TT 
        vec(2:3:(n_ell .* 3)) = source(ellbins, 4); %BB
        vec(3:3:(n_ell .* 3)) = source(ellbins, 5); %TB
      case 'TE'
        vec = zeros(n_ell .* 3, 1); 
      %Interleave the spectra in the order required for covariance matrix
      %Covariance matrix order is TT, EE, TE at each ell
        vec(1:3:(n_ell .* 3)) = source(ellbins, 1); %TT 
        vec(2:3:(n_ell .* 3)) = source(ellbins, 3); %EE
        vec(3:3:(n_ell .* 3)) = source(ellbins, 2); %TE

      case 'TEB'
      %Interleave the spectra in the order required for covariance matrix
      %Covariance matrix order is TT, EE, BB, TE, EB, TB at each ell
        vec = zeros(n_ell .* 6, 1); 
        vec(1:6:(n_ell .* 6)) = source(ellbins, 1); %TT 
        vec(2:6:(n_ell .* 6)) = source(ellbins, 3); %EE
        vec(3:6:(n_ell .* 6)) = source(ellbins, 4); %BB
        vec(4:6:(n_ell .* 6)) = source(ellbins, 2); %TE
        vec(5:6:(n_ell .* 6)) = source(ellbins, 6); %EB
        vec(6:6:(n_ell .* 6)) = source(ellbins, 5); %TB


      otherwise
        crash_not_implemented
    end
return
