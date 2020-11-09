function like = hl_like_from_final(r, ellbins, fields, n_offdiagonal, ...
  loglike, products, prep)
%
%Inputs
%  r - results structure from reduc_final_comb
%  ellbins - ell bins to include in the likelihood
%  fields - 'T' 'TE' 'B' 'TEB' etc. which CMB temperature/polarization fields to include
%  n_offdiagonal - if non-negative, only use this many offdiagonal elements of bandpower covariance matrix
%  loglike - if true, return the log likelihood
%  products, prep - supply to skip recalculation. products must have been calculated from the same r structure, with only expv being different

%%Spectra ordering is same as B1

  if ~exist('products', 'var')
    products = hl_create_data_products(r, ellbins, fields, n_offdiagonal);
  end

  switch fields
    case 'T'
      %%TT only
      C_l = r.expv(ellbins, 1);
      C_l = reshape(C_l, 1, 1, []);
    case 'E'
      %%EE only
      C_l = r.expv(ellbins, 3);
      C_l = reshape(C_l, 1, 1, []);
    case 'B'
      %%BB only
      C_l = r.expv(ellbins, 4);
      C_l = reshape(C_l, 1, 1, []);

    case 'EB'
      C_l = zeros(2, 2, length(ellbins));
      C_l(1, 1, :) = r.expv(ellbins, 3); %EE
      C_l(1, 2, :) = r.expv(ellbins, 6); %EB
      C_l(2, 1, :) = C_l(1, 2, :);
      C_l(2, 2, :) = r.expv(ellbins, 4); %BB
    case 'TB'
      C_l = zeros(2, 2, length(ellbins));
      C_l(1, 1, :) = r.expv(ellbins, 1); %TT
      C_l(1, 2, :) = r.expv(ellbins, 5); %TB
      C_l(2, 1, :) = C_l(1, 2, :);
      C_l(2, 2, :) = r.expv(ellbins, 4); %BB
    case 'TE'
      C_l = zeros(2, 2, length(ellbins));
      C_l(1, 1, :) = r.expv(ellbins, 1); %TT
      C_l(1, 2, :) = r.expv(ellbins, 2); %TE
      C_l(2, 1, :) = C_l(1, 2, :);
      C_l(2, 2, :) = r.expv(ellbins, 3); %EE

    case 'TEB'
      C_l = zeros(3, 3, length(ellbins));
      C_l(1, 1, :) = r.expv(ellbins, 1); %TT
      C_l(1, 2, :) = r.expv(ellbins, 2); %TE
      C_l(2, 1, :) = C_l(1, 2, :);
      C_l(2, 2, :) = r.expv(ellbins, 3); %EE
      C_l(1, 3, :) = r.expv(ellbins, 5); %TB
      C_l(3, 1, :) = C_l(1, 3, :);
      C_l(2, 3, :) = r.expv(ellbins, 6); %EB
      C_l(3, 2, :) = C_l(2, 3, :);
      C_l(3, 3, :) = r.expv(ellbins, 4); %BB
      

    otherwise
      crash_not_implemented
  end

  %add noise bias to theory
  C_l = C_l + products.N_l;
 
  if ~exist('prep', 'var')
    prep = hamimeche_lewis_prepare(products.C_fl, products.M_f);
  end

  like = hamimeche_lewis_likelihood(prep, products.C_l_hat, C_l, loglike);
return
