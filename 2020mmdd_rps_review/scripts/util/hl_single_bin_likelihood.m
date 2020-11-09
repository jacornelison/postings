function result = hl_single_bin_likelihood(r, iBin, iSim, iSpectrum, n_likelihood_bins, use_sysuncer)
%Input Arguments:
%r	-	result structure (e.g. from reduc_final_comb)
%iBin	-	ell bin
%iSim	-	which sim to process (0 for real data)
%iSpectrum -	TT, EE, etc.
%use_sysuncer	-	add systematic calibration uncertainty?
%n_likelihood_bins	-	number of points to calculate likelihood at
%
%Output:
%  result.likelihood_bins - model bandpower at which likelihood was calculated
%  result.likelihood - the likelihood at each model

  if iSim==0
    disp('Warning: calculating likelihood for real data');
  else
    r.real(:,:) = r.simr(:,:,iSim); %Treat sim as real data

    %Exclude this sim from covariance calculation? If so, do it later so the likelihood bins are the same for all sims

  end

  %Pick points to calculate likelihood
  mean_sim = mean(r.simr(iBin, iSpectrum, :));
  std_sim = std(r.simr(iBin, iSpectrum, :));
%  scale = max(abs([mean_sim std_sim]));
  likelihood_bins = (1:n_likelihood_bins) / n_likelihood_bins * std_sim * 20 + mean_sim - 10 * std_sim; %+/- 10 sigma around simulation mean

  result.likelihood_bins = likelihood_bins;

  %Add sys uncertainty
  if use_sysuncer & ~isfield(r,'beam_uncer')
    disp('Adding beam and abscal undertainty')
    r=add_sys_uncer(r);
  end

  result.likelihood = zeros(n_likelihood_bins, 1);

  if iSpectrum == 4
    field = 'B';
  else
    crash_not_implemented
  end

  for iLike = 1:n_likelihood_bins
  %populate r.expv
    r.expv(iBin, iSpectrum) = likelihood_bins(iLike);
    result.likelihood(iLike) = hl_like_from_final(r, iBin:iBin, field, -1, 0); 
  end

  %Normalize likelihood? But normalization for negative allowed spectra (e.g. TE) will be different
  %Do normalization elsewhere or provide option for different normalize types

end
