function r=get_bpcov(r)
% r=get_bpcov(r)
%
% reduc_final deliberately doesn't gen the single spectra cov mat and
% associated diag error bars so we can re-do it at the start of
% reduc_final_comp and chi2. This is useful when want to make sim sets
% the same size

% generate single spectra cov mat using sig+noi sims and store
% in r structure
for i=1:length(r)
  for j=1:size(r(i).sim,2)
    % fetch s+n sim bandpower
    q=squeeze(r(i).sim(:,j,:))';

    % take the single spectra covariance matrix
    C=nancov(q);
    % store in r
    r(i).cov{j}=C;

    % take the correlation coefficient also
    CC=corrcoef(q);
    r(i).corr{j}=CC;
    
    % take diag err bars from these
    r(i).derr(:,j)=sqrt(diag(r(i).cov{j}));
  end
  for j=1:size(r(i).noi,2)
    % repeat for noise-only part
    r(i).noi_err(:,j)=sqrt(diag(nancov(squeeze(r(i).noi(:,j,:))')));
  end
end

return
