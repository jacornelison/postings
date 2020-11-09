function result = calc_sigma2(Lbins, L)
%result = calc_sigma2(Lbins, L)
%
%Calculate 68% confidence interval from a likelihood.
%Lbins - the values of the model parameter for the likelihood
%L - the likelihood values

  %Based on calc_sigma2 in bb2r_new.m

  %Need to integral normalize likelihood
  L = L / sum(L);

  [m, iMax] = max(L);

  %define 2sigma=68.27%
  sigma=0.6827;

  new_contour = m;
  q=find(abs(L - new_contour) < 0.002);
  is1=q(1);is2=q(end);
  c=sum(L(is1:is2)) ;

  % repeat until we found countour that encompasses 68% of likelhood curve
  while  c < sigma
    new_contour =  new_contour - 0.001;
    q=find(abs(L - new_contour) < 0.002);
    if isempty(q)
      break
    end
    is1=q(1);is2=q(end);
    c=sum(L(is1:is2)) ;
  end

  result.iMax = iMax;
  result.iPlusSigma = is2;
  result.iMinusSigma = is1;
  result.max = Lbins(iMax);
  result.PlusSigma = Lbins(is2) - Lbins(iMax);
  result.MinusSigma = Lbins(iMax) - Lbins(is1);
return

