function [chi2,pte,sig]=likeratio2sigmas(likeratio)
% [chi2,pte,sig]=likeratio2sigmas(likeratio)
%
% Convert a likelihood ratio to effective chi2 (Wilk's theorem)
% Then convert that to a probability to exceed
% Finally convert that to Gaussian "sigmas"
%
% The first step seems to be clear here:
% http://en.wikipedia.org/wiki/Likelihood-ratio_test#Distribution:_Wilks.27s_theorem
% 
% The third step seems debatable. Below is currently two sided
% defn as apparently used in Kovac et al paper - PDG may prefer one
% sided...
% (To reproduce page 11 of
% http://arxiv.org/pdf/astro-ph/0209478v1.pdf
% -norminv(8.46e-7/2)=4.92)
%
% e.g. 
% [chi2,pte,sig]=likeratio2sigmas(like_at_val/like_at_max)

chi2=-log(likeratio)*2;
pte=1-chi2cdf(chi2,1);
sig=-norminv(pte/2);

return
