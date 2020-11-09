function rho = calc_rho(aps_BB,Brp1,Bstd,bins)
% rho = calc_rho(aps_BB,Brp1,Bstd,bins)
% 
% Simple rho estimator, i.e. value of r curve that best matches an input
% spectrum (or any fiducial model really)
% 
% INPUTS
%         aps_BB: BB points that we want to scale r to
%         Brp1:   BB curve for r = 0.1 (or whatever model you want to
%                 scale)
%         Bstd:   weights for aps_BB, can be empty for uniform weighting
%         bins:   indices of aps_BB/Brp1 to use 
%                 For standard BICEP binning, use 2:6
%                 can be empty to use all
% 
% OUTPUTS
%         rho:    single number that scales Brp1 to match aps_BB

if isempty(Bstd)
  Bstd = ones(size(aps_BB));
end
if isempty(bins)
  bins = 1:length(aps_BB);
end

w = Brp1./Bstd.^2;
Brp1sum = sum(Brp1(bins).*w(bins))/sum(w(bins));

Bcontam = sum(aps_BB(bins).*w(bins))/sum(w(bins));
rho = 0.1*Bcontam/Brp1sum; % Because model is r=0.1

return