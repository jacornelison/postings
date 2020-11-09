function [bincent Cs_l_bin C_l_bin] = aps_simple(ad,ft1,ft2,bintype)
% function [bincent Cs_l_bin C_l_bin] = aps_simple(ad,ft1,ft2,bintype)
% 
% Take the world's simplest angular power spectrum
%
% INPUTS
%        ad: standard struct for F-plane maps
%        ft1/ft2: FTs of real space maps to be squared and binned
%                 Send in the same for an auto spectrum
%        bintype: standard e.g. 'bicep_norm' or a number for equal
%                 binning up to lmax
%
% OUTPUTS
%        bincent: bin centers
%        Cs_l_bin: C_l * l(l+1)/2pi (i.e. D_l)
%        C_l: no l^2 scaling

if ischar(bintype)
  [binedges n] = get_bins(bintype);
else % If just a number, it's nbins
  lmax = max([ad.u_val{1} ad.u_val{2}])*2*pi;
  dl = lmax/bintype;
  binedges = 0:dl:lmax; % length = nbins + 1
  %bincent = dl/2 : dl : (lmax - dl/2);
end

for ii = 1:length(binedges)-1
    
  inds = find(ad.u_r*2*pi > binedges(ii) & ...
              ad.u_r*2*pi <= binedges(ii+1));

  %C_l(ii) = nanvar(ft(inds));
  C_l(ii) = nanmean(real(ft1(inds).*conj(ft2(inds))));
  bincent(ii) = nanmean(ad.u_r(inds))*2*pi;
  
end

Cs_l_bin = C_l.*bincent.*(bincent+1)./(2*pi);
C_l_bin = C_l;

return