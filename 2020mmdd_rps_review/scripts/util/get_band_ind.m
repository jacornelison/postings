function ind = get_band_ind(p,ind,band)
% ind = get_band_ind(p,band)
%  
% get the index of the band, e.g.:
%  
% get_band_ind(p,ind,95)
% >> 1
% get_band_ind(p,ind,150)
% >> 2
% get_band_ind(p,ind,220)
% >> 3

freqs = unique(p.band(ind.la));
ind = find(freqs == band);  

% if p.band is 0 or not defined, set the index to 1 by default
if isempty(ind)
  ind=1;
end
return
