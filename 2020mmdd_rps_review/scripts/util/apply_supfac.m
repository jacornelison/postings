function aps=apply_supfac(aps,supfac)
% aps=apply_supfac(aps,supfac)
%
% Apply filter/beam supression factor correction

for j=1:length(aps)
  n=size(aps(j).Cs_l,3);
  aps(j).Cs_l=aps(j).Cs_l.*repmat(supfac(j).rwf,[1,1,n]);
end
  
return
