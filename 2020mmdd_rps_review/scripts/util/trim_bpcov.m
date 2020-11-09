function r=trim_bpcov(r,noffdiag)
%  r=trim_bpcov(r)
%  trim the bpcov to noffdiag rows off the diagonal

for i=1:size(r.cov,2)
  r.cov{i}=triu(tril(r.cov{i},noffdiag),-noffdiag);
end

return