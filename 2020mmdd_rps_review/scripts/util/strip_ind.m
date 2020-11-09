function ind=strip_ind(ind,i)
% ind=strip_ind(ind,i)
%
% Strip down ind array to contain only indices in i
% Intended to strip to single rx as in:
%
% ind=strip_ind(ind,find(ismember(p.rx,[2])))

names=fieldnames(ind);

for j=1:length(names)
  d=getfield(ind,char(names(j)));
  keepind=ismember(d,i);
  ind=setfield(ind,char(names(j)),d(keepind));
end

return
