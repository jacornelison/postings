function acc=addac(acc,ac)
% function acc=addac(acc,ac)
% This function used to live in reduc_coaddpairmaps.
% add or max according to field name
%  
% Use dynamic field names. This speeds up coaddpairmaps by a huge factor and makes life
% way easier when coadding whole seasons!

% we can only add what is present in both ac structures:
fnames=intersect(fieldnames(acc),fieldnames(ac));
% removing the non-intersecting field here from acc, messes up
% the loops in which addac is used, so just keep for now
% what is present in acc but not in ac and add what is
% overlapping
%  
%  rmfnames = setdiff(fieldnames(acc),fnames);
%  for f=1:length(rmfnames)
%    acc = rmfield(acc,rmfnames{f});
%  end

for f=1:length(fnames)
  if(strfind(fnames{f},'max'))
    x1=acc.(fnames{f});
    x2=ac.(fnames{f});
    acc.(fnames{f})=max(x1,x2);
  else
    if ~iscell(ac.(fnames{f}))
      acc.(fnames{f})=acc.(fnames{f})+ac.(fnames{f});
    else
      for k=1:numel(ac.(fnames{f}))
        acc.(fnames{f}){k}=acc.(fnames{f}){k}+ac.(fnames{f}){k};
      end
    end
  end
end

return