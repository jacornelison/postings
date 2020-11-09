function ac = get_empty_ac(ac)
% function ac = get_empty_ac(ac)
% hand in an actual ac and it will return an 
% as structure that the same field names but
% all holding empty arrays, e.g. ac.wsum = [];
names = fieldnames(ac);
for ii =1:length(names)
  if iscell(ac.(names{ii}))
    for jj=1:length(ac.(names{ii}))
      ac.(names{ii}){jj}=[];
    end
  else
    ac.(names{ii})=[];
  end
end
return