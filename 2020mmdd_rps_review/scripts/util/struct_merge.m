function sf = struct_merge(s1,s2)
%  function sf = struct_merge(s1,s2)
%  merges two structs into one for all
%  fieldnames they have in common
%  mapM=struct_merge(map1,map2);
names = fieldnames(s1);
for ii =1:length(names)
  name = names{ii};
  if (~isfield(s2,name) & ~strcmp(name,'B'))
    s1 = rmfield(s1,name);
  end
end
names = fieldnames(s2);
for ii =1:length(names)
  name = names{ii};
  if (~isfield(s1,name) & ~strcmp(name,'B'))
    s2 = rmfield(s2,name);
  end
end

sf = [s1;s2];

return
