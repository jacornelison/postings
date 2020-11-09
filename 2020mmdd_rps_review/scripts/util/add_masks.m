function map=add_masks(m,map,smoothvar)
% map=add_masks(m,map,smoothvar)
%
% smoothvar default = true
%
% Places in the code where map apodization masks are defined have
% proliferated - try to do this only in one place!

if(~exist('smoothvar','var'))
  smoothvar=true;
end
  
% if only one smoothvar value given expand to all
if(numel(smoothvar)==1)
  smoothvar=repmat(smoothvar,size(map));
end

% smooth the requested Q/U maps
for i=1:numel(map)
  if(smoothvar(i))
    disp('smooth var map before use to make masks')
    maps(i)=smooth_varmaps(m,map(i));
  else
    maps(i)=map(i);
  end
end

for i=1:numel(map)
  % take Tw and Pw ap masks as traditional recipe
  
  if(~isfield(map(i),'Tw') || isempty(map(i).Tw))
    map(i).Tw=1./maps(i).Tvar;
  else
    disp('declining to replace existing Tw');
  end
  
  if(isfield(map,'Q'))
    if(~isfield(map(i),'Pw') || isempty(map(i).Pw))
      map(i).Pw=1./((maps(i).Qvar+maps(i).Uvar)/2);
    else
      disp('declining to replace existing Pw');
    end
  end
end

return
