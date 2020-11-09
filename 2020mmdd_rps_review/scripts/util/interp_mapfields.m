function map1=interp_mapfields(m1,map1,m2,map2,fields)
% map1 = interp_mapfields(m1,map1,m2,map2,fields);
%
% i.e. interp map2 on to map1's grid
% map = interp_mapfields(m,map,m2,map2,{'Tvar','Qvar','Uvar','QUcovar'})
%
% Either map1 and map2 projections must be the same, or map2's projection must be
% ra/dec.
%

if(~exist('fields','var')||issmpty(fields))
  fields=fieldnames(map2);
end

% get the grids
[x1,y1]=meshgrid(m1.x_tic,m1.y_tic);
[x2,y2]=meshgrid(m2.x_tic,m2.y_tic);

if ~strcmp(m1.proj,m2.proj)
  switch m1.proj
   case 'arc'
    [x1,y1]=arcproj_to_radec(x1,y1,m1.racen,m1.deccen);
   case 'ortho'
    [x1,y1]=sinproj_to_radec(x1,y1,m1.racen,m1.deccen);
   case 'radec'
    % nothing
  end  
end

for k=1:numel(fields)
  map1.(fields{k})=interp2(x2,y2,map2.(fields{k}),x1,y1);
end

return
