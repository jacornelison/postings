function map=jackknife_map(map,dim,fit)
% map=jackknife_map(map,dim,fit)
%
% difference maps assuming they contain the same signal
%
% dim is dimension to difference over - default is 2nd (row size), if
% dim=1 we transpose map going in and out such that diff takes place
% column wise (for frq jack)
%
% Optional arg fit means fit one map to the other and rescale to
% minimize jackknife map

if(~exist('dim','var'))
  dim=2;
end

% if requested transpose so diff over freq
if(dim==1)
  map=map';
end

% unless 2 columns of maps nothing to do
if(size(map,2)==2)

  % difference the maps
  for i=1:size(map,1)

    % if second entry is all empties, it's not really a jack
    if isempty(map(i,2).T)
      % if dmap was not created yet, use the map as it is
      if i==1
        dmap(i,1)=map(i,1);
      else
        % this also expands in the first dimension
        % but only considers field name dmap and map
        % have in common
        dmap = struct_merge(dmap,map(i,1));
      end
      continue
    end

    % we do (x-y)/2 mainly by analogy to (x+y)/2
    % for split dataset jackknifes this also leaves the noise the same
    dmap(i,1).T=(map(i,1).T-map(i,2).T)/2;
    dmap(i,1).Tvar=(map(i,1).Tvar+map(i,2).Tvar)/4;

    if(isfield(map,'Tpsub'))
      dmap(i,1).Tpsub=(map(i,1).Tpsub-map(i,2).Tpsub)/2;
    end

    if(isfield(map,'Tgsub'))
      dmap(i,1).Tgsub=(map(i,1).Tgsub-map(i,2).Tgsub)/2;
    end

    if(isfield(map,'Q'))
      dmap(i,1).Q=(map(i,1).Q-map(i,2).Q)/2;
      dmap(i,1).U=(map(i,1).U-map(i,2).U)/2;
      dmap(i,1).Qvar=(map(i,1).Qvar+map(i,2).Qvar)/4;
      dmap(i,1).Uvar=(map(i,1).Uvar+map(i,2).Uvar)/4;
      dmap(i,1).QUcovar=(map(i,1).QUcovar+map(i,2).QUcovar)/4;
    end

    if(isfield(map,'Qpsub'))
      dmap(i,1).Qpsub=(map(i,1).Qpsub-map(i,2).Qpsub)/2;
      dmap(i,1).Upsub=(map(i,1).Upsub-map(i,2).Upsub)/2;
    end

    if(isfield(map,'Qgsub'))
      dmap(i,1).Qgsub=(map(i,1).Qgsub-map(i,2).Qgsub)/2;
      dmap(i,1).Ugsub=(map(i,1).Ugsub-map(i,2).Ugsub)/2;
    end

    if(isfield(map,'Qdsub'))
      dmap(i,1).Qdsub=(map(i,1).Qdsub-map(i,2).Qdsub)/2;
      dmap(i,1).Udsub=(map(i,1).Udsub-map(i,2).Udsub)/2;
    end

    if(isfield(map,'Qd'))
      for j=1:numel(map(i,1).Qd)
        dmap(i,1).Qd{j}=(map(i,1).Qd{j}-map(i,2).Qd{j})/2;
        dmap(i,1).Ud{j}=(map(i,1).Ud{j}-map(i,2).Ud{j})/2;
      end
    end
    
  end
  
  map=dmap;
  
end

% transpose back
if(dim==1)
  map=map';
end

return
