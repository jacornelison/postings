function map=cal_coadd_maps(map,ukpervolt)
% map=cal_coadd_maps(map,ukpervolt)
%
% Apply abs cal scale factor to coadd maps

% If ukpervolt is single valued, expand to size of map first dimension (freq.)
if(size(map,1)>1 & numel(ukpervolt)==1)
  ukpervolt=repmat(ukpervolt,[size(map,1),1]);
end

% apply the cal factors
for i=1:size(map,1)
  for j=1:size(map,2)
    if(isfield(map,'nside'))
      % healpix map
      map(i,j).map=map(i,j).map*ukpervolt(i);
    else
      % ra/dec map
      if(isfield(map,'T'))
        map(i,j).T=map(i,j).T*ukpervolt(i);
      end
      
      if(isfield(map,'Tvar'))
        map(i,j).Tvar=map(i,j).Tvar*ukpervolt(i)^2;
      end
      
      if(isfield(map,'Tpsub'))
        map(i,j).Tpsub=map(i,j).Tpsub*ukpervolt(i);
        map(i,j).Tgsub=map(i,j).Tgsub*ukpervolt(i);
      end
      
      if(isfield(map,'Q'))
        map(i,j).Q=map(i,j).Q*ukpervolt(i);
        map(i,j).U=map(i,j).U*ukpervolt(i);
      end
      
      if(isfield(map,'Qvar'))
        map(i,j).Qvar=map(i,j).Qvar*ukpervolt(i)^2;
        map(i,j).Uvar=map(i,j).Uvar*ukpervolt(i)^2;
        map(i,j).QUcovar=map(i,j).QUcovar*ukpervolt(i)^2;        
      end

      if(isfield(map,'Qpsub'))
        map(i,j).Qpsub=map(i,j).Qpsub*ukpervolt(i);
        map(i,j).Upsub=map(i,j).Upsub*ukpervolt(i);
      end

      if(isfield(map,'Qgsub'))
        map(i,j).Qgsub=map(i,j).Qgsub*ukpervolt(i);
        map(i,j).Ugsub=map(i,j).Ugsub*ukpervolt(i);
      end

      if(isfield(map,'Qd'))
        for k=1:numel(map(i,j).Qd)
          map(i,j).Qd{k}=map(i,j).Qd{k}*ukpervolt(i);
          map(i,j).Ud{k}=map(i,j).Ud{k}*ukpervolt(i);
        end
      end

      if(isfield(map,'Qdsub'))
        map(i,j).Qdsub=map(i,j).Qdsub*ukpervolt(i);
        map(i,j).Udsub=map(i,j).Udsub*ukpervolt(i);
      end
      
      if(isfield(map,'D'))
        map(i,j).D=map(i,j).D*ukpervolt(i);
        map(i,j).Dvar=map(i,j).Dvar*ukpervolt(i)^2;		
      end
      
      if(isfield(map,'Dpsub'))
        map(i,j).Dpsub=map(i,j).Dpsub*ukpervolt(i);
        map(i,j).Dgsub=map(i,j).Dgsub*ukpervolt(i);
      end

      if(isfield(map,'Dd'))
        for k=1:numel(map(i,j).Dd)
          map(i,j).Dd{k}=map(i,j).Dd{k}*ukpervolt(i);
        end
      end

      if(isfield(map,'Ddsub'))
        map(i,j).Ddsub=map(i,j).Ddsub*ukpervolt(i);
      end
      
% 20130926 GPT - Tsum/Tdif deprecated, now T and D
%       if(isfield(map,'Tsum'))
%         map(i,j).Tsum=map(i,j).Tsum*ukpervolt(i);
%         map(i,j).Tdif=map(i,j).Tdif*ukpervolt(i);
%         map(i,j).Tdifvar=map(i,j).Tdifvar*ukpervolt(i)^2;		
%       end
      
    end
  end
end
  
return
