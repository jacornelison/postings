function map=rotqumaps(map,rotang)
% map=rotqumap(map,rotang)
%
% rotate Q/U maps into one another by angle theta
%
% NB: rotation angle is det angle on sky in deg

if(length(rotang)==1 && length(map)>1);
  rotang=rotang*ones(size(map));
end

% degrees to radians
rotang=rotang*pi/180;

for i=1:size(map,1)
  if(rotang(i)~=0)
    for j=1:size(map,2)
      % Q/U rotate by 2*rotang
      [map(i,j).Q,map(i,j).U]=do_qu_rotation(map(i,j).Q,map(i,j).U,rotang(i));
      % polynomial and ground-subtracted parts
      if isfield(map,'Qpsub')
        [map(i,j).Qpsub,map(i,j).Upsub]=do_qu_rotation(map(i,j).Qpsub,map(i,j).Upsub,rotang(i));
        [map(i,j).Qgsub,map(i,j).Ugsub]=do_qu_rotation(map(i,j).Qgsub,map(i,j).Ugsub,rotang(i));
      end
      % deprojected parts
      if isfield(map,'Qdsub')
        [map(i,j).Qdsub,map(i,j).Udsub]=do_qu_rotation(map(i,j).Qdsub,map(i,j).Udsub,rotang(i));
      end
      % deprojection templates
      if isfield(map,'Qd')
        for k=1:length(map(i,j).Qd)
          [map(i,j).Qd{k},map(i,j).Ud{k}]=do_qu_rotation(map(i,j).Qd{k},map(i,j).Ud{k},rotang(i));
        end
      end
      % Qvar/Uvar/QUcovar have cross terms involving 4*rotang
      if isfield(map,'Qvar')
        [map(i,j).Qvar,map(i,j).Uvar,map(i,j).QUcovar]=do_var_rotation(map(i,j).Qvar,map(i,j).Uvar,map(i,j).QUcovar,rotang(i));
      end
    end
  end
end

return

%%%%%

function [Q_out,U_out]=do_qu_rotation(Q_in,U_in,rotang)

[th,r]=cart2pol(Q_in,U_in);
[Q_out,U_out]=pol2cart(th+2*rotang,r);

return

%%%%%

function [Qvar_out,Uvar_out,QUcovar_out]=do_var_rotation(Qvar_in,Uvar_in,QUcovar_in,rotang)

c=cos(4*rotang);
s=sin(4*rotang);
Qvar_out=((1+c)/2)*Qvar_in+((1-c)/2)*Uvar_in-s*QUcovar_in;
Uvar_out=((1-c)/2)*Qvar_in+((1+c)/2)*Uvar_in+s*QUcovar_in;
QUcovar_out=(s/2)*(Qvar_in-Uvar_in)+c*QUcovar_in;

return
