function [m,map]=pad_map(m,map,p)
% [m,map]=pad_map(m,map,p)
%
% Pad the maps to increase fourier space resolution
% The boost factor here is arbitrary - the more you boost the
% higher the ell space resolution and the more correlated the
% bandpowers.

for i=1:numel(map)
  if isfield(map, 'Tvar')
    s=(p-size(map(i).Tvar))/2;
    map(i).Tvar=padarray2(map(i).Tvar,s);
  elseif isfield(map, 'kappa')
    s=(p-size(map(i).kappa))/2;
  end

  if isfield(map,'T')
    map(i).T=padarray2(map(i).T,s);
  end  
  
  if isfield(map,'Tw')
    map(i).Tw=padarray2(map(i).Tw,s);
  end
  
  if isfield(map,'Tsum')
    map(i).Tsum=padarray2(map(i).Tsum,s);
    map(i).Tdif=padarray2(map(i).Tdif,s);
    map(i).Tdifvar=padarray2(map(i).Tdifvar,s);
  end
  
  if(isfield(map,'Q'))
    map(i).Q=padarray2(map(i).Q,s);
    map(i).U=padarray2(map(i).U,s);
  end
  if(isfield(map,'Qvar'))
    map(i).Qvar=padarray2(map(i).Qvar,s);
    map(i).Uvar=padarray2(map(i).Uvar,s);
    map(i).QUcovar=padarray2(map(i).QUcovar,s);
  end
  if isfield(map,'Pw')
    map(i).Pw=padarray2(map(i).Pw,s);
  end
  
  if isfield(map,'Tgsub')
    map(i).Tgsub=padarray2(map(i).Tgsub,s);
    map(i).Tpsub=padarray2(map(i).Tpsub,s);
  end
  
  if isfield(map,'Qgsub')
    map(i).Qgsub=padarray2(map(i).Qgsub,s);
    map(i).Qpsub=padarray2(map(i).Qpsub,s);
    map(i).Ugsub=padarray2(map(i).Ugsub,s);
    map(i).Upsub=padarray2(map(i).Upsub,s);
  end

  if isfield(map,'Qdsub')
    map(i).Qdsub=padarray2(map(i).Qdsub,s);
    map(i).Udsub=padarray2(map(i).Udsub,s);
  end
    
  if(isfield(map(i),'QprojB') && ~isempty(map(i).QprojB))
    map(i).QprojB=padarray2(map(i).QprojB,s);
    map(i).UprojB=padarray2(map(i).UprojB,s);
  end
  
  if(isfield(map(i),'QprojE') && ~isempty(map(i).QprojE))
    map(i).QprojE=padarray2(map(i).QprojE,s);
    map(i).UprojE=padarray2(map(i).UprojE,s);
  end

  if(isfield(map(i),'kappa'))
    map(i).kappa=padarray2(map(i).kappa,s);
  end
  
  if isfield(map(i),'B')
    if ~isempty(map.B)
      map(i).B=padarray2(map(i).B,s);
    end
  end  
  
end

% note not all fields of m are re-calc
m.xdos=m.xdos*p/m.nx; m.ydos=m.ydos*p/m.ny;
m.nx=p; m.ny=p;

dx=m.x_tic(2)-m.x_tic(1); dy=m.y_tic(2)-m.y_tic(1);

x_tic=(1:m.nx)*dx;
x_tic=x_tic + m.x_tic(1) - x_tic(floor(s(2))+1);

y_tic=(1:m.ny)*dy;
y_tic=y_tic + m.y_tic(1) - y_tic(floor(s(1))+1);

m.x_tic=x_tic;
m.y_tic=y_tic;

% these two lines appear to be broken - they would
% put x_tic/y_tic in the last map only - and they
% were causing trouble so comment them out.
% (don't think we want x_tic/y_tic in the map struct anyway)
%map(i).y_tic=y_tic;
%map(i).x_tic=x_tic;

return

%%%%%%%%%%%%%%%%%%%%%%%%%
function x=padarray2(x,s)
% variant on default behavior of padarray which will work if pad value
% is for instance 67.5

x=padarray(x,floor(s),NaN,'pre');
x=padarray(x,ceil(s),NaN,'post');

return
