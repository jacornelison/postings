function map=mask_pntsrc(m,map,src,width,pos)
% map=mask_pntsrc(m,map,src,width)
%
% Inject gaussian shaped lumps into variance map to mask point
% sources

if(~exist('width','var'))
  width=0.2; % in deg
end

if(~exist('pos','var'))
  pos=0; 
end

[xg,yg]=meshgrid(m.x_tic,m.y_tic);

% generate empty mask
pmask=zeros(m.ny,m.nx);

% for each source
for j=1:length(src.ra)
  
  % Inject gaussian
  pmask=pmask+egauss2([1,src.ra(j),src.dec(j),width,width*cos(src.dec(j)*pi/180),0],xg,yg);
end

%limit the max to be one
%otherwise, overlapping point sources can cause final
%mask to be <0.
if(pos)
  pmask(pmask>1)=1;
end  

% flip the mask so plateau at one with dips to zero
pmask=1-pmask;

%send pmask to output
map.pmask=pmask;

% for each map
for i=1:numel(map)

  % mask used in ft is reciprocal 
  % - mult by pnt src mask
  map(i).Tvar=1./(pmask.*1./map(i).Tvar);
  
  if(isfield(map,'Qvar'))
    map(i).Qvar=1./(pmask.*1./map(i).Qvar);
    map(i).Uvar=1./(pmask.*1./map(i).Uvar);
  end
end

return
