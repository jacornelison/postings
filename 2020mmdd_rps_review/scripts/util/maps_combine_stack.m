function map=maps_combine_stack(maps)
% map=maps_combine_stack(maps)
%
% Combine over days

disp('maps_combine_stack...');

% stack the maps
for i=1:size(maps,2)
  for j=1:size(maps,3)
    map{i,j}=stack_maps_sub(maps(:,i,j));
  end
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function map=stack_maps_sub(maps)

% stack the maps vertically to perform sum over days

% check that src and dk angle of all maps supplied for combine are
% same
for i=1:numel(maps)
  src(i,:)=maps{i}.src;
  dk(i,1)=maps{i}.dk;
end
if(size(unique(src,'rows'),1)>1)
  error('Found non identical src when stacking maps');
end
if(size(unique(dk,'rows'),1)>1)
  error('Found non identical dk when stacking maps');
end

% copy first map pars
map.src=maps{1}.src;
map.dk=maps{1}.dk;
map.elc=maps{1}.elc;

map.aeoff =stack_maps_sub_sub(maps,'aeoff');

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function map=stack_maps_sub_sub(maps,mtype)

% rebin a stack of maps

% get 
for i=1:numel(maps)
  x_tic{i}=eval(sprintf('maps{%d}.%s.x_tic;',i,mtype));
  y_tic{i}=eval(sprintf('maps{%d}.%s.y_tic;',i,mtype));
end

% use same pixel size as first map
edel=y_tic{1}(2)-y_tic{1}(1);
adel=x_tic{1}(2)-x_tic{1}(1);

% select sufficient range to cover all maps
y_tic=horzcat(y_tic{:});
ly=min(y_tic)-edel/2;
hy=max(y_tic)+edel/2;
ny=round((hy-ly)/edel);
x_tic=horzcat(x_tic{:});
lx=min(x_tic)-adel/2;
hx=max(x_tic)+adel/2;
nx=round((hx-lx)/adel);

% accumulate hits and signal
nch=size(eval(sprintf('maps{%d}.%s.tot;',i,mtype)),3);
map.nhit=zeros(ny,nx,nch);
map.tot =zeros(ny,nx,nch);
for i=1:numel(maps)
  x=eval(sprintf('maps{%d}.%s.x_tic;',i,mtype));
  y=eval(sprintf('maps{%d}.%s.y_tic;',i,mtype));
  [x,y]=meshgrid(x,y);
  
  for j=1:nch    
    w=eval(sprintf('maps{%d}.%s.nhit(:,:,%d);',i,mtype,j));
    [map.x_tic,map.y_tic,nhit]=hfill2(x,y,nx,lx,hx,ny,ly,hy,w);
    map.nhit(:,:,j)=map.nhit(:,:,j)+nhit;
    
    w=eval(sprintf('maps{%d}.%s.tot(:,:,%d);',i,mtype,j));
    [map.x_tic,map.y_tic,tot]=hfill2(x,y,nx,lx,hx,ny,ly,hy,w);
    map.tot(:,:,j)=map.tot(:,:,j)+tot;
  end
end
  
return
