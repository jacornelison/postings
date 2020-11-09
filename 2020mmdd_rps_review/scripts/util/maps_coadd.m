function maps=maps_coadd(maps,ind,p,vmask)
% maps=maps_coadd(maps,ind,p,vmask)
%
% For each set of maps in cell array maps rebin channels given in ind
% to build co-added map using feed offset angles given in x/y.
%
% Optional arg vmask specifies to mask edge regions where variance is
% large
% +ve means where variance>n times minimum
% -ve means where int time<-n

if(~exist('vmask','var'))
  vmask=0;
end

disp('maps_coadd...');

% doing this way maintains shape of maps cell array
for i=1:numel(maps)
  maps{i}=coadd_maps_sub(maps{i},ind,p,vmask);
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function maps=coadd_maps_sub(maps,ind,p,vmask)

% do for each type of map
maps.rdoff=coadd_maps_sub_sub(maps,ind,p,vmask,'rdoff');

return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function map=coadd_maps_sub_sub(maps,ind,p,vmask,mtype)

% coadd maps

% get info about map from structure fields
x_tic=eval(sprintf('maps.%s.x_tic;',mtype));
y_tic=eval(sprintf('maps.%s.y_tic;',mtype));
elc=eval(sprintf('maps.elc;'));
dk=eval(sprintf('maps.dk;'));

% rotate x,y feed offset angles according to dk angle for this map
p=rotarray(p,dk);

% use same pixel size
edel=y_tic(2)-y_tic(1);
adel=x_tic(2)-x_tic(1);

% get existing map limits
lx=min(x_tic)-adel/2;
hx=max(x_tic)+adel/2;
ly=min(y_tic)-edel/2;
hy=max(y_tic)+edel/2;

% expand the map so that all feeds can be accomodated
lx=lx-1.5; hx=hx+1.5; % good to 60 deg
ly=ly-0.8; hy=hy+0.8;
nx=round((hx-lx)/adel);
ny=round((hy-ly)/edel);

% accumulate itime and signal
nch=size(eval(sprintf('maps.%s.map;',mtype)),3);
map.itime=zeros(ny,nx);
map.tot =zeros(ny,nx);

for j=ind
  [xg,yg]=meshgrid(x_tic,y_tic);  
  
  % add feed offsets accounting for fact that ra offset is deg on sky
  % but map is ra offset
  yg=yg+p.dec_off(j);
  xg=xg+p.ra_off_dos(j)./cos((elc-yg)*pi/180);
  
  w=eval(sprintf('maps.%s.itime(:,:,%d);',mtype,j));
  [map.x_tic,map.y_tic,itime]=hfill2(xg,yg,nx,lx,hx,ny,ly,hy,w);
  map.itime=map.itime+itime;
  
  w=eval(sprintf('maps.%s.map(:,:,%d).*maps.%s.itime(:,:,%d);',mtype,j,mtype,j));
  [map.x_tic,map.y_tic,tot]=hfill2(xg,yg,nx,lx,hx,ny,ly,hy,w);
  map.tot=map.tot+tot;
end

map.map=map.tot./map.itime;

rmfield(map,'tot');

% mask out noisy map edge regions if requested
if(vmask>0)
  % mask regions where variance is more than n times minimum value
  mapvar=1./map.itime; % prop to variance
  minv=min(mapvar(:));
  map.map(mapvar>vmask*minv)=NaN;
end
if(vmask<0)
  % mask regions where int time is less than threshold
  map.map(map.itimer<vmask)=NaN;
end

return
