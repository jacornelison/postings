function [mg,mapg,m,map,mapopt]=map2mapg(varargin)
% [mg,mapg,m,map,mapopt]=map2mapg(mapfile,scale)
%
% [mg,mapg,m,map,mapopt]=map2mapg(m,map,mapopt,scale)
%
% load mapfile and rotate to galactic coords
% attempts to properly account for Q/U transformation
% crudely applies emperical coordinate transformation based on euler function

% allows for variable input arguments, either filename or m,map combo
% directly.  Sets default values for pmopt and ukpervolt
[m,map,mapopt,scale]=setdef(varargin,nargin);

% create coordinate grid in l/b larger than field
mg=setmg(m,scale);
[l,b]=meshgrid(mg.x_tic,mg.y_tic);
% ra/dec points in l/b coord sys.
[ra,dec]=euler(l,b,2,0);
% to rotate Q/U, actually need to rotate vectorially  
% take the gradient of the ra/dec coordinate field with respect to
% the l,b basis, only need one really
% l/b gradients should be constants, which in the limit of a
% continuous field are equal to 1
[del.rax,del.ray]=gradient(ra);
[del.decx,del.decy]=gradient(dec);

% theta is the angle between unit vectors of RA and dec in the l,b basis
% theta is an array of vector rotation angles to go from l/b to ra/dec, 
% reverse sign to go the other way.  only theta.ra is useful, the rest
% are consistency checks.
theta.ra=atan2(del.ray,del.rax);
theta.dec=atan2(del.decy,del.decx);

% interpolate to find points relevant to ra/dec maps
% meshgrid points in RA,dec space
[x_tic,y_tic]=meshgrid(m.x_tic,m.y_tic);
% griddata takes data theta.ra, which is nonuniformly gridded in ra
% and dec, and re-grids it to x_tic,y_tic meshgrids using interpolation
theta.ra2=griddata(ra,dec,theta.ra,x_tic,y_tic);

for i=1:size(map,1)
  
  mapg(i,1).x_tic=mg.x_tic;
  mapg(i,1).y_tic=mg.x_tic;
  
  % just rotate the intensity map, no problems there, done in pixel space
  [mapg(i,1).T]=interp2(map(i).x_tic,map(i).y_tic,map(i).T,ra,dec);
  %[mapg(i,1).Tvar]=interp2(map(i).x_tic,map(i).y_tic,map(i).Tvar,ra,dec);
  %[mapg(i,1).Titime]=interp2(map(i).x_tic,map(i).y_tic,map(i).Titime,ra,dec);

  % from Stokes parameter definitions with no phase diff and zero V
  % U/Q=tan(2*chi), where chi is the rotation angle to the meas axis
  % and I^2=Q^2+U^2
  
  % Make pol vectors from Q,U maps
  % t = polarization angle
  % r = polarization magnitude
  [t,r]=qu2linpol(map(i).Q,map(i).U);
  
  t=theta.ra2+t;
  [mapg(i,1).Q]=interp2nan(map(i).x_tic,map(i).y_tic,r.*cos(2*t),ra,dec);
  [mapg(i,1).U]=interp2nan(map(i).x_tic,map(i).y_tic,r.*sin(2*t),ra,dec);
  
end

return

%%%%%%%%%%%%%%%%%%%%%%%
function mg=setmg(m,scale)

% preserve pixel sizes
mg.pixsize=m.pixsize;

% make gal map slightly larger to fit rotation
mg.xdos=scale*m.xdos; mg.ydos=scale*m.ydos;
% factor to subtract from low and high points to extend coords
dx=(scale-1)*m.xdos/2; dy=(scale-1)*m.ydos/2;

% field center
[cx,cy]=euler((m.hx+m.lx)/2,(m.hy+m.ly)/2,1,0);
if(cx>=180)
  cx=cx-360;
end

% subtract scaled factor, 
% add pixsize to l point to make even number of tics
mg.lx=cx-(mg.xdos/2)+mg.pixsize; mg.hx=cx+(mg.xdos/2);
mg.ly=cy-(mg.ydos/2)+mg.pixsize; mg.hy=cy+(mg.ydos/2);

% set coord axis
mg.x_tic=mg.lx:mg.pixsize:mg.hx; 
mg.y_tic=mg.ly:mg.pixsize:mg.hy;  

% calculate number of pixels
mg.nx=length(mg.x_tic); mg.ny=length(mg.y_tic);

return

%%%%%%%%%%%%%%%%%%%%%%%
function [m,map,mapopt,scale]=setdef(vai,nai)

switch nai  
  case 0
    load('maps_rbf/galcen_rbf_060718_filtm_weight0_jack0.mat')
  case 1    
    file=vai{1};
    load(file);
  case 2
    if(isstr(vai{1}))
      file=vai{1};
      load(file);      
      scale=vai{2};
    else
      m=vai{1};
      map=vai{2};
    end
  case 3
    m=vai{1};
    map=vai{2};
    mapopt=vai{3};
  case 4
    m=vai{1};
    map=vai{2};
    mapopt=vai{3};
    scale=vai{4};
end
    
if(~exist('mapopt','var'))
  mapopt=[];
end
if(~exist('scale','var'))
  scale=1.4;
end

return

