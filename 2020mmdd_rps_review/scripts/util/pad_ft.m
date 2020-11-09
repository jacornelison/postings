function [mo,ft]=pad_ft(m,ft,p,padval)
% [mo,ft]=pad_ft(m,ft,p)
%
% Pad the fts to increase map space resolution
% p is the size of the new FT arrays
% if the elements of p are all integer multiples of the
% corresponding dimensions of the input maps then the reverse FT
% will yield a superset of the original map points
% m is recalculated to give the proper values for the upsampled map
% - note that the "box on the sky" shifts around due to padding
% (although it stays the same size).
%
% e.g: 
% load maps/1459/real_aabd_filtp3_weight3_gs_dp1100_jack0.mat
% map=make_map(ac,m);
% map=add_masks(m,map);
% % set the max of mask to 1 as is done in deconv_map
% for i=1:numel(map)
%   map(i).Pw=map(i).Pw./max(map(i).Pw(:));
% end
% ad=calc_ad2([m.xdos,m.ydos],[m.nx,m.ny]);
% % calc the FTs as normal
% for i=1:numel(map)
%   ft(i)=calc_map_fts(map(i),ad,map(i).Pw,[2],'normal',false);
% end
% % pad them by factor 3
% n=[3,3]; [mp,ftp]=pad_ft(m,ft,[m.ny,m.nx].*n);
% % go back to map space
% for i=1:numel(map)
%   mapp(i).Q=f2i(ad,ftp(i).Q);
% end
% % plot a row of the interpolated map as red x's
% plot(mp.x_tic,mapp(2).Q(148,:),'rx-');
% hold on
% % overplot the corresponding row of the original map
% plot(m.x_tic,map(2).Q(50,:).*map(2).Pw(50,:),'b.');
% hold off

if(~exist('padval','var'))
  padval=0;
end

for i=1:numel(ft)
  if isfield(ft,'T')
    s=(p-size(ft(i).T))/2;
    ft(i).T=padarray2(ft(i).T,s,padval);
  end
  
  if isfield(ft,'Q')
    s=(p-size(ft(i).Q))/2;      
    ft(i).Q=padarray2(ft(i).Q,s,padval);
    ft(i).U=padarray2(ft(i).U,s,padval);
    ft(i).E=padarray2(ft(i).E,s,padval);
    ft(i).B=padarray2(ft(i).B,s,padval);
  end
end

% if input m empty can't recalc
if(isempty(m))
  mo=[];
  return
end

% pixel number increases simply
mo.nx=p(2); mo.ny=p(1);

% pixsize could now be non-square!
mo.pixsize(1)=m.pixsize*m.ny/mo.ny;
mo.pixsize(2)=m.pixsize*m.nx/mo.nx;

% padding always leaves the 1,1 map element unchanged
xs=m.x_tic(2)-m.x_tic(1); xs=xs*m.nx/mo.nx;
mo.x_tic=m.x_tic(1)+[0:(mo.nx-1)]*xs;
ys=m.y_tic(2)-m.y_tic(1); ys=ys*m.ny/mo.ny;
mo.y_tic=m.y_tic(1)+[0:(mo.ny-1)]*ys;

% edges of box are just pixel centers -/+ half the spacing
% - this means although it stays the same size it shifts
% around on the sky due to padding!
mo.lx=mo.x_tic(1)-xs/2;
mo.hx=mo.x_tic(end)+xs/2;
mo.ly=mo.y_tic(1)-ys/2;
mo.hy=mo.y_tic(end)+ys/2;

% degrees on sky in same way as get_map_defn 
c=m.xdos/(m.hx-m.lx);
mo.xdos=(mo.hx-mo.lx)*c;
mo.ydos=mo.hy-mo.ly;

mo.proj=m.proj;

return

%%%%%%%%%%%%%%%%%%%%%%%%%
function x=padarray2(x,s,padval)

% For an even numbered size use the trick from interpft of
% splitting the original Nyquist freq elements to maintain the
% necessary symmetry such that the FT is purely real

switch sum(mod(size(x),2)==0)
  case 0
    error('odd dims - not sure if this useful - see comments');
  
    % padding odd dim array does not seem to give a map
    % whose pixels have the superset property which makes this exercise
    % useful - if you find yourself here you will need to think more
    % about it - see interft2.m
    
    % both initial dim odd - pad symmetrically
    x=padarray(x,s,padval,'pre');
    x=padarray(x,s,padval,'post');
    
  case 1
    error('Mixed odd/even dim - if you really want this case expand pad_ft to handle it');
    
  case 2
    % both initial dim even
    % - split the previous Nyquist elements appropriately
    % the rules below were determined empirically by looking at the
    % output of interft2.m function
    x(1,2:end)=x(1,2:end)/2;
    x(end+1,2:end)=x(1,2:end);
    x(2:end-1,1)=x(2:end-1,1)/2;
    x(2:end-1,end+1)=x(2:end-1,1);
    x(1,1)=x(1,1)/4;
    x(end,1)=x(1,1);
    x(1,end)=x(1,1);
    x(end,end)=x(1,1);
    
    x=padarray(x,s,padval,'pre');
    x=padarray(x,s-1,padval,'post');
end

return
