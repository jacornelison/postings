function [r,b]=bl2beam(ad,l,Bl)
% [r,b]=bl2beam(ad,l,Bl)
%
% Convert beam window function l,Bl into real space 1d
% beam function r,b
%
% e.g.:
% [l,Bl]=get_bl('B2bbns');
% ad=calc_ad(20,512);
% [r,b]=bl2beam(ad,l,Bl);
% semilogy(r,b); grid

% interp Bl onto grid
ft=interp1(l,Bl,ad.u_r*2*pi,[],'extrap');
ft=reshape(ft,size(ad.u_r));

% do the ft
im=f2i(ad,ft);

% pull out the 1d profile
n=ad.N_pix/2+1;
r=ad.t_val_deg(n:end);
b=im(n,n:end);

% normalize to unit peak height
b=b/b(1);

return
