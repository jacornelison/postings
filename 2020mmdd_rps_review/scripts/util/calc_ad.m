function ad=calc_ad(Field_size_deg,N_pix);
% ad=calc_ad(Field_size_deg,N_pix);
%
% Make a structure containing all the axis data for an image
% in both image and fourier planes
%
% NB: Input Field_size in DEGREES
%
% Ouput structure has the following elements:
%
% Field_size_deg  (input)
% N_pix           (input)
% Field_size      (radians)
% t_del u_del     (spacings)
% t_val u_val     (axis values in radians)
% t_r u_r         (grids of radial values)
% t_val_deg       (useful for plots)

% Keep inputs in output structure too
ad.Field_size_deg=Field_size_deg;
ad.N_pix=N_pix;

% Most things in radians
ad.Field_size=Field_size_deg*(pi/180);

% Calc the spacing in t
ad.del_t=ad.Field_size/N_pix;

% Calc the spacing in u
ad.del_u=1/(ad.del_t*N_pix);

% The u values at pixel centers along the u axis
% NB: The zero value is in the N_pix/2+1 element of this even length array
ad.u_val=-ad.del_u*N_pix/2:ad.del_u:ad.del_u*(N_pix/2-1);

% The x/y positions at pixel centers
ad.t_val=-(N_pix/2)*ad.del_t:ad.del_t:(N_pix/2-1)*ad.del_t;

% Make grids of radial values in image and fourier planes
[x,y]=meshgrid(ad.t_val,ad.t_val); ad.t_r=sqrt(x.^2+y.^2);
[x,y]=meshgrid(ad.u_val,ad.u_val); ad.u_r=sqrt(x.^2+y.^2);

% x/y positions in degrees and minutes are useful for plots
ad.t_val_deg=ad.t_val*(180/pi);
ad.t_val_min=ad.t_val_deg*60;

return
