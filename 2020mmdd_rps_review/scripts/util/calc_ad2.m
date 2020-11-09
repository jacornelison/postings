function ad=calc_ad2(Field_size_deg,N_pix);
% ad=calc_ad2(Field_size_deg,N_pix);
%
% Version for non-square grid
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
ad.del_t=ad.Field_size./N_pix;

% Calc the spacing in u
ad.del_u=1./(ad.del_t.*N_pix);

% allow for odd sized dimensions
for i=1:length(N_pix)
  n=N_pix(i);
  if(rem(n,2)==0)
    m=-n/2:n/2-1;
  else
    m=-n/2+0.5:n/2-0.5;
  end
  ad.u_val{i}=m*ad.del_u(i);
  ad.t_val{i}=m*ad.del_t(i);
end

% Make grids of radial values in image and fourier planes
[x,y]=meshgrid(ad.t_val{1},ad.t_val{2}); ad.t_r=sqrt(x.^2+y.^2);
[x,y]=meshgrid(ad.u_val{1},ad.u_val{2}); ad.u_r=sqrt(x.^2+y.^2);

% x/y positions in degrees and minutes are useful for plots
ad.t_val_deg{1}=ad.t_val{1}*(180/pi); ad.t_val_deg{2}=ad.t_val{2}*(180/pi);
ad.t_val_min{1}=ad.t_val_deg{1}*60; ad.t_val_min{2}=ad.t_val_deg{2}*60;

return
