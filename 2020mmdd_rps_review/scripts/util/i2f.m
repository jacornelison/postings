function ft=i2f(ad,im)
% ft=i2f(ad,im)
%
% Clem's idea of what a ft should do.
% Remove bogus (1/N^2) normalization and make power conserving.
%
% It is Num Rec and TMS convention to use negative exponent
% in the ft operation that goes from image to frequency space.
% This corresponds to Matlab ifft function.
%
% If input 3D treats as stack of images over 3rd dimension.
%
% Generalized to work for non-square array

im=squeeze(im);
n=ndims(im);
if(n<2|n>3)
  error('Input must be 2D or 3D');
end

if(length(ad.del_t)>1)
  del_t=ad.del_t;
else
  del_t=[ad.del_t,ad.del_t];
end

if(length(ad.N_pix)>1)
  N_pix=ad.N_pix;
else
  N_pix=[ad.N_pix,ad.N_pix];
end

for i=1:size(im,3)
  % how it was until 3/18:
  %ft(:,:,i)=prod(del_t)*prod(N_pix)*fftshift(ifft2(fftshift(im(:,:,i))));
  
  % this seems to be correct (only different for odd sized arrays)
  ft(:,:,i)=prod(del_t)*prod(N_pix)*fftshift(ifft2(ifftshift(im(:,:,i))));
end

return
