function im=f2i(ad,ft)
% im=f2i(ad,ft)
%
% Clem's idea of what a ft should do.
% Include power conserving normalization.
% Make output purely real (assuming input array is invariant
% under 180 deg rotation).
%
% It is Num Rec and TMS convention to use positive exponent
% in the ft operation that goes from frequency to image space.
% This corresponds to Matlab fft function.
%
% If input 3D treats as stack of images over 3rd dimension.
%
% Generalized to work for non-square array

ft=squeeze(ft);
n=ndims(ft);
if(n<2|n>3)
  error('Input must be 2D or 3D');
end

if(length(ad.del_u)>1)
  del_u=ad.del_u;
else
  del_u=[ad.del_u,ad.del_u];
end

% Do the transform and normalize
for i=1:size(ft,3)
  % how it was until 3/18:
  %im(:,:,i)=prod(del_u)*ifftshift(fft2(ifftshift(ft(:,:,i))));
  
  % this seems to be correct (only different for odd sized arrays)
  im(:,:,i)=prod(del_u)*fftshift(fft2(ifftshift(ft(:,:,i))));
end

% Check result purely real and disguard any imag component
if(sum(rvec(abs(imag(im))))./sum(rvec(abs(real(im))))>0.01)
  warning('Significant imag components thrown away');
end  
im=real(im);
return
