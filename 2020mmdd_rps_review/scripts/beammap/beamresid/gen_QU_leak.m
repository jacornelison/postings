function [Q U] = gen_QU_leak(ad,T,beam,angles,type)
% [Q U] = gen_QU_leak(T,beam,angles,type)
%
% Given a T map and an A/B beam, generate the leaked Q/U maps
% Operates in both image and fourier planes
%
% INPUTS:
%         ad:   ad struct for T map
%         T:    input temperature map
%         beam: struct containing A/B/S/D 'i' and 'f' to be
%               convolved/multiplied with T
%         angles: set of dk angles to loop over (radians)
%         type: 'image','fourier'
%               Only 'image' guaranteed to work right now
% 
% OUTPUTS:
%         Q/U: maps (corresponding to ad) with leaked T->P

% Ensure that all points will be rotated correctly and repeatably with
% imrotate: kill any points which are more than half the length of the
% window away from center.  Of course, doing this in the F-plane will affect
% high-ell points, but for image-plane low-ell points.
indtokill = find(ad.u_r > max(ad.u_val{1})); % needs to be square
T.fp(indtokill) = 0;
beam.Af(indtokill) = 0;
beam.Bf(indtokill) = 0;
indtokill2 = find(ad.t_r > max(ad.t_val{1}));
T.im(indtokill) = 0;
beam.Ai(indtokill) = 0;
beam.Bi(indtokill) = 0;

dcos = zeros(size(T.fp));
dsin = zeros(size(T.fp));
Q_recov = zeros(size(T.fp));
U_recov = zeros(size(T.fp));
cc = 0;
ss = 0;
cs = 0;

for ii = 1:length(angles)
  
  % Rotate beam and do convolution/multiplication
  switch type
    case 'image'
      %rotbeam(ii).A = imrotate(beam.Ai,angles(ii)*180/pi,'crop');
      %rotbeam(ii).B = imrotate(beam.Bi,angles(ii)*180/pi,'crop');   
      rotbeam(ii).D = imrotate(beam.Di,angles(ii)*180/pi,'crop');

      % Should be pure real anyway...
      %T_A = conv2(real(T.im),rotbeam(ii).A,'same');
      %T_B = conv2(real(T.im),rotbeam(ii).B,'same');
      T_D = conv2(real(T.im),rotbeam(ii).D,'same');
      
      diff_T = T_D; % this appears to work...

    case 'fourier'
      rotbeam(ii).A = imrotate(beam.Af,angles(ii)*180/pi,'crop');
      rotbeam(ii).B = imrotate(beam.Bf,angles(ii)*180/pi,'crop');
      
      T_A = (T.fp).*rotbeam(ii).A;
      T_B = (T.fp).*rotbeam(ii).B;

      diff_T = (T_A - T_B)/2;
  end
  
  dcos = dcos + cos(2*angles(ii))*diff_T;
  dsin = dsin + sin(2*angles(ii))*diff_T;
  cc = cc + (cos(2*angles(ii)))^2;
  ss = ss + (sin(2*angles(ii)))^2;
  cs = cs + (cos(2*angles(ii))*sin(2*angles(ii)));

end

% World's stupidest matrix inversion
for ii = 1:size(T.fp,1)
  for jj = 1:size(T.fp,2)
      
    diffmat = [dcos(ii,jj); dsin(ii,jj)];
    obs = [cc cs; cs ss];
    recov = inv(obs)*diffmat;
    Q_recov(ii,jj) = recov(1);
    U_recov(ii,jj) = recov(2);
    
  end
end

% Return both image and f-plane Q/U maps
switch type
  case 'image'
    Q.im = Q_recov;
    U.im = U_recov;
    Q.fp = ifftshift(ifft2(ifftshift(Q_recov)));
    U.fp = ifftshift(ifft2(ifftshift(U_recov)));
    
  case 'fourier'
    Q.fp = Q_recov;
    U.fp = U_recov;
    Q.im = fftshift(fft2(fftshift(Q_recov)));
    U.im = fftshift(fft2(fftshift(U_recov)));
end

return

