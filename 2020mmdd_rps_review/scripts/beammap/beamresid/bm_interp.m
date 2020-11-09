function [beam split] = bm_interp(ad_T,ad_beam,beam_in,split_in,scalediff)
% [beam split] = bm_interp(ad_T,ad_beam,beam_in,split_in,scalediff)
% 
% Interpolate a beam map to a new ad grid -- interpolation in image space
%
% INPUTS
%        ad_T: ad struct to interpolate TO
%        ad_beam: ad struct of input beam maps
%        beam_in: real beam
%                 should have fields 'im' and 'fp' for A/B/sum/diff
%        split_in: split associated with real, empty if doesn't exist
%        scalediff: 1 to preserve normalization when interpolating
%                   (highly recommended!)  Applied to split too
% 
% OUTPUTS
%        beam: rescaled version of beam_in
%        split: rescaled version of split_in

[X_T Y_T] = meshgrid(ad_T.t_val{1});
[X_B Y_B] = meshgrid(ad_beam.t_val{1});

%A = interp2(X_B,Y_B,real(beam_in.Ai),X_T,Y_T,'spline');
%B = interp2(X_B,Y_B,real(beam_in.Bi),X_T,Y_T,'spline');
% Set to zero anything outside of original beam grid
A = interp2(X_B,Y_B,real(beam_in.Ai),X_T,Y_T,'cubic',0);
B = interp2(X_B,Y_B,real(beam_in.Bi),X_T,Y_T,'cubic',0);

if ~isempty(split_in)
  Asplit = interp2(X_B,Y_B,real(split_in.Ai),X_T,Y_T,'cubic',0);
  Bsplit = interp2(X_B,Y_B,real(split_in.Bi),X_T,Y_T,'cubic',0);
end

scalefac_A = sum(A(:));
scalefac_B = sum(B(:));

Ai = complex(A/scalefac_A);
Bi = complex(B/scalefac_B);
Si = (Ai + Bi)/2;   % sum
%Di = (Ai - Bi)/2;   % diff
Di = interp2(X_B,Y_B,real(beam_in.Di),X_T,Y_T,'cubic',0);
if scalediff
  Di = Di/ ((scalefac_A + scalefac_B)/2);
end

if ~isempty(split_in)
  Asplit_i = complex(Asplit/scalefac_A);
  Bsplit_i = complex(Bsplit/scalefac_B);
  Ssplit_i = (Asplit_i + Bsplit_i)/2;
  Dsplit_i = interp2(X_B,Y_B,real(split_in.Di),X_T,Y_T,'cubic',0);
  if scalediff
    Dsplit_i = Dsplit_i / ((scalefac_A + scalefac_B)/2);
  end
end

% Fourier plane: peak normalize
Af = ifftshift(ifft2(ifftshift(Ai)));
Af_norm = max(Af(:));
Af = Af/Af_norm;
Bf = ifftshift(ifft2(ifftshift(Bi)));
Bf_norm = max(Bf(:));
Bf = Bf/Bf_norm;
Sf = (Af + Bf)/2;
Df = (Af - Bf)/2;

% Does this work?  Hell if I know...
if ~isempty(split_in)
  Asplit_f = ifftshift(ifft2(ifftshift(Asplit_i)));
  Asplit_f = Asplit_f/Af_norm;
  Bsplit_f = ifftshift(ifft2(ifftshift(Bsplit_i)));
  Bsplit_f = Bsplit_f/Bf_norm;
  Ssplit_f = (Asplit_f + Bsplit_f)/2;
  Dsplit_f = (Asplit_f - Bsplit_f)/2;
end

beam(1).Ai = Ai;
beam(1).Bi = Bi;
beam(1).Si = Si;
beam(1).Di = Di;
beam(1).Af = Af;
beam(1).Bf = Bf;
beam(1).Sf = Sf;
beam(1).Df = Df;

if ~isempty(split_in)
  split(1).Ai = Asplit_i;
  split(1).Bi = Bsplit_i;
  split(1).Si = Ssplit_i;
  split(1).Di = Dsplit_i;
  split(1).Af = Asplit_f;
  split(1).Bf = Bsplit_f;
  split(1).Sf = Ssplit_f;
  split(1).Df = Dsplit_f;
end

return