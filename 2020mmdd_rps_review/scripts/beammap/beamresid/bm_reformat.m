function [beam split] = bm_reformat(A,B,Asplit,Bsplit)
% [beam split] = bm_reformat(A,B,Asplit,Bsplit)
%
% Get a composite beam (and split) into a format for beamresid work
% Produces integral-normalized image-plane and peak-normalized
% image-plane A/B/Sum/Diff beams
% 
% INPUTS
%          A/B/: beam maps, any size/shape (as long as ifft2 works)
%                Both must be present to make pairsum and diff
%          Asplit/Bsplit: identical operations on split maps if requested
%                         (empty okay too)
%                         But uses non-split normalization!!!
%
% OUTPUTS
%          beam/split: a single struct containing 'A/B/S/D'.i/f
%                      i.e. A, B, Sum, Diff in image and fourier planes

A_norm = sum(A(:));
B_norm = sum(B(:));

Ai = complex(A/A_norm);
Bi = complex(B/B_norm);
Si = (Ai + Bi)/2;   % sum
Di = (Ai - Bi)/2;   % diff

if ~isempty(Asplit)
  Asplit_i = complex(Asplit/A_norm);
  Bsplit_i = complex(Bsplit/B_norm);
  Ssplit_i = (Asplit_i + Bsplit_i)/2;
  Dsplit_i = (Asplit_i - Bsplit_i)/2;
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
if ~isempty(Asplit)
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

if ~isempty(Asplit)
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

