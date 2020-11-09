function beam = gen_egauss_pair(ad,sigma,dparam)
% beam = gen_egauss_pair(ad,sigma,dparam)
% 
% Given differential Gaussian parameters, generate an A, B, sum, and diff
% beam in both image and Fourier planes.
% At the moment only generates a differentially elliptic beam centered
% at zero -- need to modify to add nonzero average ellipticity
%
% INPUTS:
%          ad:     standard ad struct on which to generate beam 
%          sigma:  nominal beamwidth, degrees
%          dparam: vector containing 
%                   [dsig,    dx,     dy,     dp,       dc]
%                    degrees, arcmin, arcmin, unitless, unitless
%                  that the A/B beams should have
%
% OUTPUTS:
%          beam: struct containing A/B/S (sum)/D (diff) beams in image
%                and fourier planes

% Parse dparam
dsig = dparam(1);
dx = dparam(2);
dy = dparam(3);
dp = dparam(4);
dc = dparam(5);

% Look in e.g. ffbm_makesetofegauss for reference
[X,Y] = meshgrid(ad.t_val_deg{1},ad.t_val_deg{2});

% dx/dy are traditionally in arcminutes
a.x = (dx/2)/60;
a.y = (dy/2)/60;
b.x = -(dx/2)/60;
b.y = -(dy/2)/60;
% Use sigma/dp/dc to get s/c/p per-detector
a.s = sigma + (dsig/2);
b.s = sigma - (dsig/2);
a.c = dc/2;
b.c = -dc/2;
a.p = dp/2;
b.p = -dp/2;

% Convert input diff params into a format that egauss2 can eat
[a.fwhm_maj a.fwhm_min a.theta] = egauss2_scp2mmt(a.s,a.c,a.p);
[b.fwhm_maj b.fwhm_min b.theta] = egauss2_scp2mmt(b.s,b.c,b.p);
% This lamely gives fwhm, convert back to sigma
% Also theta is in degrees but egauss2 wants radians
beamparam_a = [1 a.x a.y a.fwhm_maj/(2*sqrt(2*log(2))) ...
               a.fwhm_min/(2*sqrt(2*log(2))) a.theta*pi/180];
beamparam_b = [1 b.x b.y b.fwhm_maj/(2*sqrt(2*log(2))) ...
               b.fwhm_min/(2*sqrt(2*log(2))) b.theta*pi/180];

% Generate map-space Gaussians
A = egauss2(beamparam_a,X,Y);
B = egauss2(beamparam_b,X,Y);

% Image plane: integral normalize and make complex 
Ai = complex(A/sum(A(:)));
Bi = complex(B/sum(B(:)));
Si = (Ai + Bi)/2;   % sum
Di = (Ai - Bi)/2;   % diff

% Fourier plane: peak normalize
Af = ifftshift(ifft2(ifftshift(Ai)));
Af = Af/max(Af(:));
Bf = ifftshift(ifft2(ifftshift(Bi)));
Bf = Bf/max(Bf(:));
Sf = (Af + Bf)/2;
Df = (Af - Bf)/2;

beam(1).Ai = Ai;
beam(1).Bi = Bi;
beam(1).Si = Si;
beam(1).Di = Di;
beam(1).Af = Af;
beam(1).Bf = Bf;
beam(1).Sf = Sf;
beam(1).Df = Df;

return

