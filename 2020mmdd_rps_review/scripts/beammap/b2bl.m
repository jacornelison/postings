function [l,B_l] = b2bl(ad,beam,nonorm)
% function [l,B_l] = b2bl(ad,beam,nonorm)
%
% Turn a beam into a b_l profile.  Does not cut off at the MTF or
% correct for the chopper aperture.
%
% INPUTS
%  
%   ad:           standard ad map struct
%   beam:         beam map matching ad (size ad.N_pix)
%                 can be 1D profile (generated from get_beamprofile)
%                 or 2D array
%   nonorm:       1 -> don't normalize to l(0) = 1
%                 default 0

if nargin < 3
  nonorm = 0;
end

% Make 2D ad struct if input beam is 1D
beam = squeeze(beam);
if isvector(beam)
  ad = calc_ad2([ad.Field_size_deg,ad.Field_size_deg],[ad.N_pix,ad.N_pix]);
end

n = floor(ad.N_pix(1)/2+1);
l = ad.u_val{1}(n:end)*2*pi;
li = 0:max(l);

if all(isnan(beam(:)))
  l = cvec(li);
  B_l = cvec(NaN(size(l)));
  return
end

% If beam is 2D, make 1D profile first
if ~isvector(beam)
  [r,profile] = get_beamprofile(ad,beam,'median','none');
else
  profile = beam;
  rmax = min(abs([ad.t_val_deg{1}(1),ad.t_val_deg{1}(end),...
	          ad.t_val_deg{2}(1),ad.t_val_deg{2}(end)]));
  r = linspace(0,rmax,round(ad.N_pix(1)/2) + 1);
end

% Convert 1-D beam profile into a 2-D beam (effectively symmetrizing it)
beam2d = beam_1d_to_2d(ad,r,profile);

% FFT this beam
bl2d = abs(i2f(ad,beam2d));

% Get 1-D B_ls
B_l = bl2d(n,n:end);

% Must have a Bl for each ell
li = 0:max(l);
B_li = interp1(l,B_l,li,'linear');

l = li;
B_l = B_li;

% Normalize - l=0 is 1 by default
if nonorm
  B_l = cvec(B_l);
  l = cvec(l);
  return
else
    if l(1) == 0
      B_l = B_l/B_l(1);
    end
return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function beam2d = beam_1d_to_2d(ad,r,profile)
  z = interp1(r,profile,ad.t_r(:)*180/pi,'linear');
  beam2d = reshape(z,size(ad.t_r));
  beam2d(ad.t_r*180/pi > max(ad.t_val_deg{1})) = 0;
return






