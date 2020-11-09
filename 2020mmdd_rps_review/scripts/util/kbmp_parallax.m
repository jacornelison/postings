function [az, el, pa] = kbmp_parallax(pos, pnt, ort, source, p, varargin)
% [az, el, pa] = kbmp_parallax(pos, pnt, ort, source, p)
%
% Calculates parallax corrected pointing for all pixels at a specified 
% source distance.
%
% First, individual pointing and orientation vectors are calculated for each 
% pixel. Next, the parallax correction is calculated by extending these 
% pointing vectors from the physical aperture location out to their 
% intersection with a sphere that is centered at the origin and has radius 
% equal to the source distance.
%
% [Arguments]
%   The first three input arguments should be arrays with shape [N,3], where 
%   the first dimension is the number of samples and the second dimension 
%   are (x,y,z) horizon coordinates. These should be obtained from 
%   kbmp_reflect.m for beam mapping with the far field flat (or directly 
%   from kbmp_mount.m if no far field flat).
%
%   pos     Position vectors for the telescope aperture.
%   pnt     Boresight pointing vectors.
%   ort     Boresight orientation vectors.
%   source  Data structure containing source position parameters.
%   p       Array data structure, from get_array_info.
%
%   Fields for source data structure:
%     source.distance  Distance, in meters, to the source.
%     source.azimuth   Azimuth, in degrees, of source as seen from Keck
%     source.height    Height, in meters, of source above floor of ground 
%                      screen
%     source.lateral_distance  Horizontal distance, in meters, between source
%                              and azimuth axis
%
% [Optional keywords]
%   'NoMirror'      By default, this function assumes that the input 
%                   boresight and orientation vectors have been reflected 
%                   off the far field flat, so the pixel positions and 
%                   orientations should be constructed with a mirror 
%                   reflection. Include this optional keyboard to get pixel 
%                   positions and orientations with normal (non-reflected) 
%                   parity.
%   'BeamCentered'  By default, this function returns apparent azimuth and 
%                   elevation. Including this optional keyword instead makes 
%                   it return coordinates in a spherical coordinate system 
%                   centered on the expected (from the mount and pointing 
%                   models) position of the main beam of each pixel. 
%                   The coordinates are:
%                     phi    Azimuthal coordinate defined so that phi = 0 is 
%                            in the direction of increasing elevation when 
%                            dk = 0.
%                     theta  Angular distance from the main beam.
%                     pa     (To be determined).
%   'BoresightCentered'  By default, this function applies pixel pointing 
%                        offsets. Including this optional keyword causes 
%                        these offsets to be ignored and the boresight 
%                        pointing is used. Thus, specifying both 
%                        'BeamCentered' and 'BoresightCentered' results in 
%                        the (r, theta) coordinate system centered on the 
%                        boresight
%		
% [Returns]
%   All three outputs have shape [N,M], where N is the number of samples 
%   and M is the number of pixels.
%
%   az  Corrected azimuth for each pixel. 
%   el  Corrected elevation for each pixel.
%   pa  Corrected parallactic angle for each pixel.

% Last update: 2014-02-08 CAB

% Parse optional arguments.
use_mirror = 1;
for i=1:length(varargin)
  if strcmp(varargin{i}, 'NoMirror')
    use_mirror = 0;
  end
end

beam_centered = 0;
for i=1:length(varargin)
  if strcmp(varargin{i}, 'BeamCentered')
    beam_centered = 1;
  end
end

boresight_centered = 0;
for i=1:length(varargin)
  if strcmp(varargin{i}, 'BoresightCentered')
    boresight_centered = 1;
  end
end


% Check size of 3-vector input arguments.
nsamp = size(pos, 1);
if (size(pos, 2) ~= 3)
  disp('[kbmp_parallax] ERROR: argument pos must have size [N,3]');
  return
end
if (size(pnt, 1) ~= nsamp) | (size(pnt, 2) ~= 3)
  disp('[kbmp_parallax] ERROR: argument pnt must have size [N,3]');
  return
end
if (size(ort, 1) ~= nsamp) | (size(ort, 2) ~= 3)
  disp('[kbmp_parallax] ERROR: argument ort must have size [N,3]');
  return
end

% Calculate the cross-product of boresight pointing and
% orientation. This will serve as the third axis for an orthonormal
% coordinate system centered around the boresight.
% Flip the sign of this vector depending on desired parity.
if use_mirror
  % Left-handed coordinates:
  %   pnt = zhat, ort = xhat, yhat = xhat CROSS zhat
  ort2 = cross(ort, pnt);
else
  % Right-handed coordinates:
  %   pnt = zhat, ort = xhat, yhat = xhat CROSS zhat
  ort2 = cross(pnt, ort);
end

% Allocate output arrays.
npix = length(p.r);
az = zeros(nsamp, npix);
el = zeros(nsamp, npix);
pa = zeros(nsamp, npix);

% Precalculate cos and sin factors.
dtor = pi / 180;
cr = cos(dtor * p.r);
sr = sin(dtor * p.r);
ct = cos(dtor * p.theta);
st = sin(dtor * p.theta);

if beam_centered
  %Calculate vector from aperture to source
  %First, calculate source position
  source_position = [ source.lateral_distance * cos(-source.azimuth *dtor), ...
    source.lateral_distance * sin(-source.azimuth *dtor), ...
    source.height ];
  aperture_source_vector = repmat(source_position, nsamp, 1) - pos; 
  aperture_source_unit_vector = aperture_source_vector ./ ...
    repmat(sqrt(dot(aperture_source_vector, aperture_source_vector, 2)), 1, 3);
end

% Loop over channels.
for i=1:npix
  if boresight_centered
    %Use boresight pointing as channel pointing
    apnt = pnt;
    aort = ort;
  else
  % Calculate pointing and orientation for this channel.
    apnt = sr(i) * ct(i) * ort - sr(i) * st(i) * ort2 + cr(i) * pnt;
    aort = (cr(i) * ct(i)^2 + st(i)^2) * ort + ...
         ((1 - cr(i)) * ct(i) * st(i)) * ort2 + ...
         (-1 * sr(i) * ct(i)) * pnt;
  end

  if beam_centered %beam centered coordinates
    aort2 = cross(apnt, aort, 2); %3rd Orthonormal vector
    %theta from pointing dot aperture - source vector
    %apnt and aperture_source_unit_vector have unit norm
    el(:, i) = acos(dot(apnt, aperture_source_unit_vector, 2)) / dtor;
    %phi from projection of aperture - source vector onto perpendicular plane
    az(:, i) = atan2(dot(aort2, aperture_source_unit_vector, 2), ...
      dot(aort, aperture_source_unit_vector, 2)) /dtor; 
    %TODO: pa
  else %normal coordinates 

    % Calculate the distance from the aperture to a sphere of radius
    % equal to the source distance, along the pointing vector.
    pos_pos = sum(pos.^2, 2);
    pos_pnt = sum(pos .* apnt, 2);
    scl = sqrt(pos_pnt.^2 - pos_pos + source.distance^2) - pos_pnt;

    % Calculate azimuth and elevation of the intersection of the
    % pointing vector with the sphere.
    for j=1:3
      bpnt(:,j) = pos(:,j) + scl .* apnt(:,j);
    end
    az(:,i) = atan2(-1 * bpnt(:,2), bpnt(:,1)) / dtor;
    el(:,i) = asin(bpnt(:,3) / source.distance) / dtor;

    % Calculate parallactic angle of the orientation vector at this
    % azimuth and elevation.
    ea(:,1) = cos(dtor * az(:,i)) .* sin(dtor * el(:,i));
    ea(:,2) = -1 * sin(dtor * az(:,i)) .* sin(dtor * el(:,i));
    ea(:,3) = -1 * cos(dtor * el(:,i));
    eb(:,1) = sin(dtor * az(:,i));
    eb(:,2) = cos(dtor * az(:,i));
    eb(:,3) = zeros(size(az(:,i)));
    pa(:,i) = atan2(sum(aort .* eb, 2), sum(aort .* ea, 2)) / dtor;
  end


  % Convert this parallactic angle to pipeline convention.
  % 1. Pipeline uses IAU convention for parallactic angle, not
  %    Healpix convention. They differ by a minus sign.
  % 2. Add p.chi (detector angle) and p.thetaref (reference angle
  %    that chi is defined with respect to).
  % See 2012-12-07 posting.
  pa(:,i) = -1 * pa(:,i) + p.chi(i) + p.chi_thetaref(i);
  % Wrap parallactic angle range to [-90, 90].
  while max(pa(:,i)) > 90
      pa(pa(:,i) > 90,i) = pa(pa(:,i) > 90,i) - 180;
  end
  while min(pa(:,i)) < -90
      pa(pa(:,i) < -90,i) = pa(pa(:,i) < -90,i) + 180;
  end
end
