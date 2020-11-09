function [pos pnt ort] = kbmp_mount(az, el, dk, mount, drumangle)
% [pos pnt ort] = kbmp_mount(az, el, dk, mount, drumangle)
%
% Calculate 3-D position and pointing for Keck mount.
%
% This function is usually not called directly. Use the wrapper function 
% keck_beam_map_pointing.m.
%
% [Arguments]
%   The input encoder arguments (az, el, dk) can be arrays
%   (i.e. timestreams), but they should all have matching length.
%
%   az         Azimuth encoder position, in degrees.
%   el         Elevation encoder position, in degrees.
%   dk         Deck encoder position, in degrees.
%   mount      Data structure containing mount parameters.
%   drumangle  Receiver drum angle, specified separately from the mount 
%              data structure.
%
%   Fields for mount data structure:
%     mount.aperture_offr  Distance between the center of the telescope 
%                          aperture and the drum axis, in meters. 
%                          Default = 0.5458 meters.
%     mount.aperture_offz  Position of the telescope aperture along the 
%                          drum axis, in meters. Default = 1.5964 meters. 
%     mount.dk_offx        Offset between the elevation and dk axes, in 
%                          meters. Default = -1.0196 meters.
%     mount.dk_offy        Offset between the azimuth and dk axes, in 
%                          meters. Default = 0 meters.
%     mount.el_tilt        Standard mount model parameter for tilt of the 
%                          elevation axis, in degrees. Default = 0 degrees. 
%     mount.el_offx        Offset between the azimuth and elevation axes, 
%                          in meters. Default = 0 meters.
%     mount.el_offz        Height of the elevation axis above the base of 
%                          the groundscreen, in meters. Default = 1.1750 
%                          meters.
%     mount.az_tilt_ha     Standard mount model parameter for tilt of the 
%                          azimuth axis in the direction of increasing hour 
%                          angle, in degrees. Default = 0 degrees. 
%     mount.az_tilt_lat    Standard mount model parameter for tilt of the 
%                          azimuth axis in the direction of increasing 
%                          latitude, in degrees. Default = 0 degrees.
%
%   NOTE on mount data structure: While the code supports azimuth and 
%   elevation tilts, I recommend using invpointing_model.m to correct the 
%   mount pointing before passing it to this function and then setting the 
%   tilt parameters to zero in the mount data structure.
%
% [Returns]
%   All three output arguments are 2-D arrays with shape [N,3],
%   where N is the number of samples for the input encoder
%   arguments. Each output sample is a 3-vector (x,y,z) in
%   horizon coordinates.
%
%   pos     Position vector for the telescope aperture, relative to a point 
%           on the azimuth axis at the height of the groundscreen floor. 
%           The units of this vector are meters.
%   pnt     Unit vector for telescope boresight pointing.
%   ort     Unit vector for telescope orientation (perpendicular to pnt).

% Last update: 2014-02-08 CAB

% Check length of input az/el/dk arguments.
nsamp = length(az);
if (length(el) ~= nsamp) | (length(dk) ~= nsamp)
  disp('[kbmp_mount] ERROR: size of input encoder arrays do not match');
  return;
end

% Initialize vectors.
pos = zeros(nsamp, 3);
pnt = repmat([0, 0, 1], nsamp, 1);
ort = repmat([1, 0, 0], nsamp, 1);

% Translate position by aperture_offr in the +x direction.
pos(:,1) = pos(:,1) + mount.aperture_offr;

% Rotate about the z-axis by dk + drum_angle - 90.
euler = zeros(nsamp, 3);
euler(:,3) = dk + drumangle - 90;
pos = rotate_3d(pos, euler);
% Don't include drum angle in orientation vector rotation, because
% it's already included in p.theta from get_array_info().
% The pointing vector is aligned with the z-axis, so this rotation
% leaves it unchanged.
euler(:,3) = dk - 90;
pnt = rotate_3d(pnt, euler);
ort = rotate_3d(ort, euler);

% Translate position by aperture_offz in the +z direction.
pos(:,3) = pos(:,3) + mount.aperture_offz;

% Translate position by dk_offx in the +x direction.
pos(:,1) = pos(:,1) + mount.dk_offx;

% Translate position by dk_offy in the +y direction.
% Included for completeness -- this parameter is zero for Keck.
pos(:,2) = pos(:,2) + mount.dk_offy;

% Rotate about the y-axis by 90 - el.
euler = zeros(nsamp, 3);
euler(:,2) = 90 - el;
pos = rotate_3d(pos, euler);
pnt = rotate_3d(pnt, euler);
ort = rotate_3d(ort, euler);

% Elevation non-perpendicularity correction.
% Rotate about the x-axis by el_tilt.
euler = [-90, mount.el_tilt, 90];
pos = rotate_3d(pos, euler);
pnt = rotate_3d(pnt, euler);
ort = rotate_3d(ort, euler);

% Translate by el_offx in the +x direction.
% Included for completeness -- this parameter is zero for Keck.
pos(:,1) = pos(:,1) + mount.el_offx;

% Translate by el_offz in the +z direction.
pos(:,3) = pos(:,3) + mount.el_offz;

% Rotate about the z-axis by -az.
euler = zeros(nsamp, 3);
euler(:,3) = -1 * az;
pos = rotate_3d(pos, euler);
pnt = rotate_3d(pnt, euler);
ort = rotate_3d(ort, euler);

% Azimuth tilt correction.
% Convert from tilt in hour angle and latitude directions to
% bearing and magnitude of the az axis tilt. 
% Bearing follows azimuth convention (i.e. left-handed).
az_tilt_bearing = atan2(-1 * sin(mount.az_tilt_ha * pi / 180), ...
    sin(mount.az_tilt_lat * pi / 180)) * 180 / pi;
az_tilt_magnitude = asin(sqrt(sin(mount.az_tilt_ha * pi / 180)^2 + ...
    sin(mount.az_tilt_lat * pi / 180)^2)) * 180 / pi;
% Now apply the rotation.
euler = [az_tilt_bearing, az_tilt_magnitude, -1 * az_tilt_bearing];
pos = rotate_3d(pos, euler);
pnt = rotate_3d(pnt, euler);
ort = rotate_3d(ort, euler);
