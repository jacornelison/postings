function [mpos mnorm] = kbmp_mirror(az, el, mount, mirror)
% [mpos mnorm] = kbmp_mirror(az, el, mount, mirror)
%
% Calculates position and orientation of Keck far field flat mirror.
%
% This function is usually not called directly. Use the wrapper function 
% keck_beam_map_pointing.m
%
% [Arguments]
%   The input encoder arguments (az, el) can be arrays (i.e. timestreams), 
%   but they should all have matching length.
%
%   az      Azimuth encoder position, in degrees.
%   el      Elevation encoder position, in degrees.
%   mount   Data structure containing mount parameters. For details, see 
%           comments for kbpm_mount.m
%   mirror  Data structure containing mirror parameters.
%
%   Fields for mirror data structure:
%     mirror.height  Distance to the mirror along the drum axis, in meters, 
%                    starting from a position even with the telescope 
%                    apertures.
%     mirror.tilt    Angle between the mirror normal vector and the drum 
%                    axis, in degrees, along the direction perpendicular to 
%                    the elevation axis. Tilt = 0 degrees corresponds to the 
%                    mirror reflecting the beams straight back to the 
%                    telescopes. Positive values of tilt correspond to the 
%                    usual mirror orientation, which reflects the beams "off 
%                    the back" of the mount (i.e. in a direction opposite
%                    of azimuth).
%     mirror.roll    Angle between the mirror normal vector and the drum 
%                    axis, in degrees, along the direction parallel to the 
%                    elevation axis. The mirror should have roll = 0 if the 
%                    carbon fiber rods are adjusted evenly on both sides of 
%                    the hexapod.
%
% [Returns]
%   Both returned values are 2-D arrays with shape [N,3], where N is the 
%   number of samples for the input (az, el) arguments. Each output sample 
%   is a 3-vector (x,y,z) in horizon coordinates. 
%
%   mpos    Position vector to the intersection between the mirror and the 
%           drum axis, relative to a point on the azimuth axis at the height 
%           of the groundscreen floor. Units are meters.
%   mnorm   Unit vector normal to the mirror surface.

% Last update: 2014-02-08 CAB

% Check length of input az/el arguments.
nsamp = length(az);
if (length(el) ~= nsamp)
  disp('[kbmp_mirror] ERROR: size of input encoder arrays do not match');
  return;
end

% Use keck_mount.m to calculate mirror position.
% 1. Set aperture_offr to zero (get on the drum axis).
% 2. Add mirror_height to aperture_offz.
mnt = mount;
mnt.aperture_offr = 0;
mnt.aperture_offz = mnt.aperture_offz + mirror.height;
[mpos a b] = kbmp_mount(az, el, 0 * az, mnt, 0);

% NOTE: Mirror normal vector is calculated using mount az/el, but
% the calculation doesn't include az tilt and el tilt
% corrections. I don't think this matters much, but it could be
% added later.

% Initialize mirror normal vector.
% This unit vector is initialized with the opposite orientation to
% the boresight pointing vector. So if the mirror orientation
% angles are set to zero, this means that the mirror is reflecting
% right back at the telescope.
mnorm = repmat([0, 0, -1], nsamp, 1);

% Apply mirror roll.
euler = zeros(nsamp, 3);
euler(:,2) = mirror.roll;
euler(:,1) = 90;
mnorm = rotate_3d(mnorm, euler);

% Apply mirror tilt.
euler = zeros(nsamp, 3);
euler(:,2) = mirror.tilt;
mnorm = rotate_3d(mnorm, euler);

% Move mount to encoder az/el position.
euler = zeros(nsamp, 3);
euler(:,2) = 90 - el;
euler(:,1) = -1 * az;
mnorm = rotate_3d(mnorm, euler);
