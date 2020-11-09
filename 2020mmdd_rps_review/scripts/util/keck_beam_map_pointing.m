function [r, theta, psi] = keck_beam_map_pointing(az, el, dk, ...
    mount, mirror, source, p, varargin)
% [r, theta, psi] = keck_beam_map_pointing(az, el, dk, mount, ...
%                                          mirror, source, p)
%
% Calculate pointing for Keck beam map analysis, including mount model,
% mirror reflection, and parallax effects. This function is a wrapper
% that calls kbmp_mount.m, kbmp_mirror.m, kbmp_reflect.m, and
% kbmp_parallax_sph.m (or kbmp_parallax.m).
%
% In normal operation, this function takes mount [az, el, dk] coordinates
% (which can be time series arrays) and returns [r, theta, psi] values
% specifying the source position relative to each detector. Detectors to
% calculate pointing for are specified by passing in the standard r and
% theta parameters (i.e. from get_array_info). The details of the pointing
% correction are controlled by the 'mount', 'mirror', and 'source'
% arguments, which are each data structures containing many fields (see
% below). This function will attempt to fill in any missing fields in the
% data structures using typical values for Keck. For the 'mount' structure
% in particular, the default values are the only ones you will ever need
% for Keck analysis.
%
% There are several flags that can be passed as optional keywords (listed
% below). It is important to note that some of these flags will completely
% change the return values of this function, so it returns [az, el, dk]
% with various transformations applied, instead of [r, theta, psi].
%
% Relevant postings that might be useful:
%   2012-01-11  rx0 Sidelobe Maps with Parallax Correction
%   2012-05-10  Pointing model for Keck beam maps
%   2012-11-20  Comparison of BICEP2 and Keck beam map pointing models
%   2012-12-07  Achieving agreement with the pipeline parallactic angle
%               convention
%   2014-02-03  Spherical coordinates for beam maps
%
% [Arguments]
%   az          Azimuth encoder coordinate. Can be a vector (timestream).
%   el          Elevation encoder coordinate. Can be a vector (timestream).
%   dk          Dk angle encoder coordinate. Can be a vector (timestream).
%   mount       Mount model data structure (see below for details).
%   mirror      Mirror model data structure (see below for details).
%   source      Source position data structure (see below for details).
%   p           Array data structure, from get_array_info.
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
%   Fields for mirror data structure:
%     mirror.height       Distance along the drum central axis from the
%                         height of the receiver aperture to the mirror, in
%                         meters. Default = 3.7186 meters.
%     mirror.tilt         Rotation of the mirror about an axis parallel to
%                         the telescope elevation axis, in degrees.
%                         Default = 40 degrees.
%     mirror.roll         Rotation of the mirror about an axis perpendicular
%                         to the tilt axis. Default = 0 degrees.
%
%   Fields for source data structure:
%     source.distance  Distance (in meters) to the source.
%     source.azimuth   Azimuth (in degrees) of source as seen from Keck.
%     source.height    Height (in meters) of source relative to the base
%                      of the Keck ground screen.
%
% [Optional keywords]
%   'NoMirror'  If this optional keyword is included, then no
%               mirror reflection will be applied (and the mirror
%               argument is ignored).
%   'NoSource'  Ignore the source argument and calculate (az, el, pa)
%               pointings for all detectors, without any parallax.
%               VERY IMPORTANT NOTE -- Using this keyword will completely
%               change the return values of this function. Instead of
%               [r, theta, psi], you will get [az, el, pa] where the
%               parallactic angle of the detector co-polar axis is measured
%               from the great circle connecting that detector pointing to
%               zenith.
%   'Legacy'    Use kbmp_parallax.m instead of kbmp_parallax_sph.m
%               VERY IMPORTANT NOTE -- Using this keyword will completely
%               change the return values of this function. Instead of
%               [r, theta, psi], you will get [az, el, pa] where the
%               parallactic angle of the detector co-polar axis is measured
%               from the great circle connecting that detector pointing to
%               zenith.
%   'BeamCentered'       The 'BeamCentered' and 'BoresightCentered' keywords
%   'BoresightCentered'  apply only if you also select 'Legacy'. If
%                        'BeamCentered' is selected, then the returned
%                        coordinates become [phi, theta, pa], which is a
%                        bit like the default [r, theta, psi] behavior,
%                        except in a different order and calculated
%                        differently.
%   'Buddy'     Get pointing of beam's buddy location instead of the beam
%               itself.  Does this by sending q.theta -> q.theta+180.
%
% [Returns]
%   All three outputs have shape [N,M], where N is the number of samples
%   for the (az, el, dk) input arguments and M is the number of detectors.
%   For details of how each coordinate is defined, read the 2014-02-03
%   posting on spherical coordinates.
%
%   r      Angular separation in degrees between the detector
%          pointing and the source.
%   theta  Direction of the source from the detector pointing,
%          relative to detector orientation axes.
%   psi    Parallactic angle of the detector co-polar axis measured at
%          the location of the source.

% Last update: 2014-02-08 CAB

% Parse optional arguments.
use_mirror = 1;
legacy_mode = 0;
find_buddy = 0;
no_chi=0;
for i=1:length(varargin)
    % 'NoMirror' keyword.
    if strcmp(varargin{i}, 'NoMirror')
        use_mirror = 0;
    end
    
    % 'Legacy' keyword.
    if strcmp(varargin{i}, 'Legacy')
        legacy_mode = 1;
    end
    
    % 'Buddy' keyword.
    if strcmp(varargin{i}, 'Buddy')
        find_buddy = 1;
    end
    
    % Check for 'NoChi' keyword.
    if strcmp(varargin{i}, 'NoChi')
        no_chi = 1;
    end
    % Other keywords are not needed in this function and are simply
    % passed on to kbmp_parallax_sph (or kbmp_parallax).
end

% Fill out mount data structure with reasonable values.
% Took these from January 11, 2012 posting and then converted from
% inches to meters.
if ~exist('mount', 'var')
    mount = [];
end
if ~isfield(mount, 'aperture_offr')
    mount.aperture_offr = 0.5458; % meters
end
if ~isfield(mount, 'aperture_offz')
    mount.aperture_offz = 1.5964; % meters
end
if ~isfield(mount, 'dk_offx')
    mount.dk_offx = -1.0196; % meters
end
if ~isfield(mount, 'dk_offy')
    mount.dk_offy = 0; % meters
end
if ~isfield(mount, 'el_tilt')
    mount.el_tilt = 0; % degrees
end
if ~isfield(mount, 'el_offx')
    mount.el_offx = 0; % degrees
end
if ~isfield(mount, 'el_offz')
    mount.el_offz = 1.1750; % meters
end
if ~isfield(mount, 'az_tilt_ha')
    mount.az_tilt_ha = 0; % degrees
end
if ~isfield(mount, 'az_tilt_lat')
    mount.az_tilt_lat = 0; % degrees
end

% Fill out mirror data structure with reasonable values.
if ~exist('mirror', 'var')
    mirror = [];
end
if ~isfield(mirror, 'height')
    % Email from Chris, 2012 May 04.
    mirror.height = 3.7186; % meters
end
if ~isfield(mirror, 'tilt')
    % Email from Chris, 2012 April 11.
    mirror.tilt = 40; % degrees
end
if ~isfield(mirror, 'roll')
    mirror.roll = 0; % degrees
end

% Fill out source data structure with reasonable values.
if ~exist('source', 'var')
    source = [];
end
if ~isfield(source, 'distance')
    % 211 meters from SPUD to DSL mast, from John's email: 24 May 2012.
    source.distance = 211; % meters
end
if ~isfield(source, 'azimuth')
    source.azimuth = 1.71;
end
if ~isfield(source, 'height')
    source.height = 13.2011;
end

% Allocate output arrays.
r = zeros(length(az), length(p.r));
theta = zeros(size(r));
psi = zeros(size(r));

% Aperture position depends on the drum angle, specified for each
% detector in p.drumangle. Rather than keeping track of multiple
% aperture positions, just loop over the unique drum angles found
% in p (i.e. the five Keck receivers).
drumangle = unique(p.drumangle);
for i=1:length(drumangle)
    % New data structure containing relevant information for these
    % pixels only.
    q.r = p.r(p.drumangle == drumangle(i));
    q.theta = p.theta(p.drumangle == drumangle(i));
    
    if find_buddy
        q.theta = q.theta + 180;
        if q.theta > 360
            q.theta = q.theta-360;
        end
    end
    
    q.chi = p.chi(p.drumangle == drumangle(i));
    q.chi_thetaref = p.chi_thetaref(p.drumangle == drumangle(i));
    q.drumangle = p.drumangle(p.drumangle == drumangle(i));
    
    % Calculate aperture position, boresight pointing, and boresight
    % orientation.
    [pos pnt ort] = kbmp_mount(az, el, dk, mount, drumangle(i));
    
    % Mirror reflection steps are optional.
    if use_mirror
        % Calculate mirror position and orientation.
        [mpos mnorm] = kbmp_mirror(az, el, mount, mirror);
        
        % Calculate reflected aperture position, boresight pointing, and
        % boresight orientation.
        [pos pnt ort] = kbmp_reflect(pos, pnt, ort, mpos, mnorm);
    end
    
    % Calculate parallax corrected pointing for all pixels.
    if legacy_mode
        [rq thetaq psiq] = kbmp_parallax(pos, pnt, ort, source, ...
            q, varargin{:});
    else
        [rq thetaq psiq] = kbmp_parallax_sph(pos, pnt, ort, source, ...
            q, varargin{:});
    end
    
    if no_chi
        r = rq;
        theta = thetaq;
        psi = psiq;
    else
        % Record output az/el/pa for these detectors.
        r(:, p.drumangle == drumangle(i)) = rq;
        theta(:, p.drumangle == drumangle(i)) = thetaq;
        psi(:, p.drumangle == drumangle(i)) = psiq;
    end
end
