function [r, theta, psi] = kbmp_parallax_sph(pos, pnt, ort, source, p, ...
    varargin)
% [r, theta, psi] = kbmp_parallax_sph(pos, pnt, ort, source, p)
%
% Calculates parallax corrected source position in (r,theta) coordinates
% relative to all detectors.
%
% In practice, this function should not be used directly. Use the wrapper
% function, keck_beam_map_pointing.m
%
% For details of the calculation, see posting on Keck analysis logbook
% from 2014-Feb-03: "Spherical coordinates for beam maps".
%
% [Arguments]
%   The first three input arguments should be arrays with shape
%   [N,3], where the first dimension is the number of samples and
%   the second dimension are (x,y,z) horizon coordinates. These
%   should be obtained from kbmp_reflect.m for beam mapping with
%   the far field flat (or directly from kbmp_mount.m if no far
%   field flat).
%
%   pos     Position vector for the telescope aperture.
%   pnt     Boresight pointing vector.
%   ort     Boresight orientation vector.
%   source  Data structure containing source position parameters.
%   p       Array data structure, from get_array_info.
%
%   Fields for source data structure:
%     source.distance  Distance (in meters) to the source.
%     source.azimuth   Azimuth (in degrees) of source as seen from Keck.
%     source.height    Height (in meters) of source relative to the base
%                      of the Keck ground screen.
%
% [Optional keywords]
%   'NoMirror'  By default, this function assumes that the input boresight
%               and orientation vectors have been reflected off the far
%               field flat, so the pixel positions and orientations should
%               be constructed with a mirror reflection. Include this
%               optional keyword to get pixel positions and orientations
%               with normal (non-reflected) parity.
%   'NoSource'  Ignore the source argument and calculate (az,el,pa)
%               pointings for all detectors, without any parallax.
%               VERY IMPORTANT NOTE -- Using this keyword will completely
%               change the return values of this function. Instead of
%               [r', theta', psi], you will get [az, el, pa] where the
%               parallactic angle of the detector co-polar axis is measured
%               from the great circle connecting that detector pointing to
%               zenith.
%   'NoChi'     Ignore detector information and calculate source location
%               and polarization in terms of boresight coords [x, y, phi].
%
% [Returns]
%   All three outputs have shape [N,M], where N is the number of
%   samples and M is the number of pixels. For details of how each
%   coordinate is defined, read the posting mentioned above.
%
%   r'      Angular separation in degrees between the detector
%          pointing and the source.
%   theta'  Direction of the source from the detector pointing,
%          relative to detector orientation axes.
%   psi    Parallactic angle of the detector co-polar axis measured at
%          the location of the source.

% Last updated: 2014-02-08 CAB

% Parse optional arguments.
no_mirror = 0;
no_source = 0;
no_chi=0;
for i=1:length(varargin)
    % Check for 'NoMirror' keyword.
    if strcmp(varargin{i}, 'NoMirror')
        no_mirror = 1;
    end
    
    % Check for 'NoSource' keyword.
    if strcmp(varargin{i}, 'NoSource')
        no_source = 1;
    end
    
    % Check for 'NoChi' keyword.
    if strcmp(varargin{i}, 'NoChi')
        no_chi = 1;
    end
end

% Check size of 3-vector input arguments.
nsamp = size(pos, 1);
if (size(pos, 2) ~= 3)
    disp('[kbmp_parallax_sph] ERROR: argument pos must have size [N,3]');
    return
end
if (size(pnt, 1) ~= nsamp) | (size(pnt, 2) ~= 3)
    disp('[kbmp_parallax_sph] ERROR: argument pnt must have size [N,3]');
    return
end
if (size(ort, 1) ~= nsamp) | (size(ort, 2) ~= 3)
    disp('[kbmp_parallax_sph] ERROR: argument ort must have size [N,3]');
    return
end

% Calculate the cross-product of boresight pointing and
% orientation. This will serve as the third axis for an orthonormal
% coordinate system centered around the boresight.
% Flip the sign of this vector depending on desired parity.
if no_mirror
    % Right-handed coordinates:
    %   pnt = zhat, ort = xhat, yhat = xhat CROSS zhat
    ort2 = cross(pnt, ort, 2);
else
    % Left-handed coordinates:
    %   pnt = zhat, ort = xhat, yhat = xhat CROSS zhat
    ort2 = cross(ort, pnt, 2);
end

% Calculate a unit vector for the direction from the aperture to
% the three-dimensional source position.
source_pos = [source.distance .* cosd(source.azimuth), ...
    -1 * source.distance .* sind(source.azimuth), ...
    source.height];
if size(source_pos, 1) == 1
    source_aperture = repmat(source_pos, nsamp, 1) - pos;
elseif size(source_pos, 1) == nsamp
    source_aperture = source_pos - pos;
else
    disp(['[kbmp_parallax_sph] ERROR: Source position should either ' ...
        'be fixed or else specified for every sample.']);
    return;
end
norm = sqrt(sum(source_aperture.^2, 2));
source_aperture = source_aperture ./ repmat(norm, 1, 3);

if no_chi
    % Calculate r and theta coordinates of the source.
    % NOTE: Need to take the real part of r for the case where
    % source_aperture DOT pnt is very close to +-1 and floating point
    % errors push us past the edge.
    r = real(acos(sum(source_aperture .* pnt, 2))) * 180 / pi;
    theta = atan2(-1 * sum(source_aperture .* ort2, 2), ...
        sum(source_aperture .* ort, 2)) * 180 / pi;
    
    % Rotate from detector pointing to location of the source.
    [sort, sort2, spnt] = rtheta_rotation(r(:,i), theta(:,i), ...
        ort, ort2, pnt);
    
    % Construct source co-polar orientation. Assumes source is
    % pointed azimuthally at the telescope, but not in elevation.
    % Also assumes plane of rotation is made by z-axis and axis
    % tangent to z and source pointing. If no rotation angle is
    % specified, asssume source pol is pointed at local zenith.
    if isfield(source,'abs_rot')
        zen_vec = repmat([0,0,1],nsamp,1);
        e_vec = cross(zen_vec,source_aperture);
        norm = sqrt(sum(e_vec.^2,2));
        e_vec = e_vec./repmat(norm,1,3);
        
        rot = source.abs_rot;
        
        ep = cosd(rot)*zen_vec-sind(rot)*e_vec;
    else
        disp('[kbmp_parallax_sph] ''NoChi'' selected, but no source rotation provided. Assuming source at zenith.')
        ep = repmat([0,0,1],nsamp,1);
    end

    % return phi as the output.
    psi = atan2(-1*dot(ep,sort2,2),dot(ep,sort,2))*180/pi;
    
    % Lastly, Calculate x and y coordinates of the source
    r = 2 * sind(r / 2) .* cosd(theta) * 180 / pi;
    theta = 2 * sind(r / 2) .* sind(theta) * 180 / pi;
else
    % Allocate output arrays.
    npix = length(p.r);
    r = zeros(nsamp, npix);
    theta = zeros(nsamp, npix);
    psi = zeros(nsamp, npix);
        
    % Loop over channels.
    for i=1:npix
        % Rotate boresight basis vectors to the pointing of this channel.
        [port, port2, ppnt] = rtheta_rotation(p.r(i), p.theta(i), ...
            ort, ort2, pnt);
        
        % If NoSource is set, just stop here and return az/el/dk in place
        % of r/theta/psi.
        if no_source
            % az instead of r.
            r(:,i) = atan2(-1 * ppnt(:,2), ppnt(:,1)) * 180 / pi;
            
            % el instead of theta.
            theta(:,i) = asin(ppnt(:,3)) * 180 / pi;
            
            % Instead of psi, calculate the parallactic angle of chi with
            % respect to zenith.
            % First, calculate copolar orientation vector.
            cc = cosd(p.chi_thetaref(i) + p.chi(i));
            sc = sind(p.chi_thetaref(i) + p.chi(i));
            copol = cc * port - sc * port2;
            % Basis vectors for parallactic angle calculation.
            eb = cross(repmat([0, 0, 1], nsamp, 1), ppnt, 2);
            norm = sqrt(sum(eb.^2, 2));
            eb = eb ./ repmat(norm, 1, 3);
            ea = cross(ppnt, eb, 2);
            % Calculate parallactic angle.
            psi(:,i) = atan2(sum(copol .* eb, 2), ...
                sum(copol .* ea, 2)) * 180 / pi;
        else
            
            % Calculate a unit vector for the direction from the aperture to
            % the three-dimensional source position.
            source_pos = [source.distance .* cosd(source.azimuth), ...
                -1 * source.distance .* sind(source.azimuth), ...
                source.height];
            if size(source_pos, 1) == 1
                source_aperture = repmat(source_pos, nsamp, 1) - pos;
            elseif size(source_pos, 1) == nsamp
                source_aperture = source_pos - pos;
            else
                disp(['[kbmp_parallax_sph] ERROR: Source position should either ' ...
                    'be fixed or else specified for every sample.']);
                return;
            end
            norm = sqrt(sum(source_aperture.^2, 2));
            source_aperture = source_aperture ./ repmat(norm, 1, 3);
            
            % Calculate r and theta coordinates of the source.
            % NOTE: Need to take the real part of r for the case where
            % source_aperture DOT ppnt is very close to +-1 and floating point
            % errors push us past the edge.
            r(:,i) = real(acos(sum(source_aperture .* ppnt, 2))) * 180 / pi;
            theta(:,i) = atan2(-1 * sum(source_aperture .* port2, 2), ...
                sum(source_aperture .* port, 2)) * 180 / pi;
            
            % Rotate from detector pointing to location of the source.
            [sort, sort2, spnt] = rtheta_rotation(r(:,i), theta(:,i), ...
                port, port2, ppnt);
            
            % Rotate orientation vectors by chi_thetaref + chi to obtain
            % detector copolar direction.
            % NOTE: Convention for chi (and chi_thetaref) is such that this
            % rotation is negative in my basis.
            sc = sind(p.chi_thetaref(i) + p.chi(i));
            cc = cosd(p.chi_thetaref(i) + p.chi(i));
            copol = cc * sort - sc * sort2;
            
            % Calculate angle between source polarization and copolar
            % direction.
            % NOTE: For now, source polarization is always assumed to be
            % the direction of increasing elevation.
            eb = cross(repmat([0, 0, 1], nsamp, 1), spnt, 2);
            norm = sqrt(sum(eb.^2, 2));
            eb = eb ./ repmat(norm, 1, 3);
            ea = cross(spnt, eb, 2);
            psi(:,i) = atan2(sum(copol .* eb, 2), ...
                sum(copol .* ea, 2)) * 180 / pi;
            
        end
    end
end


% Subfunction: rtheta_rotation
% ---
% This function takes a set of input vectors (e1, e2, e3), which
% are intended to be an orthonormal basis, and then applies the
% following solid rotation:
%   1. Rotate by theta about e3, to obtain e1', e2', e3'
%   2. Rotate by -r about e2', to obtain e1'', e2'', e3''
%   3. Rotate by -theta about e3'', to obtain e1out, e2out, e3out
function [e1out, e2out, e3out] = rtheta_rotation(r, theta, e1, e2, e3)

% Check dimensions of inputs.
%   * r and theta should either be scalar or length N vectors.
%   * e1, e2, and e3 should either have size [1,3] or [N,3].
n = max([length(r), length(theta), size(e1, 1), size(e2, 1), size(e3, 1)]);
if (size(e1, 2) ~= 3) || (size(e2, 2) ~= 3) || (size(e3, 2) ~= 3)
    disp('[rtheta_rotation] ERROR: basis vectors must be 3-vectors');
    return
end
% Expand length = 1 dimensions to N.
if length(r) == 1
    r = repmat(r, n, 1);
end
if length(theta) == 1
    theta = repmat(theta, n, 1);
end
if size(e1, 1) == 1
    e1 = repmat(e1, n, 1);
end
if size(e2, 1) == 1
    e2 = repmat(e2, n, 1);
end
if size(e3, 1) == 1
    e3 = repmat(e3, n, 1);
end

% Precompute factors.
sr = sind(r);
cr = cosd(r);
st = sind(theta);
ct = cosd(theta);

% Calculate rotated basis.
e1out = repmat(cr .* ct.^2 + st.^2, 1, 3) .* e1 + ...
    repmat((1 - cr) .* ct .* st, 1, 3) .* e2 - ...
    repmat(sr .* ct, 1, 3) .* e3;
e2out = repmat((1 - cr) .* ct .* st, 1, 3) .* e1 + ...
    repmat(cr .* st.^2 + ct.^2, 1, 3) .* e2 + ...
    repmat(sr .* st, 1, 3) .* e3;
e3out = repmat(sr .* ct, 1, 3) .* e1 - ...
    repmat(sr .* st, 1, 3) .* e2 + ...
    repmat(cr, 1, 3) .* e3;
