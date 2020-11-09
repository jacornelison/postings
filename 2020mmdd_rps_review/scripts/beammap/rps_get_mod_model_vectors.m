function [model debug] = rps_get_mod_model_vectors(aparam,eta,A,B3,B1,source,MIRROR,PLOT)
%function [model debug] = rps_get_mod_model_vectors(aparam,eta,A,B3,B1,source,MIRROR,PLOT)
% Fit RPS data using vector-based model. Intended to more accurately
% represent pointing model.
% 
% To-do: explain what all the inputs are for.

if any(isnan(aparam))
    disp('NaN''s detected in parameters, returning extrememly high values.')
    model = 1e10*ones(1,length(eta));
    return    
end

if ~exist('PLOT','var')
    PLOT = false;
end

if ~exist('MIRROR','var')
    MIRROR = true;
end

if ~isfield(source,'type')
    source.type = 'gauss';
end
if ~isfield(source,'beamwidth')
    source.beamwidth = 7.69;
end

phi = aparam(1);
eps = aparam(2);
e_collim = aparam(3); % E-plane nutation
h_collim = aparam(4); % H-plane nutation
az = aparam(5); % Az misalignment
el = aparam(6); % El misalignment
gain = aparam(7);

if MIRROR
    B2 = cross(B1,B3);
else
    B2 = cross(B3,B1);
end

% Get source position
S = [source.distance .* cosd(source.azimuth), ...
                -1 * source.distance .* sind(source.azimuth), ...
                source.height];

% Create Source pointing and orientation vectors
S3 = (S-A)/norm(S-A);

% find r/theta coordinates
r = acosd(dot(S3,B3));
theta = atan2(-dot(S3,B2),dot(S3,B1))*180/pi;

% Parallel transport boresight orientation vectors
[S1,S2,S3] = rtheta_rotation(r,theta,B1,B2,B3);

% Create Source orientation vectors:
Sphi = S1*cosd(phi)-S2*sind(phi);
Sphi2 = cross(Sphi,S3);

% Create local source axes. 2 orientation options:
% Source starts aligned locally to zenith pointing towards telescope, but
% level wrt to gravity.
% 2nd option: source is pointing directly at telescope and aligned toward
% zenith in plane tangent to pointing.
% Option 2 by default.
% xp - Source copolar axis
% zp - Source rotation axis
% yp - Source cross-polar axis.

if 0
    xp = [0,0,1]; 
    yp = cross(S,xp)/norm(cross(S,xp));
    zp = cross(xp,yp);
else
    zp = S/norm(S);
    yp = cross(zp,[0,0,1]);
    xp = cross(yp,zp);
end

% Rotate axes by nutation
r = sqrt(e_collim^2+h_collim^2);
theta = atan2(h_collim,e_collim)*180/pi;

Vmat0 = [xp' yp' zp'];
R1 = rodmatrix(theta,zp);
xyzdot = R1*Vmat0;
R2 = rodmatrix(-r,xyzdot(:,2));
xyzddot = R2*xyzdot;
R3 = rodmatrix(-theta,xyzddot(:,3));
xyzpp = R3*xyzddot;

% then misalignment
r = sqrt(az^2+el^2);
theta = atan2(el,az)*180/pi;


% Build the rotation matrices
R1 = rodmatrix(theta,zp);
xyzdot = R1*Vmat0;
R2 = rodmatrix(-r,xyzdot(:,2));
xyzddot = R2*xyzdot;
R3 = rodmatrix(-theta,xyzddot(:,3));

if PLOT
    Vmatp=R3*xyzddot;
    figure(1)
    hold off
    quiver3(0,0,0,zp(1),zp(2),zp(3),'color',[0,0,0])
    hold on
end

% Calculate the amplitude for each eta
[copol,crosspol,amp] = deal([]);
N = length(eta);
Spdb = NaN(N,3);
for i = 1:N
    R0 = rodmatrix(eta(i),zp);
    Srot = R3*R2*R1*R0*xyzpp;
    Spp = Srot(:,1);
    Stp = Srot(:,3);
    Sp2 = cross(Spp,S3)/norm(cross(Spp,S3));
    Sp = cross(S3,Sp2);
    amp(end+1) = gain.*gaussmf(acosd(dot(S3,Stp)),[source.beamwidth,0]);
    copol(end+1) = dot(Sphi,Sp).^2;
    crosspol(end+1) = dot(Sphi2,Sp).^2;
    
    if PLOT
        quiver3(0,0,0,Sp(1),Sp(2),Sp(3),'color',[0,0,0.75])
        quiver3(0,0,0,Stp(1),Stp(2),Stp(3),'color',[0.75,0,0])
        quiver3(0,0,0,Vmatp(1,3),Vmatp(2,3),Vmatp(3,3),'color',[0,0.75,0])
    end
    Spdb(i,:) = Sp;
end

if strcmp(source.type,'gauss')
    model = amp.*(copol+eps*crosspol);
else
    model = gain*(copol+eps*crosspol);
end

debug.S1 = S1;
debug.S2 = S2;
debug.S3 = S3;
debug.coax = xp;
debug.crossax = yp;
debug.rotax = zp;
debug.Sp = Spdb;
debug.Spp = Spp;
debug.Stp = Stp;
% debug.coaxp = coaxp;
% debug.crossaxp = crossaxp;
% debug.rotaxp = rotaxp;
debug.Dphi = Sphi;
% debug.e1 = e1;
% debug.e2 = e2;
debug.copol = copol;
debug.crosspol = crosspol;
debug.amp = amp;

if PLOT
   xlim([-2 2]);ylim([-2 2]);zlim([-2 2]);
   xlabel('x');ylabel('y');zlabel('z');
end


function v = rodrigues(v0,vaxis,angle)
v = [];
for i = 1:length(angle)
v(end+1,:) = v0.*cosd(angle(i)) + cross(vaxis,v0).*sind(angle(i)) + vaxis.*dot(vaxis,v0).*(1-cosd(angle(i)));
end

function R = rodmatrix(angle,vaxis)

crossk = [0 -vaxis(3) vaxis(2);...
    vaxis(3) 0 -vaxis(1);...
    -vaxis(2) vaxis(1) 0];

R = eye(3) + sind(angle)*crossk + (1-cosd(angle))*(crossk*crossk);

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

