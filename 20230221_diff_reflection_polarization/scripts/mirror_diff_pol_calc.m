function [phi_out] = mirror_diff_pol_calc(phi_in,r,theta,n,k,dk,mtilt,mroll)
% function [phi_corr] = mirror_diff_pol_correction(p,mirror,n,k)
% phi is pol angle
% r/theta is det pointing
% n is real part of index
% k is imaginary part of index
% dk is deck angle
% mirror is the mirror struct


mount = [];
mount.aperture_offr = 0;
mount.aperture_offz = 0;
mount.dk_offy = 0;
mount.dk_offx = 0;
mount.el_tilt = 0;
mount.el_offx = 0;
mount.el_offz = 0;
mount.az_tilt_ha = 0;
mount.az_tilt_lat = 0;
mount.az_offz = 0;
mount.az_tilt_org = 0; % 1.22l;
mount.el_tilt_org = 0;

Nair = 1;
Nal = n+k*1i;

phi_in = reshape(phi_in,[],1);
r = reshape(r,[],1);
theta = reshape(theta,[],1);
dk = reshape(dk,[],1);


ndet = length(r);
mirror = struct();
mirror.height = 2;

% We now have a plane defined by the mirror norm.
[phi_out] = deal(NaN(ndet,1));
for detind = 1:ndet
    % Calculate aperture position, boresight pointing, and boresight
    % orientation.
    [Apt, B3, B1] = kbmp_mount(0, 0, dk(detind), mount, 0);

    % Right-handed coordinates:
    %   pnt = B3, ort = B1, yhat = B3 cross B1
    B2 = cross(B3, B1, 2);

    % Calculate mirror position and orientation.
    %[mpos, a, mtilt] = kbmp_mount(dk*0, dk*0, dk*0, mount, 0);
    mirror.tilt = mtilt(detind);
    mirror.roll = mroll(detind);
    [mpos, mnorm] = kbmp_mirror(0, 0, mount, mirror);

    % Convert pointing to cartesian
    [D1, D2, D3] = rtheta_rotation(r(detind), theta(detind), ...
        B1, B2, B3);

    % Create a vector that represents the co-polar axis of the detector
    Dp = cosd(phi_in(detind))*D1-sind(phi_in(detind))*D2;

    % Define the normal to the plane of incidence:
    int_norm = cross(D3,mnorm)/norm(cross(D3,mnorm));
    
    
    % S polariation is the projection onto the vector normal to the
    % plane of incidence:
    Spolp = dot(Dp,int_norm)*int_norm;

    % P polarization is the projection onto the incidence plane
    % which is just the orientation vector minus the projection on the vector
    % normal to that plane (Spol[a,b,p]):
    Ppolp = Dp-Spolp;

    % Incident angle
    Ainc = (acosd(dot(mnorm,-D3)));

    % Reflection coefficient
    Atx = sqrt(1-(Nair/Nal*sind(Ainc)).^2);
    Rs = abs((Nair*cosd(Ainc)-Nal*Atx)./(Nair*cosd(Ainc)+Nal*Atx))^2;
    Rp = abs((Nair*Atx-Nal*cosd(Ainc))./(Nair*Atx+Nal*cosd(Ainc)))^2;

    % Multiply the reflection coefficients.
    Vp = Spolp*Rs+Ppolp*Rp;
    Vp = Vp/norm(Vp);

    % Calc the angles
    phi_out(detind) = real(atan2(dot(Vp,-D2),dot(Vp,D1)))*180/pi;
    
end
