function [phipair,phia,phib,Ainc,epspair] = calc_mirror_diff_pol(p,n,k,dk)

if isempty(p)
p.r = 0;
p.theta = 0;
p.drumangle = 0;
p.chi = 0;
p.chi_thetaref=0;
end
% n = 3000;
% k = 3000;
mirror = [];
mirror.height = 2;
mirror.tilt = 45;
mirror.roll = 0;

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

dk = reshape(dk,[],1);

%% Calculate the ray going onto the mirror
az = dk*0;
el = dk*0;
ndk = length(az);
ndet = length(p.r);


% Calculate aperture position, boresight pointing, and boresight
% orientation.
[Apt, B3, B1] = kbmp_mount(az, el, dk, mount, 0);

% Right-handed coordinates:
%   pnt = B3, ort = B1, yhat = B3 cross B1
B2 = cross(B3, B1, 2);

% Calculate mirror position and orientation.
%[mpos, a, mtilt] = kbmp_mount(dk*0, dk*0, dk*0, mount, 0);

[mpos, mnorm] = kbmp_mirror(0, 0, mount, mirror);
%keyboard()
% We now have a plane defined by the mirror norm.
%[x_mirr, y_mirr, h] = deal(NaN(size(x)));
%[D1, D2,D3,Spol,Ppol] = deal(NaN(ndet,ndk,3));
[Ainc,phia,phib] = deal(NaN(ndet,ndk));
for dkind = 1:ndk
    for detind = 1:ndet
        % Convert pointing to cartesian
        [D1, D2, D3] = rtheta_rotation(p.r(detind), p.theta(detind), ...
            (B1(dkind,:)), (B2(dkind,:)), (B3(dkind,:)));
        
        % s-polarized vector is just mirror norm
        Spol = dot(D1,mnorm)*mnorm;
        Spolb = dot(-D2,mnorm)*mnorm;
        
        % p-polarized vector is det pointing cross mirror normal
        e1 = cross(D1,mnorm)/norm(cross(D1,mnorm));
        e2 = cross(mnorm,e1)/norm(cross(mnorm,e1));
        Ppol = dot(D1,e2)*e2;
        
        e3 = cross(-D2,mnorm)/norm(cross(-D2,mnorm));
        e4 = cross(mnorm,e3)/norm(cross(mnorm,e3));
        Ppolb = dot(-D2,e4)*e4;
        
        % Incident angle
        Ainc(detind,dkind) = abs(acosd(dot(mnorm,-D3)));
        
        % Reflection coefficients
        A = Ainc(detind,dkind);
        Atx = sqrt(1-Nair/Nal*sind(A).^2);
        Rs = abs((Nair*cosd(A)-Nal*cosd(Atx))./(Nair*cosd(A)+Nal*cosd(Atx)))^2;
        Rp = abs((Nal*cosd(A)-Nair*cosd(Atx))./(Nal*cosd(A)+Nair*cosd(Atx)))^2;
        
        % Multiply the reflection coefficients.
        % Ideal phi A is in the direction of B1.
        Va = Rs*Spol+Rp*Ppol;
        Va = Va/norm(Va);
        
        % Pol B is the reverse of Pol A
        Vb = Rs*Spolb+Rp*Ppolb;
        Vb = Vb/norm(Vb);
        
        % Calc the angles
        % Calc the angles
        phia(detind,dkind) = real(atan2(dot(Va,-D2),dot(Va,D1)))*180/pi;
        
        phib(detind,dkind) = real(atan2(dot(Vb,-D2),dot(Vb,D1)))*180/pi;
    end
end

pa = round(phia,2);
pb = round(phib,2);
Q = (cosd(2*pa)-cosd(2*pb))/2;
U = (sind(2*pa)-sind(2*pb))/2;
phipair = 1/2*atan2(U,Q)*180/pi;
epspair = 1-sqrt(U.^2+Q.^2);


function quiv3(V1,V2,varargin)
quiver3(V1(1,1),V1(1,2),V1(1,3),V2(1,1),V2(1,2),V2(1,3),varargin{:})
