function [phiQ,aparam,bparam,rp,thp] = rps_error_prop(r,theta,dk,varargin)



fields = {'az_zero','el_zero','mntaperture_offz','mirtilt','mirroll',...
    'mirheight','srcheight','srcaz','srcdist','gridhor','srcnut1','srcnut2',...
    'srcalign1','srcalign2'};

err = [];
for i = 1:length(fields)
    err.(fields{i}) = 0;
    for j = 1:length(varargin)
        if strcmp(fields{i},varargin{j})
            err.(fields{i}) = varargin{j+1};
        end
    end
end



%% Construct pointing data

% Set up az/el/dk as if it were coming from GCP to use as input into the
% inverse pointing model.

% Actual pointing
az0 = 90;
el0 = 90.1;
dk0 = -1*dk;


err.el_tilt = 0;
err.az_tilt_ha = 0;
err.az_tilt_lat = 0;

mount.aperture_offr =  0;
mount.aperture_offz = 1;
mount.dk_offy = 0;
mount.dk_offx = 0;
mount.el_tilt = 0;
mount.el_offx = 0;
mount.el_offz = 0;
mount.az_tilt_ha = 0;
mount.az_tilt_lat = 0;

mirror.tilt = 45.0;
mirror.roll = 0.0;
mirror.height = 1.0;

source.azimuth = -90.0;
source.distance = 200.0;
source.height = 0.0;


[A0,B30,B10] = kbmp_mount(az0,el0,dk0,mount,0);
            
% Reflect with mirror
[mpos mnorm] = kbmp_mirror(az0,el0,mount,mirror);
[A0,B30,B10] = kbmp_reflect(A0, B30, B10, mpos, mnorm);
B20 = cross(B10,B30);

% Go somewhere on focal plane
[B10,B20,B30] = rtheta_rotation(r,theta,B10,B20,B30);

eta = -180:30:180;

Ap = [00 0 0 0 0 0 1];
Bp = [90 0 0 0 0 0 1];

Amodel = rps_get_mod_model_vectors(Ap,eta,A0,B30,B10,source);
Bmodel = rps_get_mod_model_vectors(Bp,eta,A0,B30,B10,source);


% Where we think we're pointing
d.antenna0.tracker.encoder_mul = [360 360 360];
d.antenna0.tracker.encoder_sign = [1 1 1];
d.antenna0.tracker.encoder_off = [0 0 0];

d.antenna0.pmac.fast_az_pos = az0;
d.antenna0.pmac.fast_el_pos = el0;
d.antenna0.pmac.fast_dk_pos = -1*dk;

d = invpointing_model(d,err);

az = d.pointing.hor.az;
el = d.pointing.hor.el;
dk = d.pointing.hor.dk;


mount.aperture_offr =  0;
mount.aperture_offz = 1 + err.mntaperture_offz;
mount.dk_offy = 0;
mount.dk_offx = 0;
mount.el_tilt = 0;
mount.el_offx = 0;
mount.el_offz = 0;
mount.az_tilt_ha = 0;
mount.az_tilt_lat = 0;

mirror.tilt = 45.0 + err.mirtilt;
mirror.roll = 0.0 +err.mirroll;
mirror.height = 1.0 + err.mirheight;

source.azimuth = -90.0 + err.srcaz;
source.distance = 200.0 + err.srcdist;
source.height = 0.0 + err.srcheight;

[A,B3,B1] = kbmp_mount(az,el,dk,mount,0);
            
% Reflect with mirror
[mpos mnorm] = kbmp_mirror(az,el,mount,mirror);
[A,B3,B1] = kbmp_reflect(A, B3, B1, mpos, mnorm);
B2 = cross(B1,B3);
rp = acosd(dot(B3,B30));
thp = atan2(-dot(B3,B20),dot(B3,B10))*180/pi;

% Go somewhere on focal plane
[B1,B2,B3] = rtheta_rotation(r,theta,B1,B2,B3);

n1 = err.srcnut1;
n2 = err.srcnut2;
a1 = err.srcalign1;
a2 = err.srcalign2;

%% Construct modulation curves for both A and B dets
eta = (-180:30:180)+err.gridhor;

Aguess = [0 0 n1 n2 a1 a2 1];
Bguess = [90 0 n1 n2 a1 a2 1];

freepar.free = [1 1 0 0 0 0 1];
freepar.lb = [-10 -1 -15 -15 -15 -15 0];
freepar.ub = [10 1 15 15 15 15 1e6];

% Estimate parameters
[aparam, aerr, agof, astat, acov] = matmin('chisq',...
    Aguess, freepar,	'rps_get_mod_model_vectors',Amodel,ones(size(Amodel))*0.1,eta,A,B3,B1,source);

freepar.free = [1 1 0 0 0 0 1];
freepar.lb = [80 -1 -15 -15 -15 -15 0];
freepar.ub = [110 1 15 15 15 15 1e6];

Bguess = [90 0 0 0 0 0 1];

[bparam, berr, bgof, bstat, bcov] = matmin('chisq',...
    Aguess, freepar,	'rps_get_mod_model_vectors',Bmodel,ones(size(Amodel))*0.1,eta,A,B3,B1,source);

Q = (cosd(2*aparam(1))-cosd(2*bparam(1)))/2;
U = (sind(2*aparam(1))-sind(2*bparam(1)))/2;
phiQ = atand(U./Q)/2;

%keyboard()

