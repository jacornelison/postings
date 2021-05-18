function kbmp_check()

%%
addpath('Z:/pipeline/util/')
load('../data/b3_array_info.mat')
%%
AZ = 49;%randi([0,359],1,1);
EL = 13;%randi([0,89],1,1);
DK = 46;%randi([0,179],1,1);

% Pointing Model
pm = [];
pm.flex_cos = 0;
pm.flex_sin = 0;
pm.az_tilt_ha = 0;
pm.az_tilt_lat = 0;
pm.el_tilt = 0;
pm.collim_x = 0;
pm.collim_y = 0;
pm.collim_mag = 0;
pm.collim_dir = 0;
pm.az_zero = 0;
pm.el_zero = 0;

% Mount struct
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

% Mirror struct
mirror = [];
mirror.height = 2;
mirror.tilt = 45;
mirror.roll = 0;

% Source Struct
source = [];
source.distance = 200;
source.azimuth = 180+AZ;
source.el = 90-EL;
source.height = source.distance*tand(source.el);

[rp,thetap,psi] = keck_beam_map_pointing(AZ,EL,DK,mount,mirror,source,p);

xp = 2 * sind(rp / 2) .* cosd(thetap) * 180 / pi;
yp = 2 * sind(rp / 2) .* sind(thetap) * 180 / pi;


[xpb,ypb,phidb] = beam_map_pointing_model(AZ,EL,DK,pm,mount,mirror,source,p);

figure(1)
clf; hold on;

plot(xp,yp,'bx',xpb,ypb,'ro')
xlabel('x'' (^o)')
ylabel('y'' (^o)')
grid on
title({'KBMP and BMPM per-detector offsets',sprintf('A=%i,E=%i,K=%i',AZ,EL,DK)})
legend({'Via KBMP','Via BMPM'});

figure(2)
clf; hold on;
plot(xp-xpb,yp-ypb,'.')
xlabel('\Deltax'' (^o)')
ylabel('\Deltay'' (^o)')
grid on
title('KBMP vs. BMPM Residuals')
