addpath('Z:\pipeline\util\')

load('../data/b3rpsfiles_2018.mat')


%%
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

AZ = 49;%randi([0,359],1,1);
EL = 13;%randi([0,89],1,1);
DK = 46;%randi([0,179],1,1);

d = [];
d.antenna0.tracker.encoder_mul = ones(length(AZ),3)*360;
d.antenna0.pmac.fast_az_pos = AZ;
d.antenna0.pmac.fast_el_pos = EL;
d.antenna0.pmac.fast_dk_pos = DK;
d.antenna0.tracker.encoder_off = repmat([0 0 0]*3.6e6,length(AZ),1);

d=invpointing_model(d,pm);
az_ideal = (d.pointing.hor.az);
el_ideal = d.pointing.hor.el;
dk_ideal = d.pointing.hor.dk;


[Ed,Ad] = reckon(el_ideal,-az_ideal,p.r,p.theta-dk_ideal-90);
Ad = -Ad;

[Adb,Edb,PAdb] = beam_map_pointing_model(AZ,EL,DK,pm,[],[],[],p);

figure(1)
clf; hold on;

plot(Ad,Ed,'bx',Adb,Edb,'ro')
xlabel('Az (^o)')
ylabel('El (^o)')
grid on
title({'Reckon and BMPM per-detector offsets',sprintf('A=%i,E=%i,K=%i',AZ,EL,DK)})
legend({'Via Reckon','Via BMPM'});

figure(2)
clf; hold on;
plot(Ad-Adb',Ed-Edb','.')
xlabel('\DeltaAz (^o)')
ylabel('\DeltaEl (^o)')
grid on
title('Reckon vs. BMPM Residuals')

