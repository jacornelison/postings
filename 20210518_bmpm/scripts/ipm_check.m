function ipm_check()

addpath(fullfile('Z:','pipeline','util'))
addpath(fullfile('Z:','dev'))
%% Init parameters
% Boresight
bs.r = 0;
bs.theta = 0;
bs.chi = 0;
bs.chi_thetaref = 0;
bs.drumangle = 0;
bs.expt = 'bicep3';

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

% Mirror
mirror = [];
mirror.height = 2;
mirror.tilt = 45;
mirror.roll = 0;

% Pointing Model
pm = [];
pm.meas_name = {'170203 02:32:23'};
pm.mjd = 5.7787e+04;
pm.flex_cos = 0;
pm.flex_sin = 0;
pm.az_tilt_ha = 0;%*0.0046;
pm.az_tilt_lat = 0;%0.0044;
pm.el_tilt = 0;%*0.0286;
pm.collim_x = 0;
pm.collim_y = 0;
pm.collim_mag = 0;%-0.5802;
pm.collim_dir = 0;%143.9816;
pm.az_zero = 0;
pm.el_zero = 0;
pm.err_flex_cos = 0;
pm.err_flex_sin = 0;
pm.err_az_tilt_ha = 0.0614;
pm.err_az_tilt_lat = 0.0610;
pm.err_el_tilt = 0.2295;
pm.err_collim_x = 0;
pm.err_collim_y = 0;
pm.err_collim_mag = 0.0585;
pm.err_collim_dir = 5.8488;
pm.err_az_zero = 0.4323;
pm.err_el_zero = 0.0604;
pm.rms_res = 0.0050;
pm.numpts = 32;

source = [];
source.distance = 200;
source.azimuth = 180;
source.el = 2.5;
source.height = source.distance*tand(source.el);

%% Check our new code against inv_pointing model
% Pointing
colormap('default');

az0 = 0:360;%
el0 = 30.9:89.9;%(lims(1):diff(lims)/(res-1):lims(2))+90;%zeros(1,res)+85;
[AZ,EL] = meshgrid(az0,el0);
AZ = reshape(AZ,[],1);
EL = reshape(EL,[],1);
DK = zeros(size(AZ))+0;
T = 1:length(AZ);

d = [];
d.antenna0.tracker.encoder_mul = ones(length(AZ),3)*360;
d.antenna0.pmac.fast_az_pos = AZ;
d.antenna0.pmac.fast_el_pos = EL;
d.antenna0.pmac.fast_dk_pos = DK;
d.antenna0.tracker.encoder_off = repmat([0 0 0]*3.6e6,length(AZ),1);

for latind = 0.005%-1:1
    for haind = 0.005%-1:1
        for elind = 0.03%-1:1
    pm.az_tilt_ha = haind;%*0.0046;
    pm.az_tilt_lat = latind;%0.0044;
    pm.el_tilt = elind;%*0.0286;
    
d=invpointing_model(d,pm);
az_ideal = wrapTo360(d.pointing.hor.az);
el_ideal = d.pointing.hor.el;
dk_ideal = wrapTo360(d.pointing.hor.dk);

[az_command, el_command, dk_command] = encoder_to_command_coords(d);

[az, el, pa] = beam_map_pointing_model(az_command,el_command,dk_command,pm,mount,[],[],bs);
az = wrapTo360(az);
dk = wrapTo360(-1*(pa+90));

fig = figure(1);
%fig = figure('Visible','off');
set(fig,'position',[200,300,1500,500]);
clf;
subplot(1,3,1)
res = medfilt1(wrapTo360(az_ideal-az+90)-90);
imagesc('XData',az0,'YData',el0,'CData',reshape(res,length(az0),length(el0)));
xlim([0,360])
ylim([30,90])
c = colorbar();
c.Label.String = '\DeltaAz (degrees)';
xlabel('Azimuth (^o)')
ylabel('Elevation (^o)')
title('AZ Residuals')
grid on
set(gca,'GridAlpha',0.25,'Layer','top');

subplot(1,3,2)
res = medfilt1(wrapTo360(el_ideal-el+90)-90);
imagesc('XData',az0,'YData',el0,'CData',reshape(res,length(az0),length(el0)));
xlim([0,360])
ylim([30,90])
c = colorbar();
c.Label.String = '\DeltaEl (degrees)';
xlabel('Azimuth (^o)')
title('EL Residuals')
grid on
set(gca,'YTickLabel',{},'GridAlpha',0.25,'Layer','top');

subplot(1,3,3)
res = medfilt1(wrapTo360(dk_ideal-dk+90)-90);
imagesc('XData',az0,'YData',el0,'CData',reshape(res,length(az0),length(el0)));
xlim([0,360])
ylim([30,90])
c = colorbar();
c.Label.String = '\DeltaDK (degrees)';
xlabel('Azimuth (^o)')
title('DK Residuals')
grid on
set(gca,'YTickLabel',{},'GridAlpha',0.25,'Layer','top');
sgtitle({'Pointing model residuals: IMP minus BMPM',sprintf('Az Tilt Lat: %i^o Az Tilt HA: %i^o El Tilt HA: %i^o',latind,haind,elind)})

fname = sprintf('pm_compare_lat_%i_ha_%i_el_%i',latind,haind,elind);
dirname = fullfile('C:','Users','James','Documents','GitHub','postings','rps_scratch','figs','');
%saveas(fig,fullfile(dirname,fname),'png')
        end
    end
end

