function moon_driver_driver_2022()
savedir = 'rps_data/2022_moon_rerun/';

load(fullfile(savedir,'moonsch.mat'))

% From Moon Analysis
mirror = struct();
mirror.height = 1.4592;
mirror.tilt = 44.887;
mirror.roll = -0.063;

moonopt = struct();
moonopt.min_dist = 0.1;
moonopt.maskind = 600;
moonopt.model = get_pointing_model(59580);
moonopt.mirror = mirror;

schind = 3;
scanind = 10;
%chind = 36;
%moontod = moon_reduc(moonsch,schind,scanind,'savedir',savedir,'moonopt',moonopt);
moonparm = moon_fit_centers(moonsch,schind,scanind,'savedir',savedir);
%moon_reduc_driver(moonsch,moonopt,savedir)
moon_farm_center_fits(moonsch,savedir)
