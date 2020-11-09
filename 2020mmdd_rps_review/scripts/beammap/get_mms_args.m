function [mount mirror source] = get_mms_args(run)
%
% function [mount mirror source] = get_mms_args(run)
%
% TSG 2018-07-26
% moved subfunction from beam map pipeline to standalone script
%
% This is a lookup function for mount/mirror/source values,
% given some FFBM run.  The are used as inputs to 
% keck_beam_map_pointing.m.  
%
% INPUTS:
%   run:      Keck runs:   'k0r4', 'k0r4t2' (tile 2 only), highbay
%                          'k083' = k7 run 2013-10, highbay
%                          'k201?' = Keck ffbm @ Pole (12,13,14,15)
%                          'k201?fsl' = Keck sidelobe @ Pole (13,14)
%             BICEP3 runs: 'hb3r5' = B3 run 5, highbay
%                          'b3r6' = B3 run 6 @ Pole, 2015-02
%                          'b3r8' = B3 run 8 @ Pole, 2016-02
%                          'b3r9' = B3 run 9 @ Pole, 2017/18-02
%
% OUTPUTS:
%   mount/mirror/source = output structs containing relevant params.
%   Any empty values means the default values in
%   keck_beam_map_pointing are appropriate for this run.
%  

switch run
  case {'k2012','k2013','k2013fsl','k2014','k2014fsl','k2015','k2016',...
        'k2017','k2018'}
    % For Keck, use mount/mirror defaults from keck_beam_map_pointing,
    mount = [];
    mirror = [];
    source.distance = 211; % meters, unchanged from keck_beam_map_pointing
  case 'hb3r5'
    mount.aperture_offr = 0;
    mount.aperture_offz = 0.997; % Aperture to el axis: 39.255"
    mount.dk_offx = 0;
    mount.dk_offy = 0;
    mount.drum_angle = 0;
    mount.el_tilt = 0;
    mount.el_offz = 0;
    mount.el_offx = 0;
    mount.az_tilt_ha = 0;
    mount.az_tilt_lat = 0;
    mirror.height = 0.954;       % Aperture to mirror: 37.561"
    mirror.tilt = 45;     
    mirror.roll = 0;
    source.distance = 70;
  case {'k0r4','k0r4t2','h083'}
    mount.aperture_offr = 0;
    mount.aperture_offz = 0.6;
    mount.dk_offx = 0;
    mount.dk_offy = 0;
    mount.drum_angle = 0;
    mount.el_tilt = 0;
    mount.el_offz = 0;
    mount.el_offx = 0;
    mount.az_tilt_ha = 0;
    mount.az_tilt_lat = 0;
    mirror.height = 0.7;
    mirror.tilt = 45;     %% Not sure about this one...
    mirror.roll = 0;
    source.distance = 70;
  case {'b3r6','b3r8','b3r9'}
    mirror.tilt = 45;
    mirror.roll = 0;
    mirror.height = 0.954;       % Aperture to mirror: 37.561"
    mount.aperture_offr = 0;
    mount.aperture_offz = 0.997; % Aperture to el axis: 39.255"
    mount.dk_offx = 0;
    mount.dk_offy = 0;
    mount.drum_angle = 0;
    mount.el_tilt = 0;
    mount.el_offz = 0;
    mount.el_offx = 0;
    mount.az_tilt_ha = 0;
    mount.az_tilt_lat = 0;
    source.distance = 190; % 210;
end

return