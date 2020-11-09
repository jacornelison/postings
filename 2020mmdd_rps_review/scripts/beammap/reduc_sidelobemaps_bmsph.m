function reduc_sidelobemaps_bmsph(t1,t2,filename, rxNum, skip, run, x_bin)
%Accumulate demodulated sidelobe data into map in x_P, y_P coordinates
%t1, t2 - start and stop time of scan
%filename - filename of output map
%rxNum - receiver 0--4 (default 0)
%skip - set to skip non-existing data files 
%run - run number e.g. 'b2' 'k2014' changes behavior of pointing model for B2/Keck. Also changes source position for different Keck observations (height may be different)
%x_bin - [Optional] Change the pixelization of the output maps (default -2:.01:2)

if ~exist('rxNum','var') || length(rxNum) == 0
  rxNum=0;
end



type='sin'; % Is this right???? should be cos?? No difference between sin and cos - 2014-01-09 IDB

% Get the files:
files=list_arc_files('arc/',t1,t2);
n_files=length(files);
tic
% Concatenate files from beammap_demod_par
for fl=1:n_files
  disp(files{fl}(1,1:15))
  loadfile = strcat('scratch/beammaps/bm_',files{fl}(1,1:15));
  if run ~= 'b2'
    loadfile = [loadfile '_rx' num2str(rxNum)];
  end
  loadfile = [loadfile '.mat']; 
  if skip & ~exist(loadfile, 'file')
    disp(['Expected ' loadfile ' does not exist. SKIP']);
    continue;
  end
  e=load(loadfile);
  switch type % Throwout unnecessary fields to save memory
    case 'sin'
      e=rmfield(e.d,{'cos','fb'});
    case 'cos'
      e=rmfield(e.d,{'sin','fb'});      
    case 'quad'
      e=rmfield(e.d,'fb'); 
    case 'raw'
      e=rmfield(e.d,{'sin','cos'});
    case 'moon'
      e=e.d;
  end
  d(fl)=e;
end
toc

clear e

d=structcat(1,d);

nsamp=length(d.pos(:,1));

%[p ind]=get_array_info('20131231','obs');
[p ind]=get_array_info(files{1}(1,1:8),'obs');

%As of 2014-02-03 kbmp needs the following source information:
%    source.distance — The horizontal distance to the source from the zero point of the mount model (mount azimuth axis). For Keck observing the DSL mast, this is 211 meters.
%    source.azimuth — The azimuth of the source.
%    source.height — The height of the source, relative to the zero point of the mount model (floor of the Keck groundscreen). It is often convenient to calculate this from the source elevation: source.height = source.distance * tan(source_elevation). For the DSL mast, I have in the past used height = 13.2 meters (elevation = 3.58°).
%Note for B2 the origin of the mount model is the location of the focal plane

%source_lat = repmat(61.2,nsamp,1);
%source_lon = repmat(336,nsamp,1);
%Height of MAPO mast FROM base of groundshield = 470 in.
%Extra height of 12.5 + 96 in. from extension rod
%Total height = 578.5 in.
%horizontal Distance from azimuth axis to mast = 525 in.
%distance = 781.2 in. = 19.8

switch run
  case 'k2014'
% Source height for 2014 FSL mapping:
%  55' from source horn to MAPO roof. 54' from source horn to height of groundshield base 54' = 16.46 m
source.distance = 19.8;
source.azimuth = -18; % -22 deg. From http://bicep.caltech.edu/~spuder/keck_analysis_logbook/analysis/20120111_sidelobe_maps/ %Modified to center the source by eye
source.height = 16.46;
source.lateral_distance = 13.335; 
  case 'b2'
    %From util/sidelobe_parallax.m the source "distance" (approx from telescope to source) is 9.1 m
    %From Randol's BICEP2 reduc_sidelobemaps.m the source elevation is 61.2 deg.
    %So source.distance ~ 9.1 m * cos(61.2 deg)
    %Could get more accurate distances from DSL building drawings
    source.distance = 9.1 .* 0.48175;
    source.azimuth = 336; %From Randol's BICEP2 reduc_sidelobemaps.m
    source.height = 9.1 .* 0.87631;    
  otherwise
    error('unknown_run');
end

%Set mount model parameters
switch run
  case 'k2014'
%Keck: do nothing. kbmp has Keck mount as default.
    mount = [];
  case 'b2'
    %B2 mount parameters from http://bmode.caltech.edu/~spuder/keck_analysis_logbook/analysis/20121120_beam_map_pointing/
    mount.aperture_offz = 0.6;
    mount.aperture_offr = 0;
    mount.dk_offx = 0;
    mount.dk_offy = 0;
    mount.el_tilt = 0;
    mount.el_offx = 0;
    mount.el_offz = 0;
    mount.az_tilt_ha = 0;
    mount.az_tilt_lat = 0;
  otherwise
    error('unknown_run');
end

%x_bin=-180:1:180;
%y_bin=0:1:90;
%[xg yg]=meshgrid(x_bin,y_bin);
%[X Y Z]=sph2cart(xg*pi/180,yg*pi/180,ones(size(xg)));
if ~exist('x_bin','var') || length(x_bin) == 0
  x_bin = -2:.01:2;
end
y_bin = x_bin;

%keyboard

for ii=1:528
  ip = ii + 528 * rxNum; %Index into p
  if ~isnan(p.r(ip))
    pp=structcut(p,ip);
%    [phi theta pa_ap] =  keck_beam_map_pointing(d.pos(:,1), d.pos(:,2), ...
%      d.pos(:,3), [], [], source, pp,  'NoMirror', 'BeamCentered');
    [r_prime theta_prime psi] =  keck_beam_map_pointing(d.pos(:,1), d.pos(:,2), ...
      d.pos(:,3), mount, [], source, pp,  'NoMirror');
    [x_P, y_P] = rtheta_to_xbyb(r_prime, theta_prime);
%    [map(:,:,ii) map_std(:,:,ii) xi yi count_here] = (grid_map(phi.',theta.',d.sin(:,ii),x_bin,y_bin));
    [map(:,:,ii) map_std(:,:,ii) xi yi count_here] = (grid_map(x_P.',y_P.',d.sin(:,ii),x_bin,y_bin));
    map_count(:,:,ii) = reshape(count_here, length(yi), length(xi));
% [ZI SI XI YI I] = GRID_MAP(...) also returns the bin index
    [pvr(:,ii) S r]=binned_mean(r_prime.',d.sin(:,ii),200);
  end
end

if isempty(strmatch('sidelobemaps/',filename))
  filename=['sidelobemaps/' filename];
end

save(filename,'pvr','map','map_std','map_count', 'r', 'x_bin', 'y_bin', 'p', 'ind', 'mount', 'source', 'type', 'files')
