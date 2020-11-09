function bm_winfitmap(bmopt)
% bm_winfitmap(bmopt)
%
% Function to take the full binned map (primary output of bm_makemap),
% window around the main beam, and perform a 2D Gaussian fit.
%
% 20180614: in the beam map pipeline upgrade, this code remains
% mostly the same.  Optimized to take in maps in 
% xp/yp or ap_az/el coords and window without interpolation.
% In the case of xp/yp, have the option to move beam centroid
% to (0,0) after fitting.  This requires interpolation but is
% necessary to ensure beams are composited properly.
%
% We could in principle window in x/y, but that would add
% complications in not only this code but further downstream in
% the pipeline as well (each det would have different bins).
% We will never composite in x/y (distortion effects) and 
% we can extract pointing information from the unshifted
% xpyp_ideal window maps, so I don't see a need to add this
% functionality.
%
% INPUTS (should all be passed in with bmopt)
%
%   expt:         'bicep2','keck','bicep3';
%   fullmaptime:  timestamp of full map, e.g. '20121113_043509'
%                 Get this from get_bm_info
%
% OPTIONAL INPUTS
%  
%   mapcoord:   'azel_ap' -- window drawn around source_az,el
%               'xpyp' -- window around (0,0), apply origin shift 
%                         if requested (default)
%   shift:      'none' do not shift, keep ideal pointings
%               'ab_centroid' (default) shift pair centroid to origin
%               'detector' shift per-det beam center to origin
%   source_az:  only necessary for mapcoord = 'azel_ap'
%               Beam center in apparent az
%   source_el:  only necessary for mapcoord = 'azel_ap'
%               Beam center in apparent el
%   fullmapdir: directory of full map, default 'beammaps/maps'
%   demodtype:  'square' (default), 'choplet'
%   winmapsize: radius from center to window edge in deg, default 8
%   winmapdir:  directory to save in, default 'beammaps/maps_windowed'
%   buddy:      = 1 to load/window maps with buddy beam at center
%                 Only works with mapcoord = 'xpyp' so far
%   suffix:     optional string to add to map name (default '') 

% Parse bmopt
expt = bmopt.expt;
fullmaptime = bmopt.fullmaptime;
if ~isfield(bmopt,'mapcoord')
  mapcoord = 'xpyp';
  bmopt.mapcoord = mapcoord;
else
  mapcoord = bmopt.mapcoord;
end
if ~isfield(bmopt,'shift')
  shift = 'ab_centroid';
  bmopt.shift = shift;
else
  shift = bmopt.shift;
end
if ~isfield(bmopt,'source_az') && strcmp(mapcoord,'azel_ap')
  error(['Need source_az/el for windowing in apparent az,el'])
end
if ~isfield(bmopt,'source_el') && strcmp(mapcoord,'azel_ap')
  error(['Need source_az/el for windowing in apparent az,el'])
end
if ~isfield(bmopt,'fullmapdir')
  fullmapdir = 'beammaps/maps/';
  bmopt.fullmapdir = fullmapdir;
else
  fullmapdir = bmopt.fullmapdir;
end
if ~isfield(bmopt,'demodtype')
  demodtype = 'square';
  bmopt.demodtype = demodtype;
else
  demodtype = bmopt.demodtype;
end
if ~isfield(bmopt,'winmapsize')
  winmapsize = 8;
  bmopt.winmapsize = winmapsize;
else
  winmapsize = bmopt.winmapsize;
end
if ~isfield(bmopt,'winmapdir')
  winmapdir = 'beammaps/maps_windowed';
  bmopt.winmapdir = winmapdir;
else
  winmapdir = bmopt.winmapdir;
end
if ~isfield(bmopt,'buddy')
  buddy = 0;
  bmopt.buddy = buddy;
else
  buddy = bmopt.buddy;
  if buddy == 1 && ~strcmp(mapcoord,'xpyp')
    error('Can only make full buddy beam maps with mapcoord = xpyp.')
  end
end
if ~isfield(bmopt,'suffix')
  suffix = '';
  bmopt.suffix = suffix;
else
  suffix = bmopt.suffix;
end

bmopt

% Load up the map
filename = [fullmapdir '/map_' fullmaptime '_' demodtype ...
            '_' mapcoord suffix];
disp(['Loading ' filename '...'])
w = load(filename);

% Get map (masked and unmasked), bins, dk from fullmap file
map = double(w.bm.map);
map_masked = double(w.bm.map_masked);
if strcmp(expt,'keck')
  map_mirrormask = double(w.bm.map_mirrormask);
end
x_bin = w.bm.x_bin;
y_bin = w.bm.y_bin;
dk = w.bm.dk;
p = w.bm.p;
ind = w.bm.ind;

N_det = length(p.r);

% Find the map resolution
res = nanmedian(diff(x_bin)); % better be the same for az/el

% In bm pipeline re-write, much thought was put into binning strategy.
%
% X PRIME, Y PRIME
% Full map in xpyp are already all in the same bins, with a bin 
% centered at (0,0).  If 'xpyp_ideal' chosen, we just draw a window
% around origin and we're done, no interpolation.
% But for real 'xpyp' maps we need
% the beam centroid at the origin, so we'll do the fit, interp the 
% fit beam center to the origin, then window and re-fit.
%
% APPARENT AZ, EL
% Draw window around user-provided source location.  No interpolation. 

switch mapcoord
  case 'xpyp'
    cent_az = 0;
    cent_el = 0;
  case 'azel_ap'
    cent_az = source_az;
    cent_el = source_el;
end

%az_cut = (-winmapsize : res : winmapsize) + source_az;
%el_cut = (-winmapsize : res : winmapsize) + source_el;

% For Keck, find which dets have main beams masked out by mirror.
% Mark them, do not fit them.
% Later, will update compositing code such that these maps with
% sidelobe info but no main beam can be combined with other maps.
if strcmp(expt,'keck')
  onmirror = ...
      findbeamsonmirror(N_det,x_bin,y_bin,map_mirrormask,cent_az,cent_el);
else
  onmirror = ones(N_det,1);
end

% Get x_bin, y_bin from simple window around center, window map
l_bin = round(winmapsize * 2 / res);
az_cut = find(abs(x_bin-cent_az) <= winmapsize);
el_cut = find(abs(y_bin-cent_el) <= winmapsize);
x_bin = x_bin(az_cut);
y_bin = y_bin(el_cut);
mapcut = map(el_cut,az_cut,:);
mapcut_masked = map_masked(el_cut,az_cut,:);
  
% Fit the map
% Also only fit on inner 2 degrees because it takes FOREVER otherwise
% Don't bother fitting maps where beam was marked as off the mirror
disp('Fitting inner 2 deg of maps...')
A_preshift = fit_beamparams(x_bin,y_bin,mapcut,...
                            cent_az,cent_el,buddy,l_bin,onmirror);

% Recenter if desired, and re-fit.  
if ~strcmp(shift,'none')
  % Code adopted from ffbm_maskground:
  % Take a map of a source centered at arbitrary az/el and use the fits 
  % to shift the axes so (0,0) is at the pair centroid or beam
  % center
  mapcutnew = NaN(size(mapcut));
  for ii = 1:length(ind.la)
    fita = A_preshift(:,ind.la(ii));
    fitb = A_preshift(:,ind.lb(ii));
    % Don't shift if fit center is too far away (bad fit/beam)
    % or if any fits are NaN
    if sqrt(fita(2).^2 + fita(3).^2) > 8 || ...
        sqrt(fitb(2).^2 + fitb(3).^2) > 8 || ...
        any(isnan(fita)) || any(isnan(fitb))
      mapcutnew(:,:,ind.la(ii)) = mapcut(:,:,ind.la(ii));
      mapcutnew(:,:,ind.lb(ii)) = mapcut(:,:,ind.lb(ii));
      continue
    end
    [xx_a yy_a xx_b yy_b] = shiftaxes(fita,fitb,shift,x_bin,y_bin);
    % Interpolate map to new center
    [X,Y] = meshgrid(x_bin,y_bin);
    [XX_a,YY_a] = meshgrid(xx_a,yy_a);
    [XX_b,YY_b] = meshgrid(xx_b,yy_b);
    mapcutnew(:,:,ind.la(ii)) = ...
        interp2(XX_a,YY_a,mapcut(:,:,ind.la(ii)),X,Y);
    mapcutnew(:,:,ind.lb(ii)) = ...
        interp2(XX_b,YY_b,mapcut(:,:,ind.lb(ii)),X,Y);
  end
  % Re-fit
  mapcut = mapcutnew;
  mapcutnew_masked = mapcutnew;
  mapcutnew_masked(isnan(mapcut_masked(:))) = NaN;
  mapcut_masked = mapcutnew_masked;
  A = fit_beamparams(x_bin,y_bin,mapcutnew,...
                     cent_az,cent_el,buddy,l_bin,onmirror);
else
  A = A_preshift;
end

% Prep output structure.  ffbm_makesinglecompfile wants:
% map, x_bin, y_bin, A
bm.map = mapcut;
bm.map_masked = mapcut_masked;
bm.x_bin = x_bin;
bm.y_bin = y_bin;
bm.A = A;
if strcmp(mapcoord,'xpyp') && ~strcmp(shift,'none')
  bm.A_preshift = A_preshift;
end
bm.onmirror = onmirror;
bm.dk = dk;

% Save
if strcmp(mapcoord,'xpyp') && strcmp(shift,'none')
  mapcoord = 'xpyp_ideal';
end
sizestr = [strrep(num2str(winmapsize),'.','p') 'deg'];
savename = [winmapdir '/mapwin_' fullmaptime '_' demodtype ...
            '_' mapcoord '_' sizestr suffix];
disp(['Saving windowed map:' savename]);
save(savename,'bm');

return


%%%%%%%%%%%%%%%%%%%%%% SUB-FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%

function onmirror = ...
    findbeamsonmirror(N_det,x_bin,y_bin,map,cent_az,cent_el)
% Take masked windowed map and find which dets have main beams
% that were cut by the mirror mask.  
% N_det length vector, = 1 on mirror, = 0 off mirror
onmirror = ones(N_det,1);
temp_az_cut = find(abs(x_bin-cent_az) <= 1.5);
temp_el_cut = find(abs(y_bin-cent_el) <= 1.5);
for idet = 1:N_det
  thismap = map(temp_el_cut,temp_az_cut,idet);
  % If it's all NaN, it's probably an off detector 
  if numel(find(isnan(map(:,:,idet)))) == numel(map(:,:,idet))
    continue
  elseif any(isnan(thismap(:)))
    onmirror(idet) = 0;
  end
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function A = ...
    fit_beamparams(x_bin,y_bin,map,cent_az,cent_el,buddy,l_bin,onmirror)

% Turn off all warnings which swamp the output file
warning('off','all') 

% Only fit inner 2 degrees
fit_az_cut = find(abs(x_bin-cent_az) <= 2);
fit_el_cut = find(abs(y_bin-cent_el) <= 2);

A = NaN(7,size(map,3));
% Don't use initial guesses provided by default in normfit2d,
% since they depend on window size.  But funny stuff happens
% if any guesses are identically zero, so make the epsilon instead.
switch buddy
  case 0
    A0 = [300 cent_az cent_el 0.15 0.15 1e-5 0];
  case 1
    A0 = [15 cent_az cent_el 0.15 0.15 1e-5 0];
end
A0(A0==0) = A0(A0==0) + 1e-5;

% Use optimset to increase tolerance during fitting
opts = optimset('TolFun',1e-8,'TolX',1e-8,'MaxFunEvals',7000,'MaxIter',7000);

for ii = 1:size(map,3)
  thismap = map(:,:,ii);
  thismap = thismap(fit_el_cut,fit_az_cut);
  % If map is all NaN or if main beam is off mirror, skip the fitting
  if numel(find(isnan(thismap))) == numel(thismap) || ...
        onmirror(ii) == 0
    continue
  else
    thismap = inpaint_nans(thismap);
    xfit = x_bin(fit_az_cut);
    yfit = y_bin(fit_el_cut);
    disp(ii)
    A(:,ii) = normfit2d(xfit,yfit,thismap,'A0',A0,'opts',opts);
  end
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [xx_a yy_a xx_b yy_b] = shiftaxes(fita,fitb,cent,x_bin,y_bin)
% Take a map of a source centered at arbitrary az/el and use the fits to
% shift the axes so (0,0) is at the pair centroid or beam center

% Find centroid or beam center
% Do A/B both exist? Just check X component of fit
if ~isnan(fita(2)) & ~isnan(fitb(2)) % Both good, use centroid
  switch cent
    case 'ab_centroid'
      xcent_a = (fita(2) + fitb(2))./2;
      xcent_b = (fita(2) + fitb(2))./2;
      ycent_a = (fita(3) + fitb(3))./2;
      ycent_b = (fita(3) + fitb(3))./2;
    case 'detector'
      xcent_a = fita(2);
      xcent_b = fitb(2);
      ycent_a = fita(3);
      ycent_b = fitb(3);
  end
else % For these, 'cent' doesn't matter
  if ~isnan(fita(2)) % Just A exists, use A
    xcent_a = fita(2);
    xcent_b = fita(2);
    ycent_a = fita(3);
    ycent_b = fita(3);
  elseif ~isnan(fitb(2)) % Just B exists, use B
    xcent_a = fitb(2);
    xcent_b = fitb(2);
    ycent_a = fitb(3);
    ycent_b = fitb(3);
  else % Neither exists!  
    xcent_a = 0;
    xcent_b = 0;
    ycent_a = 0;
    ycent_b = 0;
  end
end

% Slide the axes over
xx_a = x_bin - xcent_a;
xx_b = x_bin - xcent_b;
yy_a = y_bin - ycent_a;
yy_b = y_bin - ycent_b;

% Move fit values too
%fita(2) = fita(2) - xcent_a;
%fitb(2) = fitb(2) - xcent_b;
%fita(3) = fita(3) - ycent_a;
%fitb(3) = fitb(3) - ycent_b;

return
