function bm_makemap(bmopt)
% bm_makemap(bmopt)
% 
% Function to take the demodulated beam map timestreams
% and bin them into maps, with masking calculated from
% bm_timestream.
%
% Note that the coordinate systems are now calculated upstream in 
% bm_timestream (horizon, apparent az/el, x/y, xp/yp).  Choice of binning
% must be one of those four.
% 
% Note that maps are binned into pixels 2*el_res (deg)
%
% INPUTS (should all be passed in with bmopt)
%
%   expt:        'bicep2','keck','bicep3'
%   year:        Bicep2: 2012
%                Keck: 2012+
%                Bicep3: 2016+
%   number:      which run # from this FFBM campaign? From runfile
%   run:         BICEP2 runs: 'r5','r8'
%                Keck runs:   'k0r4', 'k0r4t2' (tile 2 only), highbay
%                             'k083' = k7 run 2013-10, highbay
%                             'k201?' = Keck ffbm @ Pole 
%                             'k201?fsl' = Keck sidelobe @ Pole (13,14)
%                BICEP3 runs: 'hb3r5' = B3 run 5, highbay
%                             'b3r6' = B3 run 6 @ Pole, 2015-02
%                             'b3r8' = B3 run 8 @ Pole, 2016-02
%                             'b3r9' = B3 run 9 @ Pole. 2017/18-02
%
% OPTIONAL
%
%   demodtype:   'square' (default), 'choplet' 
%                (for filenames of timestream data to load)
%   demodshift:  if used manual demodshift in bm_timestream, put that here
%                if used a default value, leave this empty
%                (for filenames of timestream data to load)
%   timestreamdir: directory from which to load timestream data
%                  default 'beammaps/timestream'
%   fullmapdir:  directory in which to save beammaps
%                default beammaps/maps
%   el_res:      Resolution of elevation steps (deg), 0.05 default
%                Pixel size of maps is this x 2
%   xpypmax:     For maps binned into xpyp, this sets the range on 
%                the binning.  Bins will be
%                -xpypmax : el_res*2 : xpypmax (default 40)
%                This forces the existence of a bin centered at
%                (0,0) as long as xpypmax is divisible by el_res*2
%   component:   Specify demodulation component to map. options are:
%                'cos':  cos demodulation, default
%                'sin':  sin demodulation
%                'comp': complex: re((cos+i*sin)*exp(i*theta*(pi/180))) - 
%                        send in theta or defaults to 0, i.e. cos
%                'quad': quadrature sum, sqrt( cos.^2 + sin.^2 )
%                'phase':arctan(sin./cos)
%                'raw':  raw fb, no demodulation
%   theta:       angle by which to rotate complex map, 0 default (degrees)
%   mapcoord:    Which coordinate system to use for binning
%                Note: finding the coords is done in bm_timestream
%                so it must be one of these options from that step:
%                'azel_ap': apparent az/el per det, behind mirror
%                'hor': mount horizon coords 
%                'xy': focal plane centered coordinates
%                'xpyp': coords centered on ideal det center
%   buddy:       = 1 if you want to make map in xp/yp coords
%                    with buddy beam at (0,0) instead of main beam.
%                    Only works with mapcoord = "xpyp".
%                = 0 default
%   pairdiffsum: = 0 do individual detector maps (default)
%                = 1 load pair diffed/summed timestreams and bin those
%   suffix:      optional string to add to end of map name (default '')

% Parse bmopt 
expt = bmopt.expt;
year = bmopt.year;
number = bmopt.number;
run = bmopt.run;

pp = get_bm_info(year);
t1 = pp.t1{number};
t2 = pp.t2{number};
if ~isfield(bmopt,'demodtype')
  demodtype = 'square';
  bmopt.demodtype = demodtype;
else
  demodtype = bmopt.demodtype;
end
if ~isfield(bmopt,'demodshift')
  demodshift = [];
  bmopt.demodshift = demodshift;
else
  demodshift = bmopt.demodshift;
end
if ~isfield(bmopt,'timestreamdir')
  timestreamdir = 'beammaps/timestream';
  bmopt.timestreamdir = timestreamdir;
else
  timestreamdir = bmopt.timestreamdir;
end
if ~isfield(bmopt,'fullmapdir')
  fullmapdir = 'beammaps/maps';
  bmopt.fullmapdir = fullmapdir;
else
  fullmapdir = bmopt.fullmapdir;
end
if ~isfield(bmopt,'el_res')
  el_res = 0.05;
  bmopt.el_res = el_res;
else
  el_res = bmopt.el_res;
end
if ~isfield(bmopt,'xpypmax')
  xpypmax = 40;
  bmopt.xpypmax = xpypmax;
else
  xpypmax = bmopt.xpypmax;
end
if ~isfield(bmopt,'component')
  component = 'cos';
  bmopt.component = component;
else
  component = bmopt.component;
end
if ~isfield(bmopt,'theta')
  theta = 0;
  bmopt.theta = theta;
else
  theta = bmopt.theta;
end
if ~isfield(bmopt,'mapcoord')
  mapcoord = 'xpyp';
  bmopt.mapcoord = mapcoord;
else
  mapcoord = bmopt.mapcoord;
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
if ~isfield(bmopt,'pairdiffsum')
  pairdiffsum = 0;
  bmopt.pairdiffsum = 0;
else
  pairdiffsum = bmopt.pairdiffsum;
end
if ~isfield(bmopt,'suffix')
  suffix = '';
else
  suffix = bmopt.suffix;
end

if strcmp(mapcoord,'xpyp') && mod(xpypmax,el_res*2) > 0
  error(['For xpyp maps we need xpypmax divisible by el_res*2 to' ...
        ' enforce a bin a (0,0)'])
end

bmopt

% Get the files:
files = list_arc_files('arc/',t1,t2);
n_files = length(files);

[p ind] = get_array_info(files{1}(1:8),'ideal');

% If there is a user-input demod shift, make sure the demodded arcfile has
% it in the filename
if ~isempty(demodshift)
  demodstr = ['_shift' int2str(demodshift)];
else
  demodstr = [];
end

% If doing pair summed/diffed maps, include in filename
if pairdiffsum
  pairstr = 'pair_';
else
  pairstr = '';
end

% Choose the first arcfile as the beammap name
kk = 0;
while ~exist([timestreamdir,'/bm_',pairstr,files{kk+1}(1,1:15),...
              '_',demodtype,demodstr,'.mat'],'file')
  kk = kk + 1;
end
filename = files{kk+1}(1,1:15);

% Warn user if map already exists
flstr = [fullmapdir,'/','map_',filename,'_',demodtype,'_',...
         mapcoord,suffix,demodstr,'.mat'];
if exist_file(flstr) 
  disp(sprintf('Map %s already exists, will be overwritten ',flstr));
end

% Concatenate timestream files
disp('### Concatenating timestream files')
tic
for fl = 1:n_files
  disp(files{fl}(1,1:15))
  % Does the file exist?  
  timestreamname = [timestreamdir,'/bm_',pairstr,files{fl}(1,1:15),...
               '_',demodtype,demodstr,'.mat'];
  if exist(timestreamname,'file')
    e = load([timestreamdir,'/bm_',pairstr,files{fl}(1,1:15),...
	      '_',demodtype,demodstr,'.mat']);
  
    % Keep only the component we're plotting
    switch component
      case 'sin'
        if pairdiffsum
          e = rmfield(e.d,{'cos_sum','cos_diff','fb_sum','fb_diff'});
        else
          e = rmfield(e.d,{'cos','fb'});
        end
      case 'cos'
        if pairdiffsum
          e = rmfield(e.d,{'sin_sum','sin_diff','fb_sum','fb_diff'});
        else
          e = rmfield(e.d,{'sin','fb'});
        end
      case {'comp','quad','phase'}
        if pairdiffsum
          e = rmfield(e.d,{'cos_sum','fb_sum','fb_diff'});
        else
          e = rmfield(e.d,'fb');  
        end
      case 'raw'
        if pairdiffsum
          e = rmfield(e.d,{'cos_sum','cos_diff','sin_sum','sin_diff'});
        else
          e = rmfield(e.d,{'sin','cos'});
        end
    end
    d(fl) = e;
  end
end
disp(['Done concatenating (' num2str(toc) ' seconds)']);

clear e
d = structcat(1,d);
d.pointing = structcat(1,d.pointing);

% Avg dk angle through schedule.  B3 needs sign flip.
switch expt
  case 'bicep3'
    dk = -nanmean(d.pointing.dk);
  otherwise
    dk = nanmean(d.pointing.dk);
end

% Gather sin demod, cos demod, or quadrature sum:
componentstr = ['_' component]; 
[d N_pts N_maps] = getamp(d,component,theta,pairdiffsum);

% Which coord do we want to use?
mapperrx = N_maps / length(unique(p.rx));
switch mapcoord
  case 'azel_ap'
    az_coord = d.pointing.az_ap;
    el_coord = d.pointing.el_ap;
  case 'xy'
    % For x/y or xp/yp: Loading the beam map then plotting
    %   >> imagescnan(bm.x_bin,bm.y_bin,bm.map(:,:,ii))
    %   >> set(gca,'YDir','normal','XDir','reverse')
    % should give the standard parity beam
    switch expt
      case 'keck'
        az_coord = NaN([N_pts N_maps]);
        el_coord = NaN([N_pts N_maps]);
        for ii = 1:length(unique(p.rx))
          az_coord(:,(1:mapperrx)+(ii-1)*mapperrx) = ...
              repmat(d.pointing.y(:,ii),1,mapperrx);
          el_coord(:,(1:mapperrx)+(ii-1)*mapperrx) = ...
              repmat(d.pointing.x(:,ii),1,mapperrx);
          maskmirror(:,(1:mapperrx)+(ii-1)*mapperrx) = ...
              repmat(d.maskmirror(:,ii),1,mapperrx);
        end
      otherwise
        az_coord = repmat(d.pointing.y,1,N_maps);
        el_coord = repmat(d.pointing.x,1,N_maps);
    end
  case 'xpyp'
    if buddy == 1
      az_coord = d.pointing.yp_buddy;
      el_coord = d.pointing.xp_buddy;
    else
      az_coord = d.pointing.yp;
      el_coord = d.pointing.xp;
    end
  case 'hor'
    az_coord = repmat(d.pointing.az,1,N_maps); 
    el_coord = repmat(d.pointing.el,1,N_maps); 
end
% Make mirror mask same size as d.amp matrix
maskmirror = NaN([N_pts N_maps]);
for ii = 1:length(unique(p.rx))
  maskmirror(:,(1:mapperrx)+(ii-1)*mapperrx) = ...
      repmat(d.maskmirror(:,ii),1,mapperrx);
end

% Make masks double and set 0 --> NaN
maskmirror(maskmirror==0) = NaN;
maskground = double(d.maskground);
maskground(maskground==0) = NaN;

%%%%%%%%%%%%%%%%
% MAKE THE MAP
disp('### Making the map...');
if pairdiffsum
  % Sum
  amp = d.amp_sum;
  d.amp = amp;
  [map_sum x_bin y_bin] = ...
      create_map(d,ind,el_res,az_coord,el_coord,mapcoord,xpypmax);
  % Sum, ground+mirror mask
  d.amp = amp .* maskground .* maskmirror;
  [map_sum_masked x_bin y_bin] = ...
      create_map(d,ind,el_res,az_coord,el_coord,mapcoord,xpypmax);
  % Sum, mirror mask only
  if strcmp(expt,'keck')
    d.amp = amp .* maskmirror;
    [map_sum_mirrormask x_bin y_bin] = ...
      create_map(d,ind,el_res,az_coord,el_coord,mapcoord,xpypmax);
  end
  % Diff
  amp = d.amp_diff;
  d.amp = amp;
  [map_diff x_bin y_bin] = ...
      create_map(d,ind,el_res,az_coord,el_coord,mapcoord,xpypmax);
  % Diff, ground+mirror mask
  d.amp = amp .* maskground .* maskmirror;
  [map_diff_masked x_bin y_bin] = ...
      create_map(d,ind,el_res,az_coord,el_coord,mapcoord,xpypmax);
  % Diff, mirror mask only
  if strcmp(expt,'keck')
    d.amp = amp .* maskmirror;
    [map_diff_mirrormask x_bin y_bin] = ...
      create_map(d,ind,el_res,az_coord,el_coord,mapcoord,xpypmax);
  end 
  bm.map_sum = map_sum;
  bm.map_diff = map_diff;
  bm.map_sum_masked = map_sum_masked;
  bm.map_diff_masked = map_diff_masked;
  if strcmp(expt,'keck')
    bm.map_sum_mirrormask = map_sum_mirrormask;  
    bm.map_diff_mirrormask = map_diff_mirrormask;
  end
else % per detector maps
  [map x_bin y_bin] = ...
      create_map(d,ind,el_res,az_coord,el_coord,mapcoord,xpypmax);
  % Per-det, ground + mirror mask
  amp = d.amp;
  d.amp = amp .* maskground .* maskmirror;
  [map_masked x_bin y_bin] = ...
      create_map(d,ind,el_res,az_coord,el_coord,mapcoord,xpypmax);
  % Per-det, mirror mask only
  if strcmp(expt,'keck')
    d.amp = amp .* maskmirror;
    [map_mirrormask x_bin y_bin] = ...
        create_map(d,ind,el_res,az_coord,el_coord,mapcoord,xpypmax);
  end
  bm.map = map;
  bm.map_masked = map_masked;
  if strcmp(expt,'keck')
    bm.map_mirrormask = map_mirrormask;
  end
end

%%%%%%%%%%%%%%%%
disp(['Done binning (' num2str(toc) ' seconds)']);

% Finish output struct
bm.x_bin = x_bin;
bm.y_bin = y_bin;
%bm.source_az = source_az;
%bm.source_el = source_el;
bm.dk = dk;
bm.bmopt = bmopt;
bm.p = p;
bm.ind = ind;

disp('### Saving the full map...');
savename = [fullmapdir,'/map_',pairstr,filename,'_',demodtype,...
            '_',mapcoord,suffix,demodstr];
save(savename,'bm','-v7.3');
disp(sprintf('Saved full map %s',savename));  

system(['chmod g+w ' savename '.mat']);

return


%%%%%%%%%%%%%%%%%%%%%% SUB-FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%

function [d N_pts N_maps] = getamp(d,component,theta,pairdiffsum)

% Keep only the component we want.  Different case 
% for pair diff/sum vs single det.
% Also extract # points in the timestream as well as number
% of maps (which will be # dets or # light pairs)
if pairdiffsum
  switch component
    case 'sin' 
      d.amp_sum = d.sin_sum;
      d.amp_diff = d.sin_diff;
      d = rmfield(d,{'sin_sum','sin_diff'});
    case 'cos'
      d.amp_sum = d.cos_sum;
      d.amp_diff = d.cos_diff;
      d = rmfield(d,{'cos_sum','cos_diff'});
    case 'quad'
      d.amp_sum = sqrt(d.cos_sum.^2 + d.sin_sum.^2);
      d.amp_diff = sqrt(d.cos_diff.^2 + d.sin_diff.^2);
      d = rmfield(d,{'cos_sum','cos_diff','sin_sum','sin_diff'});
    case 'comp'
      d.amp_sum = real(complex(d.cos_sum,d.sin_sum).*exp(1i*theta*(pi/180)));
      d.amp_diff = real(complex(d.cos_diff,d.sin_diff).*exp(1i*theta*(pi/180)));
      d = rmfield(d,{'cos_sum','cos_diff','sin_sum','sin_diff'});
    case 'phase'
      d.amp_sum = atan2(d.sin_sum,d.cos_sum);
      d.amp_diff = atan2(d.sin_diff,d.cos_diff);
      d = rmfield(d,{'cos_sum','cos_diff','sin_sum','sin_diff'});
    case 'raw'
      d.amp_sum = d.fb_sum;
      d.amp_diff = d.fb_diff;
      d = rmfield(d,{'fb_sum','fb_diff'});
    case 'moon'
      d.amp_sum = d.fb_sum;
      d.amp_diff = d.fb_diff;
      source_az = 0;
      source_el = 0;
      d = rmfield(d,{'fb_sum','fb_diff'});    
  end
  N_pts = size(d.amp_sum,1);
  N_maps = size(d.amp_sum,2);
else
  switch component
    case 'sin' 
      d.amp = d.sin;
      d = rmfield(d,'sin');
    case 'cos'
      d.amp = d.cos;
      d = rmfield(d,'cos');
    case 'quad'
      d.amp = sqrt(d.cos.^2 + d.sin.^2);
      d = rmfield(d,{'cos','sin'});
    case 'comp'
      d.amp = real(complex(d.cos,d.sin).*exp(1i*theta*(pi/180)));
      d = rmfield(d,{'cos','sin'});
    case 'phase'
      d.amp = atan2(d.sin,d.cos);
      d = rmfield(d,{'cos','sin'});
    case 'raw'
      d.amp = d.fb;
      d = rmfield(d,'fb');
    case 'moon'
      d.amp = d.fb;
      source_az = 0;
      source_el = 0;
      d = rmfield(d,'fb');    
  end
  N_pts = size(d.amp,1);
  N_maps = size(d.amp,2);
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [map x_bin y_bin] = ...
    create_map(d,ind,el_res,az_coord,el_coord,mapcoord,xpypmax)

% Choose az/el bins
[az_coord,el_coord] = remove_outliers(az_coord,el_coord,ind);
switch mapcoord
  case 'xpyp'
    x_bin = -xpypmax : el_res*2 : xpypmax;
    y_bin = -xpypmax : el_res*2 : xpypmax;    
  otherwise
    x_bin = min(min(az_coord(:,ind.l))) : el_res*2 : ...
            max(max(az_coord(:,ind.l)));
    y_bin = min(min(el_coord(:,ind.l))) : el_res*2 : ...
            max(max(el_coord(:,ind.l)));    
end

map = NaN(length(y_bin),length(x_bin),size(d.amp,2)); 

for ii = 1:size(d.amp,2); 
  % Generate x and y map coordinates
  x = az_coord(:,ii);
  y = el_coord(:,ii);
  % Actual binning
  map(:,:,ii) = grid_map(x,y,d.amp(:,ii),x_bin,y_bin);
end

% Convert double --> single, save space
map = single(map);

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [az el] = remove_outliers(az, el, ind)
% Simple function to exclude outliers from az and el coordinates
% This can be removed once the skipped samples are excluded in bm_demod
% This assumes outliers are common to all the channels
% so we choose 1st ind.l channel

ch = ind.l(1);
q =  find(abs(az(:,ch) - median(az(:,ch))) > 5*std(az(:,ch)) | ...
    abs(el(:,ch) - median(el(:,ch))) > 5*std(el(:,ch)));
disp(q)
if ~isempty(q)
  az(q,:) = NaN;
  el(q,:) = NaN;
end

return

