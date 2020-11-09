function ffbm_makesinglecompfile(compopt)
% function ffbm_makesinglecompfile(compopt)
%
% Function to save all of the component maps we want to use into a single
% file (instead of having to load up 50+ maps).  Used as input to
% make composite maps.
%
% Beam map pipeline upgrade: the windowed maps have already been
% masked and are already in xp/yp coordinates.  Now we can go
% straight from component --> composite maps.  
%
% Option to include only masked or only unmasked maps in
% component file (both in one file would make file way too large). 
% Masked maps will be used in compositing.
%
% Must have maps windowed in xp/yp.  Saves in xp/yp.
% 
% INPUTS (should all be sent in with compopt)
%
%   expt:          'bicep2','keck','bicep3'
%   year:          bicep2: 2012
%                  keck: 2012-2018
%                  bicep3: 2016-2018
%
% OPTIONAL INPUTS 
%
%   mask:            = 'masked' (default) include masked maps only
%                    = 'unmasked' only include unmasked maps
%   demodtype:       'square' (default), 'choplet'
%   mapsize:         map size in deg (8 default), used for map filename
%   suffix:          optional string to add to map name, '' default
%   bmdir:           directory from which to load beammaps
%                    (default beammaps/maps_windowed)
%   mapstoload:      vector of beam map run numbers to use, i.e. 1:10
%                    (default is everything in bmrunlist.csv)
%   componentdir:    directory in which to save output map struct
%                    (default beammaps/maps_component)
%   componentfile:   name of file to save in componentdir
%                    (default ffbm_year_allcomp_size_un/masked.mat)

% Parse compopt
expt = compopt.expt;
year = compopt.year;
if ~isfield(compopt,'mask');
  mask = 'masked';
else
  mask = compopt.mask;
end
if ~isfield(compopt,'demodtype');
  demodtype = 'square';
else
  demodtype = compopt.demodtype;
end
if ~isfield(compopt,'mapsize');
  sizestr = '8deg';
else
  sizestr = [strrep(num2str(compopt.mapsize),'.','p') 'deg'];
end
if ~isfield(compopt,'suffix')
  suffix = '';
else
  suffix = compopt.suffix;
end
if ~isfield(compopt,'bmdir')
  bmdir = 'beammaps/maps_windowed/';
else
  bmdir = compopt.bmdir;
end
% Which maps do we want?
% Get beam map run info
bm = get_bm_info(year);
if ~isfield(compopt,'mapstoload')
  compopt.mapstoload = bm.number;
else
  bm = structcut(bm,compopt.mapstoload);
end
if ~isfield(compopt,'componentdir')
  compopt.componentdir = 'beammaps/maps_component';
  componentdir = compopt.componentdir;
else
  componentdir = compopt.componentdir;
end
if ~isfield(compopt,'componentfile')
  compopt.componentfile = ['ffbm_' num2str(year) '_allcomp'...
	'_' sizestr '_' mask suffix];
  componentfile = compopt.componentfile;
else
  componentfile = compopt.componentfile;
end

[p ind] = get_array_info([num2str(year) '0201']);

% Load up each beam map.  Only use windowed map!
for ii = 1:length(bm.number)

  filename = bm.filename{ii};
  disp(['Loading map ' num2str(ii) ': ' filename]);
  wmapname = [bmdir,'/mapwin_',filename,'_',demodtype,...
	'_xpyp_',sizestr,suffix];
  w = load(wmapname);
  w = w.bm;
  
  % Save the maps in a larger array
  if strcmp(mask,'masked')
    component{ii} = w.map_masked;
  elseif strcmp(mask,'unmasked')
    component{ii} = w.map;
  else
    error('compopt.mask needs to be ''masked'' or ''unmasked''')
  end
  fit{ii} = w.A;
  x_bin{ii} = w.x_bin;
  y_bin{ii} = w.y_bin;
  
end

% Rearrange so we're in a more logical structure: for each detector, have an
% array of each map taken and the fit values
for ii = 1:length(p.gcp)
  for jj = 1:length(component)
    map{ii}.component{jj} = component{jj}(:,:,ii);
    map{ii}.fit{jj} = fit{jj}(:,ii);
  end
end

field_size_x = range(w.x_bin) + (w.x_bin(2) - w.x_bin(1));
field_size_y = range(w.y_bin) + (w.y_bin(2) - w.y_bin(1));

% Sometimes the field sizes are slightly off (16.0999, 16.10001) etc.
% This can screw with the ad calculation - here we use the world's
% jankiest method of rounding to 2 decimal places (since we are never
% going higher resolution than that), because we're stuck with R2009 - if
% we had 2014 or higher we could use round(x,n)...
field_size_x = str2num(sprintf('%.2f',field_size_x));
field_size_y = str2num(sprintf('%.2f',field_size_y));

% Create ad struct.  
% comp.ad.t_val_deg{1} = increasing x prime 
% comp.ad.t_val_deg{2} = decreasing y prime
% To make our standard-parity beam plot, we do:
%   >> imagesc(ad.t_val_deg{2},ad.t_val_deg{1},map);
%   >> set(gca,'ydir','normal','xdir','reverse')
%   >> xlabel('y prime'); ylabel('x prime')
% When we used to bin in apparent az/el then rotate, we had
% to be careful to get this exactly right.  Binning in xp/yp
% eliminates this confusion; ad struct as is will match above convention.

comp.ad = calc_ad2([field_size_x field_size_y],...
    [length(w.x_bin) length(w.y_bin)]);
comp.map = map;
comp.bm = bm;
comp.x_bin = x_bin;
comp.y_bin = y_bin;
comp.compopt = compopt;

disp('Saving ...')
savename = [componentdir '/' componentfile];
save(savename,'comp','-v7.3');
disp(['Saved single component file: ' savename]);

return