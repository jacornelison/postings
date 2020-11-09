function ffbm_makebl(compopt)
% function ffbm_makebl(compopt)
%
% Take composite maps and 
% - Make a stacked beam per-rx or per-frequency
% - Calculate the radial profile
% - Calculate the B_l window function
% - Save and plot the above  (stacked map contains map/radial profile)
%
% This works best with composites made out to large (e.g. 8 deg) radii.
%
% INPUTS (should all be sent in with compopt)
%
%  expt:          'keck','bicep3'
%  year:          keck: 2014, 2015
%                 bicep3: 2015
%
% OPTIONAL INPUTS
%
%   compositedir:  directory from which to load composite maps
%                  (default beammaps/maps_composite)
%   compositefile: file in compositedir to load
%                  (default ffbm_year_all_8deg_xpyp)
%   suffix:        optional string to add to map to load, '' default
%   stackmapdir:   directory in which to save stacked map
%                  (default beammaps/averaged_maps)
%   stackmapfile:  filename of stacked map
%                  (default averagedmap_year_freq/rx.mat)
%   bldir:         directory in which to save B_l file
%                  (default beammaps/beamparams)
%   blfile:        filename to save
%                  (default beamfile_year0101_sum_(freq/rx).fits)
%   perwhat:       how to separate detectors
%                  'perfreq' (default),'perrx'
%   onlychflag:    1 (default) to only coadd dets which pass chflags
%   method:        how to calc stacked map - for each pixel in the map, use
%                  'median' (default),'mean'
%   weight:        how to weight the detectors (for 'mean' above)
%                  'equal' (default),'std','cmb' (if they exist)
%   chopper:       chopper diameter (inches) for B_l correction
%                  0 (default for pre-2013)
%                  18 (default for 2013-2015)
%                  24 (default for 2016-2017)
%   norm:          how to make the beam and radial profile plots
%                  'peak' (default), 'dBi'
%   deoutlier:     1 (default) to remove spikes that may have survived composite
%   plotbeam:      1 to plot the stacked beam (default)
%   plotprofile:   1 to plot the radial profile (default)
%   plotbl:        1 to plot the B_l

% Parse compopt
expt = compopt.expt;
year = compopt.year;
if ~isfield(compopt,'suffix')
  suffix = '';
else
  suffix = compopt.suffix;
end
if ~isfield(compopt,'compositedir')
  compopt.compositedir = 'beammaps/maps_composite';
  compositedir = compopt.compositedir;
else
  compositedir = compopt.compositedir;
end
if ~isfield(compopt,'compositefile')
  compopt.compositefile = ['ffbm_' num2str(year) '_all'...
	'_8deg_xpyp' suffix];
  compositefile = compopt.compositefile;
else
  compositefile = compopt.compositefile;
end
if ~isfield(compopt,'stackmapdir')
  compopt.stackmapdir = 'beammaps/averaged_maps';
  stackmapdir = compopt.stackmapdir;
else
  stackmapdir = compopt.stackmapdir;
end
if ~isfield(compopt,'stackmapfile')
  compopt.stackmapfile = ['averagedmap_' num2str(year) '_'];
  stackmapfile = compopt.stackmapfile;
else
  stackmapfile = compopt.stackmapfile;
end
if ~isfield(compopt,'bldir')
  compopt.bldir = 'beammaps/bl';
  bldir = compopt.bldir;
else
  bldir = compopt.bldir;
end
if ~isfield(compopt,'blfile')
  compopt.blfile = ['beamfile_' num2str(year) '0101_sum_' ];
  blfile = compopt.blfile;
else
  blfile = compopt.blfile;
end
if ~isfield(compopt,'perwhat')
  compopt.perwhat = 'perfreq';
  perwhat = compopt.perwhat;
else
  perwhat = compopt.perwhat;
end
if ~isfield(compopt,'onlychflag')
  compopt.onlychflag = 1;
  onlychflag = compopt.onlychflag;
else
  onlychflag = compopt.onlychflag;
end
if ~isfield(compopt,'method')
  compopt.method = 'median';
  method = compopt.method;
else
  method = compopt.method;
end
if ~isfield(compopt,'weight')
  compopt.weight = 'equal';
  weight = compopt.weight;
else
  weight = compopt.weight;
end
if ~isfield(compopt,'chopper')
  switch year
    case {2016,2017}
      compopt.chopper = 24;
    case {2013,2014,2015}
      compopt.chopper = 18;
    otherwise
      compopt.chopper = 0;
  end
  chopper = compopt.chopper;
else
  chopper = compopt.chopper;
end
if ~isfield(compopt,'norm')
  compopt.norm = 'peak';
  norm = compopt.norm;
else
  norm = compopt.norm;
end
if ~isfield(compopt,'deoutlier')
  compopt.deoutlier = 1;
  deoutlier = compopt.deoutlier;
else
  deoutlier = compopt.deoutlier;
end
if ~isfield(compopt,'plotbeam')
  compopt.plotbeam = 1;
  plotbeam = compopt.plotbeam;
else
  plotbeam = compopt.plotbeam;
end
if ~isfield(compopt,'plotprofile')
  compopt.plotprofile = 1;
  plotprofile = compopt.plotprofile;
else
  plotprofile = compopt.plotprofile;
end
if ~isfield(compopt,'plotbl')
  compopt.plotbl = 1;
  plotbl = compopt.plotbl;
else
  plotbl = compopt.plotbl;
end

% Get chflags if required
chflags = [];
if onlychflag
  chflags = get_default_chflags(expt,year);
end

[p ind] = get_array_info([num2str(year) '0201'],...
                         [],[],[],[],chflags);
p = rmfield(p,'expt'); % for structcut to work

% Load up composite maps
filename = [compositedir '/' compositefile];
load(filename);

% Choose the right detectors
switch perwhat
  case 'perfreq'
    bands = unique(p.band);
    bands = bands(find(bands)); % Sometimes 0 is in bands...
    bands = bands(~isnan(bands)); % and Nans too...
    
    for ii = 1:length(bands)
      % Choose relevant dets
      keepme = p.band == bands(ii);
      p0 = structcut(p,keepme);
      ind0 = make_ind(p0);
      map0 = map(keepme);

      % Make the stacked map
      stack = stackmaps(map0,ad,method,weight,deoutlier);
      if plotbeam
        plot_a_beam(stack,ad,expt,year,norm,...
                    [num2str(bands(ii)) ' GHz']);
      end

      % Make the radial profile
      [r,profile] = get_beamprofile(ad,stack,'median',norm);
      if plotprofile
        plot_a_profile(r,profile,expt,year,norm,...
                       [num2str(bands(ii)) ' GHz']);
      end
      
      % Save the stacked map/profile
      savefile = [stackmapdir '/' stackmapfile ...
                  num2str(bands(ii)) '.mat'];
      if exist(savefile,'file')
        s = input([savefile ' exists.  Overwrite?'],'s');
        if strcmp(s,'n') || strcmp(s,'no')
          return
        end
      end
      save(savefile,'stack','ad','r','profile');
      
      % Make/save the B_l 
      savefile = [bldir '/' blfile num2str(bands(ii)) '.fits'];
      [l,B_l] = make_save_bl(ad,stack,expt,...
                             bands(ii),chopper,savefile);
      if plotbl
        plot_a_bl(l,B_l,expt,year,...
                  [num2str(bands(ii)) ' GHz']);
      end
    end
    
  case 'perrx'
    rxNum = unique(p.rx);
    
    for ii = 1:length(rxNum)
      % Choose relevant dets
      keepme = p.rx == rxNum(ii);
      p0 = structcut(p,keepme);
      ind0 = make_ind(p0);
      map0 = map(keepme);
      
      % Make the stacked map
      stack = stackmaps(map0,ad,method,weight,deoutlier);
      if plotbeam
        plot_a_beam(stack,ad,expt,year,norm,...
                    ['rx' num2str(rxNum(ii))]);
      end
      
      % Make the radial profile
      [r,profile] = get_beamprofile(ad,stack,'median',norm);
      if plotprofile
        plot_a_profile(r,profile,expt,year,norm,...
                       ['rx' num2str(rxNum(ii))]);
      end
      
      % Save the stacked map/profile
      savefile = [stackmapdir '/' stackmapfile ...
                  'rx' num2str(rxNum(ii)) '.mat'];
      if exist(savefile,'file')
        s = input([savefile ' exists.  Overwrite?'],'s');
        if strcmp(s,'n') || strcmp(s,'no')
          return
        end
      end
      save(savefile,'stack','ad','r','profile');
      
      % Make/save the B_l
      savefile = [bldir '/' blfile 'rx' num2str(rxNum(ii)) '.fits'];
      [l,B_l] = make_save_bl(ad,stack,expt,...
                             nanmedian(p0.band),chopper,savefile);
      if plotbl
        plot_a_bl(l,B_l,expt,year,...
                  ['rx' num2str(rxNum(ii))]);
      end
    end
    
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function stack = stackmaps(map,ad,method,weight,deoutlier)

% Make a map array which we can easily manipulate
mapstack = NaN(ad.N_pix(1),ad.N_pix(2),length(map));

for ii = 1:length(map)
  if ~isempty(map(ii).T)
    mapstack(:,:,ii) = map(ii).T;
    % Peak normalize each map - should be okay since these are already
    % composites and should be largely artifact-free
    mapstack(:,:,ii) = mapstack(:,:,ii)./nanmax(map(ii).T(:));
  end
end

% Stack according to method
switch method
  case 'median'
    stack = nanmedian(mapstack,3);
  case 'mean'
    % Weight each map
    switch weight
      case 'equal'
	stack = nanmean(mapstack,3);
      case 'std'
      case 'cmb'
    end
end

if deoutlier
  % Occasionally, the stacked map will contain outliers which are
  % likely due to the masking (they tend to look like spokes).  These
  % are definitely not part of the beam and could probably be removed
  % with more careful component mapmaking - however, it's difficult to
  % figure out a priori how to do this.  Here if requested, find the
  % median and std in each annular bin of the map and replace any
  % values more than N sigma away from the median with a value
  % consistent with the distribution of the rest of the annulus.
  dr = 0.1; % degree annuli
  rad = ad.t_r * 180/pi; % easy-to-use radial bins
  rbins = 0:dr:max(rad(:));
  
  for ii = 1:length(rbins)-1
    inds = find( rad > rbins(ii) & rad <= rbins(ii+1) );
    if ~isempty(inds)
      vals = stack(inds);
      med = nanmedian(vals);
      std = nanstd(vals);
      outliers = find( abs(vals - med) > 3*std );
      if ~isempty(outliers)
        vals(outliers) = randn(size(outliers))*0.5*std + med;
        stack(inds) = vals;
      end
      clear vals med std outliers
    end
    clear inds
  end
  
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [l,B_l] = make_save_bl(ad,map,expt,band,chopper,savefile)

% Take basic B_l
[l,B_l] = b2bl(ad,map,0);

% Correct for chopper
diam = chopper * 0.0254; % inches to meters
switch expt
  case {'b2','bicep3'}
    dist = 195;
  case 'keck'
    dist = 211; % distance from MAPO to DSL, meters
end

radius = atan((diam/2) / dist); % radius in radians
arg = l * radius;
bl_anal = 2*besselj(1,arg)./arg;
% For some reason the first element is NaN...
bl_anal(1) = 1;
% Don't apply correction if chopper = 0 (get NaNs)
if chopper > 0
  B_l = B_l./bl_anal;
end

% Cut off at the MTF
% See 20151116_k2015bl posting for brief summary
switch expt
  case 'keck'
    switch band
      case 100
        l_c = 590; % 534;  
      case 150
        l_c = 933; % 830;  
      case 210 
	l_c = 1399; % 1078;
      case 220
        l_c = 1438; % 1139; -- used in previous analysis
      case 270
	l_c = 1714; % 1272;
    end
  case 'b2'
    l_c = 933; % 830;
  case 'bicep3'
    l_c = 1170; % 900;
end

B_l(l > l_c) = 0;

% Write to fits file
if savefile
  if exist(savefile,'file')
    s = input([savefile ' exists.  Overwrite?'],'s');
    if strcmp(s,'n') || strcmp(s,'no')
      return
    end
  end
  %write_beam_file(savefile,l,B_l);
  % Deprecated subfunction, use pipeline standard instead
  write_fits_bls(savefile,l,B_l); 
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function write_beam_file(savefile,l,Bl)

dl = l(2:end) - l(1:end-1);
if any(find(dl ~= 1))
  error('B_l must be monotonic in ell and delta_ell must be one');
end

% Load a sample beam window funciton
filename = '/n/home03/rwa/wmap7yr/bl_v_wmap_7yr.fits';
data = fitsread(filename,'BinTable');
info = fitsinfo(filename);

% Write file
data_out = {cvec(Bl)};

info.BinaryTable.FieldFormat{1} = 'D';
info.BinaryTable.Keywords{10,2} = '1D';

info.BinaryTable.Rows = numel(l);
info.BinaryTable.DataSize = 8*numel(l);
info.BinaryTable.Keywords{5,2} = numel(l);

fitswrite_bintable(data_out,info,savefile)

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plot_a_beam(map,ad,expt,year,norm,descr)

figure
clf
setwinsize(gcf,800,800)

% This assumes the map is in x'/y' as documented in ffbm_rotatemaps -
% i.e. dimension 1 = increasing x prime
%      dimension 2 = decreasing y prime
%      ad.t_val_deg{1} = increasing x prime
%      ad.t_val_deg{2} = decreasing y prime 

switch norm
  case 'peak'
    imagesc(ad.t_val_deg{2},ad.t_val_deg{1},...
            10*log10(abs(map./nanmax(map(:)))));
    cb = colorbar();
    caxis([-60 0])
    xlabel(cb,'dB');
  case 'dBi'
    % Total map power is 1
    map = map / nansum(map(:));
    isomap = ones(size(map));
    isomap = isomap / nansum(isomap(:));
    % Spread out the isotropic radiator over the whole sky
    area = ad.Field_size_deg(1)*ad.Field_size_deg(2);
    isomap = isomap * (area / (4*pi*(180/pi)^2));
    imagesc(ad.t_val_deg{2},ad.t_val_deg{1},...
            10*log10(abs(map./isomap)));
    cb = colorbar();
    caxis([0 60]);
    xlabel(cb,'dBi');
end

set(gca,'YDir','normal')
set(gca,'XDir','reverse')
xlabel('y'', Degrees')
ylabel('x'', Degrees')
colormap('gray')
axis equal

axis([-8 8 -8 8])

switch expt
  case 'keck'
    titlestr = ['Keck ' num2str(year) ' ' descr]; 
    title(titlestr);
    savename = [titlestr '_avgbeam_' norm];
    savename = strrep(savename,' ','_');
    print('-depsc2',savename);
  case 'bicep3'
    titlestr = ['BICEP3 ' num2str(year)];
    title(titlestr);
    savename = [titlestr '_avgbeam_' norm];
    savename = strrep(savename,' ','_');
    print('-depsc2',savename);
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plot_a_profile(r,profile,expt,year,norm,descr)

figure
clf
setwinsize(gcf,400,200)

switch norm
  case 'peak'
    plot(r,10*log10(abs(profile./nanmax(profile))),'LineWidth',2);
    ylim([-60 0]);
    set(gca,'YTick',[-60 -50 -40 -30 -20 -10 0]);
    grid on
    ylabel('dB')
  case 'dBi'
    % This assumes you asked for 'dBi' out of get_beamprofile, so you
    % don't have to manipulate at all
    plot(r,10*log10(profile),'LineWidth',2);
    ylim([-10 60]);
    set(gca,'YTick',[-10 0 10 20 30 40 50 60]);
    grid on
    ylabel('dBi');
end

xlabel('Degrees from Beam Center')
xlim([0 8]); % This is to match the beams paper plot

switch expt
  case 'keck'
    titlestr = ['Keck ' num2str(year) ' ' descr];
    title(titlestr);
    savename = [titlestr '_profile_' norm];
    savename = strrep(savename,' ','_');
    print('-depsc2',savename);
  case 'bicep3'
    titlestr = ['BICEP3 ' num2str(year)];
    title(titlestr);
    savename = [titlestr '_profile_' norm];
    savename = strrep(savename,' ','_');
    print('-depsc2',savename);
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plot_a_bl(l,B_l,expt,year,descr)
% Do both linear/log plots

figure
clf
setwinsize(gcf,400,300)

plot(l,B_l,'LineWidth',2);
xlabel('ell')
ylabel('B_l')
set(gca,'YTick',[0 0.25 0.5 0.75 1]);
grid on
xlim([0 1600])
ylim([0 1])

switch expt
  case 'keck'
    titlestr = ['Keck ' num2str(year) ' ' descr];
    title(titlestr);
    savename = [titlestr '_Bl_lin'];
    savename = strrep(savename,' ','_');
    print('-depsc2',savename);
  case 'bicep3'
    titlestr = ['BICEP3 ' num2str(year)];
    title(titlestr);
    savename = [titlestr '_Bl_lin'];
    savename = strrep(savename,' ','_');
    print('-depsc2',savename);
end

figure
clf
setwinsize(gcf,400,300)

semilogy(l,B_l,'LineWidth',2);
xlabel('ell')
ylabel('B_l')
grid on
xlim([0 1600])
ylim([1e-4 1])

switch expt
  case 'keck'
    titlestr = ['Keck ' num2str(year) ' ' descr];
    title(titlestr);
    savename = [titlestr '_Bl_log'];
    savename = strrep(savename,' ','_');
    print('-depsc2',savename);
  case 'bicep3'
    titlestr = ['BICEP3 ' num2str(year)];
    title(titlestr);
    savename = [titlestr '_Bl_log'];
    savename = strrep(savename,' ','_');
    print('-depsc2',savename);
end

return
