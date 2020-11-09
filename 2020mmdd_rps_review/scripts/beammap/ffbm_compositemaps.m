function ffbm_compositemaps(compopt)
% function ffbm_compositemaps(compopt)
%
% Take rotated/masked beam maps (from ffbm_makesingledatafile) and composite 
% them together.  Should apply cuts generated in ffbm_makecuts.
% Makes either standard composites (using all good maps) or split maps:
% half 1, half 2, half sum, and half diff (from two distinct halves)
% Makes both median and mean maps
%
% After 2018 BM pipeline upgrade, the masking (and rotating in the
% case of dk0 maps) is done in ffbm_makesingledatafile, so we load
% those data products for compositing.
%
% INPUTS (should all be sent in with compopt)
%
%   expt:         'bicep2','keck','bicep3'
%   year:         bicep2: 2012
%                 keck: 2012,...,2018
%                 bicep3: 2016,...,2018
%
% OPTIONAL INPUTS
%
%   mapsize:       map size in deg (8 default), used for map filename
%                  This refers to the radius from beam center - masked
%                  maps are actually square, so in this function we cut
%                  down to a circle with this radius
%   coord:         Coordinate system (used for filename)
%                  'dk0' (default) is the standard input to reduc_makesim
%                  'xpyp' for normal visualization (x'/y')
%   xyear:         Flag to use 'xyear' rotated maps (0 default)
%                  If 1, turns off applycuts (since this was already
%                  applied in the xyear rotated maps)
%   componentdir:  directory from which to load component maps
%                  (default beammaps/maps_component)
%   componentfile:   file in rotateddir to load
%                  (default ffbm_year_allcomp_size_masked_coord)
%   compositedir:  directory in which to save composite maps
%                  (default beammaps/maps_composite)
%   compositefile: file in compositedir to save to
%                  (default ffbm_year_all)
%   applycuts:     0: use all maps!  (please don't)
%                  1: use only maps which pass cuts (default)
%   cutdir:        directory from which to load the cut structure
%                  (default beammaps/cuts)
%   cutfile:       file in cutdir to load
%                  (default ffbm_cuts_year)
%   cutlist:       criteria on which to cut 
%                  (default everything in cutfile)
%   onlycommon_ab: 0: use all schedules passing cuts (default)
%                  1: only keep schedules where A/B were both measured
%   normmethod:    'integral' (default), 'peak' (needs non-NaN fit)
%   weight:        'equal' (default),'std'
%                  Factor by which to multiply each individual map
%   prepforsim:    0: just make 'composite' struct for plotting
%                  1: also save 'map','ad' vars for reduc_makesim
%                     Sets all NaNs to zero, zeroes outsize a radius
%                     Only saves pairs which have at least one good A and B
%                     map  (apply channel cuts later)
%   normforsim:    Which composites to use for sim preparation
%                  'median' (default) takes median through all runs per det
%                  'mean' takes equal-weighted mean
%   radtozero:     radius outside of which to zero (default is size, above)
%   fillnans:      how to fill in NaNs within r < radtozero
%                  'zero' (default) just replaces all NaNs with zero
%                  'noise' inserts noise consistent with that radius
%   makesplits:    0: just make normal composites (default)
%                  1: make two equal-weight splits, sum and diff
%   splitdir:      directory in which to save split maps (for makesplits=1)
%                  (default beammaps/maps_split)
%   splitfile:     file in splitdir to save to (adds 1/2 to end of name)
%                  (default ffbm_year_half[1/2/sum/diff])
%   suffix:        optional string to add to map name, '' default

% Parse compopt
expt = compopt.expt;
year = compopt.year;
if ~isfield(compopt,'mapsize');
  sizestr = '8deg';
  mapsize = 8;
  compopt.mapsize = mapsize;
else
  sizestr = [strrep(num2str(compopt.mapsize),'.','p') 'deg'];
  mapsize = compopt.mapsize;
end
if ~isfield(compopt,'coord')
  compopt.coord = 'dk0';
  coord = compopt.coord;
else
  coord = compopt.coord;
end
if ~isfield(compopt,'xyear')
  compopt.xyear = 0;
  xyear = compopt.xyear;
else
  xyear = compopt.xyear;
end
if compopt.xyear
  xyearstr = '_xyear';
else
  xyearstr = '';
end
if ~isfield(compopt,'suffix')
  suffix = '';
else
  suffix = compopt.suffix;
end
if ~isfield(compopt,'componentdir')
  compopt.componentdir = 'beammaps/maps_component';
  componentdir = compopt.componentdir;
else
  componentdir = compopt.componentdir;
end
if ~isfield(compopt,'componentfile')
  compopt.componentfile = ['ffbm_' num2str(year) '_allcomp_' ...
                      sizestr '_masked_' coord suffix];
  componentfile = compopt.componentfile;
else
  componentfile = compopt.componentfile;
end
if ~isfield(compopt,'applycuts')
  compopt.applycuts = 1;
  applycuts = compopt.applycuts;
else
  applycuts = compopt.applycuts;
end
% Turn off cuts manually if xyear
if compopt.xyear
  applycuts = 0;
end
if ~isfield(compopt,'cutdir')
  compopt.cutdir = 'beammaps/cuts';
  cutdir = compopt.cutdir;
else
  cutdir = compopt.cutdir;
end
if ~isfield(compopt,'cutfile')
  compopt.cutfile = ['ffbm_cuts_' num2str(year)];
  cutfile = compopt.cutfile;
else
  cutfile = compopt.cutfile;
end
if ~isfield(compopt,'cutlist')
  compopt.cutlist = {'mirror','notlight','nanfit','peak',...
      'sigma','ellip','hand'};
  cutlist = compopt.cutlist;
else
  cutlist = compopt.cutlist;
end
if ~isfield(compopt,'compositedir')
  compopt.compositedir = 'beammaps/maps_composite';
  compositedir = compopt.compositedir;
else
  compositedir = compopt.compositedir;
end
if ~isfield(compopt,'compositefile')
  compopt.compositefile = ['ffbm_' num2str(year) '_all'...
	xyearstr '_' sizestr '_' coord suffix];
  compositefile = compopt.compositefile;
else
  compositefile = compopt.compositefile;
end
if ~isfield(compopt,'onlycommon_ab')
  compopt.onlycommon_ab = 0;
  onlycommon_ab = compopt.onlycommon_ab;
else
  onlycommon_ab = compopt.onlycommon_ab;
end
if ~isfield(compopt,'normmethod')
  compopt.normmethod = 'integral';
  normmethod = compopt.normmethod;
else
  normmethod = compopt.normmethod;
end
if ~isfield(compopt,'weight')
  compopt.weight = 'equal';
  weight = compopt.weight;
else
  weight = compopt.weight;
end
if ~isfield(compopt,'prepforsim')
  compopt.prepforsim = 1;
  prepforsim = compopt.prepforsim;
else
  prepforsim = compopt.prepforsim;
end
if ~isfield(compopt,'normforsim')
  compopt.normforsim = 'median';
  normforsim = compopt.normforsim;
else
  normforsim = compopt.normforsim;
end
if ~isfield(compopt,'radtozero')
  compopt.radtozero = compopt.mapsize;
  radtozero = compopt.radtozero;
else
  radtozero = compopt.radtozero;
end
if ~isfield(compopt,'fillnans')
  compopt.fillnans = 'zero';
  fillnans = compopt.fillnans;
else
  fillnans = compopt.fillnans;
end
if ~isfield(compopt,'makesplits')
  compopt.makesplits = 0;
  makesplits = compopt.makesplits;
else
  makesplits = compopt.makesplits;
end
if makesplits
  if ~isfield(compopt,'splitdir')
    compopt.splitdir = 'beammaps/maps_split';
    splitdir = compopt.splitdir;
  else
    splitdir = compopt.splitdir;
  end
  if ~isfield(compopt,'splitfile')
    compopt.splitfile = ['ffbm_' num2str(year) '_half'...
	  xyearstr '_' sizestr '_' coord suffix '_'];
    splitfile = compopt.splitfile;
  else
    splitfile = compopt.splitfile;
  end
end
  
[p ind] = get_array_info([num2str(year) '0201']);

% Load up component maps
filename = [componentdir '/' componentfile];
disp(['Loading component maps: ' filename]);
load(filename);

% Load cuts
if applycuts
  cutname = [cutdir '/' cutfile];
  load(cutname);
  cut = cuts.cuts;
end

% Initialize map cells and structs to proper length
map_med_1 = cell(1,length(p.gcp));
map_mean_1 = cell(1,length(p.gcp));
map_std_1 = cell(1,length(p.gcp));
map_fit_1 = cell(1,length(p.gcp));

map_med_2 = cell(1,length(p.gcp));
map_mean_2 = cell(1,length(p.gcp));
map_std_2 = cell(1,length(p.gcp));
map_fit_2 = cell(1,length(p.gcp));

map_1(length(p.gcp)).T = [];
map_2(length(p.gcp)).T = [];

for ii = 1:length(ind.la)
  
  if mod(ii,10) == 0
    disp(['Compositing pair ' num2str(ii) ' / ' ...
	  num2str(length(ind.la))]);
  end
  % Make a total cut mask
  totmask_a = ones(length(comp.bm.number),1);
  totmask_b = ones(length(comp.bm.number),1);
  
  if applycuts
    for jj = 1:length(comp.bm.number);
      a = structcut(cut(jj),ind.la(ii));
      b = structcut(cut(jj),ind.lb(ii));
      totcut_a = 0;
      totcut_b = 0;
      for kk = 1:length(cutlist)
	totcut_a = totcut_a + getfield(a,cutlist{kk});
	totcut_b = totcut_b + getfield(b,cutlist{kk});
      end
      % Handcut == 1 means the fit may suck, but it's probably okay for
      % compositing
      if totcut_a == 0 | a.hand == 1 
	totmask_a(jj) = 0;
      end
      if totcut_b == 0 | b.hand == 1
	totmask_b(jj) = 0;
      end
    end
  else
    totmask_a = zeros(length(comp.bm.number),1);
    totmask_b = zeros(length(comp.bm.number),1);
  end
  
  % Indices of good maps
  goodinds_a = find(totmask_a == 0);
  goodinds_b = find(totmask_b == 0);
  
  % If we want to keep only schedules where A/B were both measured
  if onlycommon_ab
    goodinds_a = intersect(goodinds_a,goodinds_b);
    goodinds_b = goodinds_a;
  end

  goodinds_a_1 = goodinds_a;
  goodinds_b_1 = goodinds_b;

  % Choose indices for half1/2 splits if required (overwrite _1)
  if makesplits
    % Require an even number of maps (otherwise we risk split leakage)
    % If odd, discard last map
    if mod(length(goodinds_a),2)
      goodinds_a = goodinds_a(1:end-1);
    end
    if mod(length(goodinds_b),2)
      goodinds_b = goodinds_b(1:end-1);
    end
    % I used to do odd indices for map1 and even for map2.  However, if
    % you have a very even beam mapping cadence and there are repeatable,
    % systematic differences between (say) dk angles, then these could
    % show up in the half1 - half2 difference when you don't want them
    % to.  So instead, use a random permutation of good indices.  No
    % reason to have the same indices for A and B...
    perm_a = randperm(length(goodinds_a));
    perm_b = randperm(length(goodinds_b));
    % NOW choose even/odd indices of the permutation for 1 and 2
    % Odd indices for 1
    goodinds_a_1 = perm_a(1:2:end); %goodinds_a(1:2:end);
    goodinds_b_1 = perm_b(1:2:end); %goodinds_b(1:2:end);
    % Even indices for 2
    goodinds_a_2 = perm_a(2:2:end); %goodinds_a(2:2:end);
    goodinds_b_2 = perm_b(2:2:end); %goodinds_b(2:2:end);
  end
      
  % Final gridding 
  map_sum_a_1 = NaN([length(comp.ad.t_val_deg{1}) ...
	length(comp.ad.t_val_deg{2}) length(goodinds_a_1)]);
  map_sum_b_1 = NaN([length(comp.ad.t_val_deg{1}) ...
	length(comp.ad.t_val_deg{2}) length(goodinds_b_1)]);
  if makesplits
    map_sum_a_2 = NaN([length(comp.ad.t_val_deg{1}) ...
	  length(comp.ad.t_val_deg{2}) length(goodinds_a_2)]);
    map_sum_b_2 = NaN([length(comp.ad.t_val_deg{1}) ...
	  length(comp.ad.t_val_deg{2}) length(goodinds_b_2)]);
  end
  
  % Accumulate good component maps
  for jj = 1:length(goodinds_a_1)
    fit = comp.map{ind.la(ii)}.fit{goodinds_a_1(jj)};
    % Normalize each beam - if integral, use only inner 2x2 deg
    map_final = norm_beam(comp.map{ind.la(ii)}.component{goodinds_a_1(jj)},...
                         normmethod,fit,2,comp.ad);
    switch weight
      case 'equal'
	map_sum_a_1(:,:,jj) = map_final;
      case 'std'
	weightmap = 1./nanstd(map_final(:));
	map_sum_a_1(:,:,jj) = map_final * weightmap;
    end
  end % Component map loop
  for jj = 1:length(goodinds_b_1)
    fit = comp.map{ind.lb(ii)}.fit{goodinds_b_1(jj)};
    % Normalize each beam - if integral, use only inner 2x2 deg
    map_final = norm_beam(comp.map{ind.lb(ii)}.component{goodinds_b_1(jj)},...
                          normmethod,fit,2,comp.ad);
    switch weight
      case 'equal'
	map_sum_b_1(:,:,jj) = map_final;
      case 'std'
	weightmap = 1./nanstd(map_final(:));
	map_sum_b_1(:,:,jj) = map_final * weightmap;
    end
  end % Component map loop
  
  if makesplits
    for jj = 1:length(goodinds_a_2)
      fit = comp.map{ind.la(ii)}.fit{goodinds_a_2(jj)};
      % Normalize each beam - if integral, use only inner 2x2 deg
      map_final = norm_beam(comp.map{ind.la(ii)}.component{goodinds_a_2(jj)},...
                            normmethod,fit);
      switch weight
	case 'equal'
	  map_sum_a_2(:,:,jj) = map_final;
	case 'std'
	  weightmap = 1./nanstd(map_final(:));
	  map_sum_a_2(:,:,jj) = map_final * weightmap;
      end
    end % Component map loop
    for jj = 1:length(goodinds_b_2)
      fit = comp.map{ind.lb(ii)}.fit{goodinds_b_2(jj)};
      % Normalize each beam - if integral, use only inner 2x2 deg
      map_final = norm_beam(comp.map{ind.lb(ii)}.component{goodinds_b_2(jj)},...
                            normmethod,fit);
      switch weight
	case 'equal'
	  map_sum_b_2(:,:,jj) = map_final;
	case 'std'
	  weightmap = 1./nanstd(map_final(:));
	  map_sum_b_2(:,:,jj) = map_final * weightmap;
      end
    end % Component map loop
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Make the composites

  map_med_1{ind.la(ii)} = nanmedian(map_sum_a_1,3);
  map_mean_1{ind.la(ii)} = nanmean(map_sum_a_1,3);
  map_std_1{ind.la(ii)} = nanstd(map_sum_a_1,[],3);
  map_med_1{ind.lb(ii)} = nanmedian(map_sum_b_1,3);
  map_mean_1{ind.lb(ii)} = nanmean(map_sum_b_1,3);
  map_std_1{ind.lb(ii)} = nanstd(map_sum_b_1,[],3);

  if makesplits
    map_med_2{ind.la(ii)} = nanmedian(map_sum_a_2,3);
    map_mean_2{ind.la(ii)} = nanmean(map_sum_a_2,3);
    map_std_2{ind.la(ii)} = nanstd(map_sum_a_2,[],3);
    map_med_2{ind.lb(ii)} = nanmedian(map_sum_b_2,3);
    map_mean_2{ind.lb(ii)} = nanmean(map_sum_b_2,3);
    map_std_2{ind.lb(ii)} = nanstd(map_sum_b_2,[],3);
  end
  
  % From above potential maps, choose which to normalize, fit, and send
  % into sims 'T' struct
  switch normforsim
    case 'median'
      maptofit_a_1 = map_med_1{ind.la(ii)};
      maptofit_b_1 = map_med_1{ind.lb(ii)};
      if makesplits
	maptofit_a_2 = map_med_2{ind.la(ii)};
	maptofit_b_2 = map_med_2{ind.lb(ii)};
      end
    case 'mean'
      maptofit_a_1 = map_mean_1{ind.la(ii)};
      maptofit_b_1 = map_mean_1{ind.lb(ii)};
      if makesplits
	maptofit_a_2 = map_mean_2{ind.la(ii)};
	maptofit_b_2 = map_mean_2{ind.lb(ii)};
      end
  end
  
  % reduc_makesim wants a 'map' struct, there should be no Nans, 
  % and it should be integral normalized
  % 
  % First, set everything outside radtozero = 0
  % Next, integral normalize (but only using stuff within 2 degrees)
  % Normalized map can still have NaNs in it
  if prepforsim
    maptofit_a_1 = zeroradius(maptofit_a_1,comp.ad,radtozero);
    maptofit_a_1 = norm_beam(maptofit_a_1,'integral',[],2,comp.ad);
    maptofit_b_1 = zeroradius(maptofit_b_1,comp.ad,radtozero);
    maptofit_b_1 = norm_beam(maptofit_b_1,'integral',[],2,comp.ad);

    % Now choose whether to insert noise consistent with the beam maps,
    % or to set NaNs to zero
    if ~isempty(goodinds_a_1) & ~isempty(goodinds_b_1)
      maptofit_a_1 = denanify(maptofit_a_1,fillnans);
      maptofit_b_1 = denanify(maptofit_b_1,fillnans);
      map_1(ind.la(ii)).T = maptofit_a_1;
      map_1(ind.lb(ii)).T = maptofit_b_1;
    else
      map_1(ind.la(ii)).T = NaN(size(maptofit_a_1));
      map_1(ind.lb(ii)).T = NaN(size(maptofit_b_1));
    end
  end

  % Now fit - only do the inner 2 degrees!
  map_fit_1{ind.la(ii)} = fit_inner(maptofit_a_1,comp.ad,2);
  map_fit_1{ind.lb(ii)} = fit_inner(maptofit_b_1,comp.ad,2);

  if makesplits
    if prepforsim
      maptofit_a_2 = zeroradius(maptofit_a_2,comp.ad,radtozero);
      maptofit_a_2 = norm_beam(maptofit_a_2,'integral',[],2,comp.ad);
      maptofit_b_2 = zeroradius(maptofit_b_2,comp.ad,radtozero);
      maptofit_b_2 = norm_beam(maptofit_b_2,'integral',[],2,comp.ad);
      
      % Now choose whether to insert noise consistent with the beam maps,
      % or to set NaNs to zero
      if ~isempty(goodinds_a_2) & ~isempty(goodinds_b_2)
        maptofit_a_2 = denanify(maptofit_a_2,fillnans);
        maptofit_b_2 = denanify(maptofit_b_2,fillnans);
	map_2(ind.la(ii)).T = maptofit_a_2;
	map_2(ind.lb(ii)).T = maptofit_b_2;
      else
	map_2(ind.la(ii)).T = NaN(size(maptofit_a_2));
	map_2(ind.lb(ii)).T = NaN(size(maptofit_b_2));
      end
    end
    map_fit_2{ind.la(ii)} = fit_inner(maptofit_a_2,comp.ad,2);
    map_fit_2{ind.lb(ii)} = fit_inner(maptofit_b_2,comp.ad,2);
  end
  
end % Lights loop

% Prep output struct
composite.bm = comp.bm;
composite.ad = comp.ad;
composite.compopt = compopt;
ad = comp.ad;

composite.map_med = map_med_1;
composite.map_mean = map_mean_1;
composite.map_std = map_std_1;
composite.map_fit = map_fit_1;

map = map_1;

switch makesplits
  case 0
    savename = [compositedir '/' compositefile];    
    if prepforsim
      save(savename,'composite','map','ad');
    else
      save(savename,'composite');
    end
    disp(['Saved composite map file: ' savename])
  case 1
    savename = [splitdir '/' splitfile '1'];
    if prepforsim
      save(savename,'composite','map','ad');
    else
      save(savename,'composite');
    end
    disp(['Saved split map file: ' savename])
    
    composite.map_med = map_med_2;
    composite.map_mean = map_mean_2;
    composite.map_std = map_std_2;
    composite.map_fit = map_fit_2;
    
    map = map_2;
    
    savename = [splitdir '/' splitfile '2'];
    if prepforsim
      save(savename,'composite','map','ad');
    else
      save(savename,'composite');
    end
    disp(['Saved split map file: ' savename])
    
    % Make sum
    for ii = 1:length(map_1)
      composite.map_med{ii} = (map_med_1{ii} + map_med_2{ii})/2;
      composite.map_mean{ii} = (map_mean_1{ii} + map_mean_2{ii})/2;
      map(ii).T = (map_1(ii).T + map_2(ii).T) / 2;
    end

    composite.map_std = [];
    composite.map_fit = [];

    savename = [splitdir '/' splitfile 'sum'];
    if prepforsim
      save(savename,'composite','map','ad');
    else
      save(savename,'composite');
    end
    disp(['Saved split map file: ' savename])
    
    % Make diff
    for ii = 1:length(map_1)
      composite.map_med{ii} = (map_med_1{ii} - map_med_2{ii})/2;
      composite.map_mean{ii} = (map_mean_1{ii} - map_mean_2{ii})/2;
      map(ii).T = (map_1(ii).T - map_2(ii).T) / 2;
    end

    composite.map_std = [];
    composite.map_fit = [];

    savename = [splitdir '/' splitfile 'diff'];
    if prepforsim
      save(savename,'composite','map','ad');
    else
      save(savename,'composite');
    end
    disp(['Saved split map file: ' savename])
end

return % Main function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function zeroed_map = zeroradius(map,ad,rad)
% Zero out all pixels greater than rad

[xx yy] = meshgrid(ad.t_val_deg{1},ad.t_val_deg{2});
rr = sqrt(xx.^2 + yy.^2);
map(rr > rad) = 0;
zeroed_map = map;

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function norm_map = norm_beam(map,normmethod,fit,deg,ad)
% Normalizes beam either by integral or peak, based on fits in the
% component map stage
% If integral normalizing, only use map values within a square 
% deg x deg around the beam, so that far-away artifacts don't matter
% To get the normalization factor, inpaint NaNs  so that there is at
% least a semblance of consistency - but keep NaN pixels in the final map

switch normmethod
  case 'integral'
    if ~exist('deg','var') % No 'deg' argument, just use whole map
      norm_map = map./nansum(map(:));
    else
      ax1ind = find(abs(ad.t_val_deg{1}) < deg);
      ax2ind = find(abs(ad.t_val_deg{2}) < deg);
      mapcut = map(ax1ind,ax2ind);
      mapcut = inpaint_nans(mapcut);
      % Presumably mapcut, which has no NaNs, is a decent estimate of the
      % total power in the beam
      norm_map = map./sum(mapcut(:));
      % However, map is allowed to nave NaNs in it
    end
  case 'peak'
    norm_map = map./fit(1);
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function fit = fit_inner(map,ad,deg)
% Fits inner portion of a map - just does square window
% Inpaint nans so the fitter doesn't fail

ax1ind = find(abs(ad.t_val_deg{1}) < deg);
ax2ind = find(abs(ad.t_val_deg{2}) < deg);
mapcut = map(ax1ind,ax2ind);
mapcut = inpaint_nans(mapcut);

% Use optimset to increase tolerance during fitting
opts = optimset('TolFun',1e-8,'TolX',1e-8,'MaxFunEvals',7000,'MaxIter',7000);

fit = normfit2d(ad.t_val_deg{1}(ax1ind),...
                ad.t_val_deg{2}(ax2ind),mapcut,'opts',opts);

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function denaned_map = denanify(map,fillnans);
% Fills NaNs in map with zeros or noise

switch fillnans
  case 'zero'
    inds = find(isnan(map));
    map(inds) = 0;
  case 'noise'
end

denaned_map = map;

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      year:             2012,2013,2014
%      rxNum:            keck only
%      dk:               dkangles to include when adding
%                        B2: [0 90 -90 -180] 
%                        K,rx1: [22 -50 22 -86 -14] PHASE OUT!!! 
%      sc:               schedules to include  PHASE OUT!!!
%                        B2: [26 27 28] 
%                        Keck: [0 1 2] 
%      filename:         filename to save as
%      crop:             0 (do not crop, super large map)
%                        1 (crop to 8x8 deg)
%      savecomp:         0 (default)
%                        1 (save the component maps in a stack)
%      cutmaps:          0 (default, no hand cuts)
%                        1 (hand cuts)
%      plotcomp:         0 (default, no plotting component maps)
%                        1 (plot component maps)
%      subdir:           use with plotcomp
%      choplet:          0 (default)
%                        1 (look for choplet demodded maps)
%      cent:             1 (default, center on common AB centroid)
%                        0 (per-detector centering)
%
% % B2 example
% plotcomp = 1; mapsourcetype = 'uberchopper'; dk = [0 90 -90 -180];
% sc = [26 27 28];
% subdir = ['beammapplots/' mapsourcetype '_' num2str(sc) '_cut']
% subdir = 'beammapplots/uberchopper_all';
% cutmaps = 1; plotcomp = 1;
% compositemaps(experiment,mapsourcetype,year,rxNum,dk,sc,...
%               cutmaps,plotcomp,subdir)
%
% % Keck Rx1
% experiment = 'keck'; mapsourcetype = ''; year = 2012; rxNum = 1;
% dk = [22 -50 22 -86 -14]; sc = [1]; subdir = ''; cutmaps = 0;
% plotcomp = 1; crop = 0; savecomp = 0;
% compositemaps(experiment,mapsourcetype,year,rxNum,dk,sc,...
%               filename,crop,savecomp,cutmaps,plotcomp,subdir)


%{
[time dkangle sched] = ffbm_findbmruntimes(experiment,mapsourcetype,year,rxNum);

switch experiment
  case 'bicep2'
    [pp ind] = get_array_info();    
  case 'keck'
    switch year
      case 2012
	[pp ind] = get_array_info('20120201');	
      case 2013
	[pp ind] = get_array_info('20130201');	
      case 2014
	[pp ind] = get_array_info('20140201');
    end
    cutind = (pp.rx == rxNum);
    pp = structcut(pp,cutind);
    ind = make_ind(pp);    
end

% Load in maps
disp('Loading maps');
kk = 1;

for hh = 1:length(sc)
  for jj = 1:length(dk)
    whichtime = find(dkangle == dk(jj) & sched == sc(hh));
    for ii = 1:length(whichtime)
      switch experiment
	case 'bicep2'
	  ww{kk} = load(['beammaps/maps_component/' char(time{whichtime(ii)}) ...
		'_b2_rot' num2str(dk(jj)) '.mat']);
	case 'keck'
	  if choplet
	    ww{kk} = load(['beammaps/maps_component/' char(time{whichtime(ii)}) ...
		  '_' experiment '_cl_rx' num2str(rxNum) '_rot' ...
		  num2str(dk(jj)) '.mat']); 
	  else
	    try
	      ww{kk} = load(['beammaps/maps_component/' char(time{whichtime(ii)}) ...
		  '_' experiment '_rx' num2str(rxNum) '_rot' num2str(dk(jj)) ...
		  '.mat']); 
	      disp(['Loaded map ' char(time{whichtime(ii)}) '_rx' num2str(rxNum)])
	    catch err
	      disp('no map')
	    end
	    ww{kk}.time = time{whichtime(ii)};
	end
      end
      ww{kk}.dk = dk(jj);
      ww{kk}.sc = sc(hh);
      kk = kk + 1;
    end
  end
end
    
% Final map binning
if crop
  x_bin = -4:0.1:4;
  y_bin = -4:0.1:4;
else
  x_bin = -20:0.1:20;
  y_bin = -20:0.1:20;
end

y_bin = fliplr(y_bin);
[X,Y] = meshgrid(x_bin,y_bin);

% Remove from ww struct if the map didn't load
for ii = 1:length(ww)
  maphere(ii) = isfield(ww{ii},'map');
end
ww = ww(maphere);

disp('Recentering')
for ii = 1:length(ww)
  
  map_cent = NaN(length(y_bin),length(x_bin),528);
  map_w = NaN(length(y_bin),length(x_bin),528);
  Afit_cent = ww{ii}.Afit;
  bmcut = ones(528,1);
  
  % Find common center of A and B
  if cent % Per-pair centroiding -- same for A and B
    xcent_a = (ww{ii}.Afit(2,ind.a) + ww{ii}.Afit(2,ind.b))./2;
    xcent_b = (ww{ii}.Afit(2,ind.a) + ww{ii}.Afit(2,ind.b))./2;
    ycent_a = (ww{ii}.Afit(3,ind.a) + ww{ii}.Afit(3,ind.b))./2;
    ycent_b = (ww{ii}.Afit(3,ind.a) + ww{ii}.Afit(3,ind.b))./2;
  else    % Per-detector centroiding
    xcent_a = ww{ii}.Afit(2,ind.a);
    xcent_b = ww{ii}.Afit(2,ind.b);
    ycent_a = ww{ii}.Afit(3,ind.a);
    ycent_b = ww{ii}.Afit(3,ind.b);
  end
  
  % Keck-specific thresholds
  switch experiment
    case 'bicep2'
      bada = [];
      badb = [];
    case 'keck'
      % NaN out the pixels without good beams using the peak of the fit
      switch year
	case 2012
	  bada = find(ww{ii}.Afit(1,ind.a) < 10 | ww{ii}.Afit(1,ind.a) > 300);
	  badb = find(ww{ii}.Afit(1,ind.b) < 10 | ww{ii}.Afit(1,ind.b) > 300);
	case 2013
	  switch rxNum
	    case 0
	      bada = find(ww{ii}.Afit(1,ind.a) < 10 | ww{ii}.Afit(1,ind.a) > 400);
	      badb = find(ww{ii}.Afit(1,ind.b) < 10 | ww{ii}.Afit(1,ind.b) > 400);		    
	    case 1
	      bada = find(ww{ii}.Afit(1,ind.a) < 10 | ww{ii}.Afit(1,ind.a) > 500);
	      badb = find(ww{ii}.Afit(1,ind.b) < 10 | ww{ii}.Afit(1,ind.b) > 500);
	    case 2
	      bada = find(ww{ii}.Afit(1,ind.a) < 0 | ww{ii}.Afit(1,ind.a) > 500);
	      badb = find(ww{ii}.Afit(1,ind.b) < 0 | ww{ii}.Afit(1,ind.b) > 500);
	    case 3
	      bada = find(ww{ii}.Afit(1,ind.a) < 10 | ww{ii}.Afit(1,ind.a) > 400);
	      badb = find(ww{ii}.Afit(1,ind.b) < 10 | ww{ii}.Afit(1,ind.b) > 400);	
	    case 4
	      bada = find(ww{ii}.Afit(1,ind.a) < 10 | ww{ii}.Afit(1,ind.a) > 600);
	      badb = find(ww{ii}.Afit(1,ind.b) < 10 | ww{ii}.Afit(1,ind.b) > 600);
	  end % rxNum switch
	case 2014
	  switch rxNum
	    case 0
	      bada = find(ww{ii}.Afit(1,ind.a) < 10 | ww{ii}.Afit(1,ind.a) > 200);
	      badb = find(ww{ii}.Afit(1,ind.b) < 10 | ww{ii}.Afit(1,ind.b) > 200);
	    case 1
	      bada = find(ww{ii}.Afit(1,ind.a) < 10 | ww{ii}.Afit(1,ind.a) > 500);
	      badb = find(ww{ii}.Afit(1,ind.b) < 10 | ww{ii}.Afit(1,ind.b) > 500);
	    case 2
	      bada = find(ww{ii}.Afit(1,ind.a) < 10 | ww{ii}.Afit(1,ind.a) > 200);
	      badb = find(ww{ii}.Afit(1,ind.b) < 10 | ww{ii}.Afit(1,ind.b) > 200);
	    case 3
	      bada = find(ww{ii}.Afit(1,ind.a) < 10 | ww{ii}.Afit(1,ind.a) > 400);
	      badb = find(ww{ii}.Afit(1,ind.b) < 10 | ww{ii}.Afit(1,ind.b) > 400);	
	    case 4
	      bada = find(ww{ii}.Afit(1,ind.a) < 10 | ww{ii}.Afit(1,ind.a) > 600);
	      badb = find(ww{ii}.Afit(1,ind.b) < 10 | ww{ii}.Afit(1,ind.b) > 600);
	  end % rxNum switch
      end % year switch
  end  % experiment switch
  
  % Remove NaN fits
  tmp = find(isnan(ww{ii}.Afit(1,ind.a))); 
  bada = [bada tmp];
  tmp = find(isnan(ww{ii}.Afit(1,ind.b)));
  badb = [badb tmp];
    
  % Remove if the beam center is nowhere close to where it should be
  % Threshold: 2 degrees from median
  xmed = nanmedian(ww{ii}.Afit(2,:));
  ymed = nanmedian(ww{ii}.Afit(3,:));
  tmp = find(abs(ww{ii}.Afit(2,ind.a) - xmed) > 2 | ...
      abs(ww{ii}.Afit(3,ind.a) - ymed) > 2); 
  bada = [bada tmp];
  tmp = find(abs(ww{ii}.Afit(2,ind.b) - xmed) > 2 | ...
      abs(ww{ii}.Afit(3,ind.b) - ymed) > 2); 
  badb = [badb tmp];
  
  % NaN out bad maps 
  % ww{ii}.map(:,:,ind.a(bada))=NaN;
  % ww{ii}.map(:,:,ind.b(badb))=NaN;
  
  % 1 is good, NaN is bad
  bmcut(ind.a(bada)) = NaN;
  bmcut(ind.b(badb)) = NaN;
  
  % Maybe have the end statement here?  Moved for B2
  
  % Center the beams of detectors that don't have working pairs - 
  % i.e. if A is bad and B is good, center both on B.  
  % If both are bad, center on nothing.
  badboth = intersect(bada,badb);
  bada = setdiff(bada,badboth);
  badb = setdiff(badb,badboth);
  
  xcent_a(bada) = ww{ii}.Afit(2,ind.b(bada));
  xcent_b(bada) = ww{ii}.Afit(2,ind.b(bada));
  ycent_a(bada) = ww{ii}.Afit(3,ind.b(bada));
  ycent_b(bada) = ww{ii}.Afit(3,ind.b(bada));
  xcent_a(badb) = ww{ii}.Afit(2,ind.a(badb));
  xcent_b(badb) = ww{ii}.Afit(2,ind.a(badb));
  ycent_a(badb) = ww{ii}.Afit(3,ind.a(badb));
  ycent_b(badb) = ww{ii}.Afit(3,ind.a(badb));
  
  xcent_a(badboth) = NaN;
  xcent_b(badboth) = NaN;
  ycent_a(badboth) = NaN;
  ycent_b(badboth) = NaN;
  
  % Move the maps to a new map with a common center of (0,0) using interp2
  for jj = 1:length(ind.a)
    
    % Separately for A and B
    if ~isnan(xcent_a(jj)) & ~isnan(xcent_b(jj)) % A and B both exist
      xx_a = ww{ii}.x_bin - xcent_a(jj);
      yy_a = ww{ii}.y_bin - ycent_a(jj);
      xx_b = ww{ii}.x_bin - xcent_b(jj);
      yy_b = ww{ii}.y_bin - ycent_b(jj);
    else 
      if ~isnan(ww{ii}.Afit(2,ind.a(jj))) % A fit exists, use A 
	xx_a = ww{ii}.x_bin - ww{ii}.Afit(2,ind.a(jj));
	yy_a = ww{ii}.y_bin - ww{ii}.Afit(3,ind.a(jj));
	xx_b = ww{ii}.x_bin - ww{ii}.Afit(2,ind.a(jj));
	yy_b = ww{ii}.y_bin - ww{ii}.Afit(3,ind.a(jj));
      elseif ~isnan(ww{ii}.Afit(2,ind.b(jj))) % B fit exists, use B
	xx_a = ww{ii}.x_bin - ww{ii}.Afit(2,ind.b(jj));
	yy_a = ww{ii}.y_bin - ww{ii}.Afit(3,ind.b(jj));
	xx_b = ww{ii}.x_bin - ww{ii}.Afit(2,ind.b(jj));
	yy_b = ww{ii}.y_bin - ww{ii}.Afit(3,ind.b(jj));
      else % None exist
	xx_a = ww{ii}.x_bin;
	yy_a = ww{ii}.y_bin;
	xx_b = ww{ii}.x_bin;
	yy_b = ww{ii}.y_bin;
      end
    end
    
    [XX_a,YY_a] = meshgrid(xx_a,yy_a);
    [XX_b,YY_b] = meshgrid(xx_b,yy_b);
    
    map_cent(:,:,ind.a(jj)) = interp2(XX_a,YY_a,ww{ii}.map(:,:,ind.a(jj)),X,Y);
    map_cent(:,:,ind.b(jj)) = interp2(XX_b,YY_b,ww{ii}.map(:,:,ind.b(jj)),X,Y);
    
    map_cent(:,:,ind.a(jj)) = map_cent(:,:,ind.a(jj))/ww{ii}.Afit(1,ind.a(jj));
    map_cent(:,:,ind.b(jj)) = map_cent(:,:,ind.b(jj))/ww{ii}.Afit(1,ind.b(jj));
  end
  
  % Move the Afit values to match recentered map
  Afit_cent(2,ind.a) = ww{ii}.Afit(2,ind.a) - xcent_a;
  Afit_cent(3,ind.a) = ww{ii}.Afit(3,ind.a) - ycent_a;
  Afit_cent(2,ind.b) = ww{ii}.Afit(2,ind.b) - xcent_b;
  Afit_cent(3,ind.b) = ww{ii}.Afit(3,ind.b) - ycent_b;
  
  winaz = (x_bin > -2) & (x_bin <= 2);
  winel = (y_bin > -2) & (y_bin <= +2);
  
  [maskaz maskel] = meshgrid(winaz,winel);
  mask = maskaz & maskel;
  
  % Find weights
  for kk = 1:528
    tmap = map_cent(:,:,kk);
    tmap(mask) = NaN;
    weights(kk) = 1./nanstd(nanstd(tmap));
    tmap1 = map_cent(:,:,kk);
    tmap1(~isnan(tmap1)) = 1;
    % tmap1(isnan(tmap1)) = 0;
    map_w(:,:,kk) = tmap1*weights(kk);
  end
  
  % Add maps with equal weights and count maps that go into 1 pixel
  for kk = 1:528
    tmap1 = map_cent(:,:,kk);
    tmap1(~isnan(tmap1)) = 1;
    map_w(:,:,kk) = tmap1;
  end
  
  % Put everything in a new struct
  qq{ii}.Afit = Afit_cent;
  qq{ii}.map = map_cent;
  qq{ii}.x_bin = x_bin;
  qq{ii}.y_bin = y_bin;
  qq{ii}.weightmask = map_w;
  qq{ii}.dk = ww{ii}.dk;
  qq{ii}.sc = ww{ii}.sc; 
  if isfield(ww{ii},'time')
    qq{ii}.time = ww{ii}.time;
  end
  
  qq{ii}.badmapcut = bmcut;
end % Maps loop

clear ww;

% Hand cuts
if cutmaps
  switch experiment
    case 'bicep2'
      % Made for B2 but not used in end
      mapstocut = makehandcut(sc);
      mapstocut = load(handcut);
    case 'keck'
      switch year
	case 2013
	  load(['~clwong/20140506_keckbm/handcut_rx' num2str(rxNum)]);
	  % ~clwong/20140506_keckbm/makehandcut_keck2013.m makes
	  % handcut_rxnum.mat
	case 2014
	  hc = load(['mask_ffbm_2014/mask_rx' num2str(rxNum)]);
	  % Turn this into CLW's mapstocut
	  mapstocut = zeros(length(qq),528);
	  for ii = 1:length(qq)
	    for jj = 1:length(hc.mask)
	      if strcmp(qq{ii}.time,hc.mask(jj).time)
		mapstocut(ii,:) = hc.mask(jj).handcut;
	      end
	    end
	  end
      end % year switch    
  end % experiment switch
  
  for ii = 1:length(qq)
    indtocut = find(mapstocut(ii,:));
    for jj = 1:length(indtocut)
      qq{ii}.weightmask(:,:,indtocut(jj)) = NaN;
    end
  end
end % Hand cuts

% If Keck, cut out maps in a bad position
if strcmp(experiment,'keck')
  for ii = 1:length(qq)
    switch year
      case 2012
	[mask100 mask150] = ffbm_maskMaker('back');
	mask = mask150;
	mask(isnan(mask)) = 1;
      case 2013
	[mask100 mask150] = ffbm_maskMaker('back');
	mask = mask150;
	mask(isnan(mask)) = 1;
      case 2014
	tmp = qq{ii}.time;
	mirror_change_date = '2014-02-22 02:00:00';
	bm_date = sprintf('%s-%s-%s %s:%s:%s',...
	    tmp(1:4),tmp(5:6),tmp(7:8),tmp(10:11),...
	    tmp(12:13),tmp(14:15));
	if datenum(bm_date) < datenum(mirror_change_date)		
	  [mask100 mask150] = ffbm_maskMaker('back');
	else
	  [mask100 mask150] = ffbm_maskMaker('front');
	end
	switch rxNum
	  case {0,2}
	    mask = mask100;
	  case {1,3,4}
	    mask = mask150;
	end
	mask(isnan(mask)) = 1;
    end
    
    %  if year == 2012 | year == 2013
    %    [mask100 mask150] = ffbm_maskMaker('back');
    %    mask = mask150;
    %    mask(isnan(mask)) = 1; % Doesn't like NaNs below
    %  end
    
    
    qq{ii}.poscut = NaN(1,528);
    % Figure out where the rx is
    rxpos = dk_rx_pos_ffbm(rxNum,qq{ii}.dk);
    % pos3 is the base of drum.  Here 0 is good and nonzero is bad, so
    % reverse it here
    poscut = ~mask(rxpos,:);
    % -> Now 0 is bad, so turn 0s into NaNs
    qq{ii}.poscut(poscut ~= 0) = poscut(poscut ~= 0);
  end
end

disp('Adding maps')
map = NaN(length(y_bin),length(x_bin),528);
map_std = NaN(length(y_bin),length(x_bin),528);

for jj = 1:528
  map_sum = NaN(length(y_bin),length(x_bin),length(qq));
  sum_w = NaN(length(y_bin),length(x_bin),length(qq));
  % Sum over how many maps 
  for ii = 1:length(qq)
    switch experiment
      case 'bicep2'
	map_sum(:,:,ii) = qq{ii}.map(:,:,jj).*qq{ii}.weightmask(:,:,jj);
      case 'keck'
	map_sum(:,:,ii) = qq{ii}.map(:,:,jj).*qq{ii}.weightmask(:,:,jj).*qq{ii}.poscut(jj)*qq{ii}.badmapcut(jj);
    end
    sum_w(:,:,ii) = qq{ii}.weightmask(:,:,jj);    
  end
  
  % map(:,:,jj)=nansum(map_sum,3)./nansum(sum_w,3);
  map(:,:,jj) = nanmedian(map_sum,3);
  wmap(:,:,jj) = nansum(sum_w,3);
  map_std(:,:,jj) = nanstd(map_sum,[],3);
end

% Crop if needed
if crop
  sizeindeg_x = 8.1; 
  sizeindeg_y = 8.1; 
  Field_size_deg = [sizeindeg_x sizeindeg_y ];
  N_pix = [length(x_bin) length(y_bin)];
  ad = calc_ad2(Field_size_deg,N_pix);
  save(['beammaps/maps_cropped/' filename],'map','ad','map_std','qq')
elseif cutmaps
  save(['beammaps/maps_coadded/' filename '_cut'],'map','x_bin','y_bin','wmap','map_std')
else
  save(['beammaps/maps_coadded/' filename],'map','x_bin','y_bin','wmap','map_std')
end

if plotcomp
  if cutmaps
    ffbm_plotCompMaps(experiment,year,rxNum,x_bin,y_bin,map,map_std,qq,subdir,mapstocut)
  else
    ffbm_plotCompMaps(experiment,year,rxNum,x_bin,y_bin,map,map_std,qq,subdir)
  end
end

% Save out stack of files
if savecomp 
  sizeindeg = 8;
  
  % The ad structure looks like this
  winaz = x_bin >= -sizeindeg/2 & x_bin <= sizeindeg/2;
  winel = y_bin >= -sizeindeg/2 & y_bin <= sizeindeg/2;
  
  win_x = find(winaz);
  win_x_bin = x_bin(win_x);
  win_y = find(winel);
  win_y_bin = y_bin(win_y);
  
  sizeindeg_x = 8.1; 
  sizeindeg_y = 8.1; 
  Field_size_deg = [sizeindeg_x sizeindeg_y ];
  N_pix = [length(win_x) length(win_y)];
  ad = calc_ad2(Field_size_deg,N_pix);
  
  map = NaN(81,81,528,length(qq));
  poscut = NaN(528,length(qq));
  dk = NaN(1,length(qq));
  
  for ii = 1:length(qq)
    map(:,:,:,ii) = qq{ii}.map(win_y(1):win_y(end),win_x(1):win_x(end),:);
    dk(ii) = qq{ii}.dk;
    poscut(:,ii) = qq{ii}.poscut;
    badmapcut(:,ii) = qq{ii}.badmapcut;
  end
  
  if ~cutmaps   
    save(['beammaps/maps_cropped/' filename '_comp'],'map','ad','dk','poscut','badmapcut') 
  else 
    save(['beammaps/maps_cropped/' filename '_comp'],'map','ad','dk','poscut','badmapcut','mapstocut')
  end
end
  
return


      
%}
