function ffbm_xyearcombine(compopt)
% function ffbm_xyearcombine(compopt)
%
% Take rotated component maps from all years and generate a new set of
% rotated maps (for each year) containing the full set of components
% which should be used in the composite, so that we use all available
% maps for the highest S/N composites.  The output should be compatible
% with the input to ffbm_compositemaps, except that all of the cuts will
% be implemented here instead and ffbm_compositemaps simply
% coadds/medians the maps.
%
% INPUTS (should all be sent in with compopt)
% 
%   expt:        'keck','b3'
%   year:        keck: 2012,2013,2014,2015
%                b3: 2016
%
% OPTIONAL INPUTS
%
%   mapsize:       map size in deg (2 default), used for map filename
%                  This refers to the radius from beam center - masked
%                  maps are actually square, so in this function we cut
%                  down to a circle with this radius
%   coord:         Coordinate system (used for filename)
%                  'dk0' (default) is the standard input to reduc_makesim
%                  'xpyp' for normal visualization (x'/y')
%   rotateddir:    directory from which to load rotated maps
%                  (default beammaps/maps_rotated)
%   rotatedfile:   file in rotateddir to load
%                  (default ffbm_year_allcomp_rotated_mapsizedeg_coord)
%   xyearfile:     file in rotateddir to save to
%                  (default ffbm_year_allcomp_rotated_xyear_mapsizedeg_coord)
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
%   suffix:        string to add to end

% Parse compopt
expt = compopt.expt;
year = compopt.year;
if ~isfield(compopt,'mapsize');
  sizestr = '2deg';
  mapsize = 2;
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
if ~isfield(compopt,'suffix')
  suffix = '';
else
  suffix = compopt.suffix;
end
if ~isfield(compopt,'rotateddir')
  compopt.rotateddir = 'beammaps/maps_rotated';
  rotateddir = compopt.rotateddir;
else
  rotateddir = compopt.rotateddir;
end
if ~isfield(compopt,'xyearfile')
  compopt.xyearfile = ['ffbm_' num2str(year) ...
                      '_allcomp_rotated_xyear_' ...
                      sizestr '_' coord suffix];
  xyearfile = compopt.xyearfile;
else
  xyearfile = compopt.xyearfile;
end
if ~isfield(compopt,'applycuts')
  compopt.applycuts = 1;
  applycuts = compopt.applycuts;
else
  applycuts = compopt.applycuts;
end
if ~isfield(compopt,'cutdir')
  compopt.cutdir = 'beammaps/cuts';
  cutdir = compopt.cutdir;
else
  cutdir = compopt.cutdir;
end
if ~isfield(compopt,'cutlist')
  compopt.cutlist = {'mirror','notlight','nanfit','peak',...
      'sigma','ellip','hand'};
  cutlist = compopt.cutlist;
else
  cutlist = compopt.cutlist;
end
if ~isfield(compopt,'onlycommon_ab')
  compopt.onlycommon_ab = 0;
  onlycommon_ab = compopt.onlycommon_ab;
else
  onlycommon_ab = compopt.onlycommon_ab;
end

[p ind] = get_array_info([num2str(year) '0201']);

% Get the master list and figure out which rotated maps files/cuts we
% want to load
[P,K] = ParameterRead('beammaps/bmcombine_master.csv');

% Cut down to this year
P = structcut(P,P.year == year);
% Find the relevant years
yearstoload = unique([P.map1 P.map2 P.map3 P.map4 P.map5 P.map6]);
yearstoload = yearstoload(~isnan(yearstoload));

for ii = 1:length(yearstoload)
  
  yearstoload(ii)
  % Get the right rotated maps and cut file
  rotatedfile = ['ffbm_' num2str(yearstoload(ii)) ...
                 '_allcomp_rotated_' sizestr '_' coord suffix];
  cutfile = ['ffbm_cuts_' num2str(yearstoload(ii))];

  filename = [rotateddir '/' rotatedfile];
  load(filename);
  
  if applycuts
    cutname = [cutdir '/' cutfile];
    load(cutname);
    cut = cuts.cuts;
  end
  
  % For all pairs, make a total cut mask
  for jj = 1:length(ind.la)
    
    totmask_a = ones(length(comp.bm.number),1);
    totmask_b = ones(length(comp.bm.number),1);
    
    for kk = 1:length(comp.bm.number)
      
      a = structcut(cut(kk),ind.la(jj));
      b = structcut(cut(kk),ind.lb(jj));
      totcut_a = 0;
      totcut_b = 0;
      for ll = 1:length(cutlist)
        totcut_a = totcut_a + getfield(a,cutlist{ll});
        totcut_b = totcut_b + getfield(b,cutlist{ll});
      end
      % Handcut == 1 means the fit may suck, but it's probably okay for
      % compositing
      if totcut_a == 0 | a.hand == 1 
	totmask_a(kk) = 0;
      end
      if totcut_b == 0 | b.hand == 1
	totmask_b(kk) = 0;
      end
    end
    
    % Here we go a slightly different route: to make the struct sizes
    % easier to deal with, instead of keeping only good maps, we NaN out
    % all bad maps so that in ffbm_compositemaps we don't have to worry
    % about cuts.  Sacrifice a bit on file size but gain a lot on sanity.
    
    % Indices of bad maps
    badinds_a = find(totmask_a ~= 0);
    badinds_b = find(totmask_b ~= 0);
    
    % Automatically do common A/B
    badinds_a = union(badinds_a,badinds_b);
    badinds_b = badinds_a;
    
    % Now check if this year's maps are even appropriate for the
    % requested rx/year combo!
    thisrx = structcut(P,P.rx == p.rx(ind.la(jj)));
    okyears = unique([thisrx.map1 thisrx.map2 thisrx.map3 ...
                      thisrx.map4 thisrx.map5 thisrx.map6]);
    okyears = okyears(~isnan(okyears));
    % If not on the acceptable year list, choose all indices as bad 
    if ~ismember(yearstoload(ii),okyears)
      badinds_a = 1:length(comp.bm.number);
      badinds_b = 1:length(comp.bm.number);
    end
    % Check for rx1 exception
    donotuse = rx1_exception(year,yearstoload(ii),...
                             p.rx(ind.la(jj)),p.tile(ind.la(jj)));
    if donotuse
      badinds_a = 1:length(comp.bm.number);
      badinds_b = 1:length(comp.bm.number);
    end
    
    % Go through rotated component maps and NaN out both the map and the
    % fit for all bad maps
    for kk = 1:length(badinds_a)
      comp.map{ind.la(jj)}.component{badinds_a(kk)} = ...
          NaN(size(comp.map{ind.la(jj)}.component{badinds_a(kk)}));
      comp.map{ind.lb(jj)}.component{badinds_a(kk)} = ...
          NaN(size(comp.map{ind.lb(jj)}.component{badinds_a(kk)}));
      comp.map{ind.la(jj)}.fit{badinds_a(kk)} = ...
          NaN(size(comp.map{ind.la(jj)}.fit{badinds_a(kk)}));
      comp.map{ind.lb(jj)}.fit{badinds_a(kk)} = ...
          NaN(size(comp.map{ind.lb(jj)}.fit{badinds_a(kk)}));
    end
    
  end % channel loop
    
  % Save the maps/etc  struct
  xyear{ii} = comp;
  
end 

clear comp

% Now smash them all together in an unholy mess

for ii = 1:length(ind.e)

  map{ii}.component = xyear{1}.map{ii}.component;
  map{ii}.fit = xyear{1}.map{ii}.fit;
  if length(xyear) > 1
      for jj = 2:length(xyear)
        map{ii}.component = ...
            [map{ii}.component xyear{jj}.map{ii}.component];
        map{ii}.fit = ...
            [map{ii}.fit xyear{jj}.map{ii}.fit];
      end
  end

end

% Don't need az_ap, el_ap fields
% But we do need a new bm struct with the total number
for ii = 1:length(xyear)
  len(ii) = length(xyear{ii}.az_ap);
end
comp.bm.number = 1:sum(len);

comp.map = map;
comp.ad = xyear{1}.ad;
comp.coord = xyear{1}.coord;

savename = [rotateddir '/' xyearfile];
save(savename,'comp','-v7.3');

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function donotuse = rx1_exception(xyear,thisyear,rx,tile)
% Rx1 Tile 1 was replaced in between 2012 and 2013

donotuse = 0;

% Turn on donotuse in 2 cases
switch rx
  case 1
    switch tile
      case 1
        switch xyear
          case 2012
            if (thisyear == 2013) | (thisyear == 2014)
              donotuse = 1;
            end
          case {2013,2014}
            if thisyear == 2012
              donotuse = 1;
            end
        end 
    end
end

return