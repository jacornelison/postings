function ffbm_maskground(compopt)
% function ffbm_maskground(compopt)
%
% Function to take windowed component map, shift pair to center of map, and
% mask out the ground/SPT.
% Mask everything in the component struct, ignoring cuts/etc.
%
% INPUTS (should all be sent in with compopt)
%
%  expt:           'bicep2','keck','bicep3'
%  year:           bicep2: 2012
%                  keck: 2012-2016
%                  bicep3: 2015
%
% OPTIONAL INPUTS
%
%  mapsize:        map size in deg (2 default), used for map filename
%  componentdir:   directory from which to load map struct
%                  (default beammaps/maps_component)
%  componentfile:  file in componentdir to load
%                  (default ffbm_year_allcomp)
%  maskeddir:      directory in which to save masked maps
%                  (default beammaps/maps_masked)
%  maskedfile:     file in maskeddir to save to
%                  (default ffbm_year_allcomp_masked)
%  cent:           center either on common A/B centroid or per-det
%                  'ab_centroid' (default), 'detector'
%  suffix:          optional string to add to map name, '' default

% Parse compopt
expt = compopt.expt;
year = compopt.year;
if ~isfield(compopt,'mapsize');
  sizestr = '2deg';
else
  sizestr = [strrep(num2str(compopt.mapsize),'.','p') 'deg'];
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
  compopt.componentfile = ['ffbm_' num2str(year) '_allcomp'...
	'_' sizestr suffix];
  componentfile = compopt.componentfile;
else
  componentfile = compopt.componentfile;
end
if ~isfield(compopt,'maskeddir')
  compopt.maskeddir = 'beammaps/maps_masked';
  maskeddir = compopt.maskeddir;
else
  maskeddir = compopt.maskeddir;
end
if ~isfield(compopt,'maskedfile')
  compopt.maskedfile = ['ffbm_' num2str(year) '_allcomp_masked'...
	'_' sizestr suffix];
  maskedfile = compopt.maskedfile;
else
  maskedfile = compopt.maskedfile;
end
if ~isfield(compopt,'cent')
  compopt.cent = 'ab_centroid';
  cent = compopt.cent;
else
  cent = compopt.cent;
end

[p ind] = get_array_info([num2str(year) '0201']);

% Load up map file
filename = [componentdir '/' componentfile];
load(filename);

% Go through all light pairs and do yo thang
for ii = 1:length(ind.la)
  if mod(ii,10) == 0
    disp(['Masking pair ' num2str(ii) ' / ' num2str(length(ind.la))])
  end
  for jj = 1:length(comp.bm.number)

    fita = comp.map{ind.la(ii)}.fit{jj};
    fitb = comp.map{ind.lb(ii)}.fit{jj};
    x_bin = comp.az_ap{jj};
    y_bin = comp.el_ap{jj};
    % Take a map of a source centered at arbitrary az/el and use the fits 
    % to shift the axes so (0,0) is at the pair centroid or beam center
    [xx_a yy_a xx_b yy_b fita fitb] = shiftaxes(fita,fitb,cent,x_bin,y_bin);
    
    % Use the re-zeroed grids to interpolate to a new map with (0,0) 
    % at the center
    % 
    % HERE WE DEFINE THE AD STRUCT
    % ad.t_val_deg{1} = increasing apparent az, zero-centered
    % ad.t_val_deg{2} = increasing apparent el, zero-centered

    [X,Y] = meshgrid(comp.ad.t_val_deg{1},comp.ad.t_val_deg{2});
    [XX_a,YY_a] = meshgrid(xx_a,yy_a);
    [XX_b,YY_b] = meshgrid(xx_b,yy_b);

    % Recenter
    comp.map{ind.la(ii)}.component{jj} = interp2(XX_a,YY_a,...
	comp.map{ind.la(ii)}.component{jj},X,Y);
    comp.map{ind.lb(ii)}.component{jj} = interp2(XX_b,YY_b,...
	comp.map{ind.lb(ii)}.component{jj},X,Y);
    
    % Again, to get a standard beam map, we now do:
    %   imagesc(ad.t_val_deg{1},ad.t_val_deg{2},map)
    %   set(gca,'ydir','normal')
    % and the result is a map with apparent az increasing left to right
    % and apparent el increasing going up
    % But the parity is opposite that of x'/y'
    
    % Now that the maps are shifted, we NO LONGER USE AZ_AP/EL_AP
    % even though ad.t_val_deg is in the same directions
    
    % Move fit values to match map;
    comp.map{ind.la(ii)}.fit{jj} = fita;
    comp.map{ind.lb(ii)}.fit{jj} = fitb;
    
    % In principle we would make separate masks for both detectors if we
    % center on the common AB centroid, since we calculate masks by
    % distance from the main beam.  However, A/B offsets are so small (<1
    % arcminute worst-case scenario) that it's not worth it.
    
    % Mask out ground completely.  True (1) means we cut it out
    % az/el center are (0,0)
    azlimcalc = zeros(size(comp.ad.t_val_deg{1}));
    ellimcalc = comp.ad.t_val_deg{2} < -1.5;
    [calcmaskaz calcmaskel] = meshgrid(azlimcalc,ellimcalc);
    groundmask = calcmaskaz | calcmaskel;
    comp.map{ind.la(ii)}.component{jj}(groundmask) = NaN;
    comp.map{ind.lb(ii)}.component{jj}(groundmask) = NaN;
    
    % If Keck, mask out SPT
    switch expt
      case 'keck'
	azlimcalc = comp.ad.t_val_deg{1} > -7 & comp.ad.t_val_deg{1} < -2;
	ellimcalc = comp.ad.t_val_deg{1} < 1;
	[calcmaskaz calcmaskel] = meshgrid(azlimcalc,ellimcalc);
	% Here we need both to be true
	sptmask = calcmaskaz & calcmaskel;
	comp.map{ind.la(ii)}.component{jj}(sptmask) = NaN;
	comp.map{ind.lb(ii)}.component{jj}(sptmask) = NaN;
    end
    
    % Mask the mast
    % How much to eat into the main beam depends on the beam size
    % Semi-arbitrary, based on looking at the maps
    % The mast is mostly obvious in Keck 220 maps - hard to see in B3
    % 2016.  However let's be safe!  It has a similar beamwidth to Keck
    % 220 so we'll go with that distance.
    % Keck 270s -- match 220s for now, but may move mask closer in later.
    switch expt
      case 'keck'
        switch p.band(ind.la(ii))
          case 100
            ellimcalc = comp.ad.t_val_deg{2} < -1.2;
          case 150
            ellimcalc = comp.ad.t_val_deg{2} < -1;
          case {210,220,270}
            ellimcalc = comp.ad.t_val_deg{2} < -0.5;
        end
      case 'bicep3'
        ellimcalc = comp.ad.t_val_deg{2} < -0.5;
    end
    azlimcalc = comp.ad.t_val_deg{1} > -0.3 & comp.ad.t_val_deg{1} < 0.3;
    [calcmaskaz calcmaskel] = meshgrid(azlimcalc,ellimcalc);
    % Here we need both to be true
    mastmask = calcmaskaz & calcmaskel;
    comp.map{ind.la(ii)}.component{jj}(mastmask) = NaN;
    comp.map{ind.lb(ii)}.component{jj}(mastmask) = NaN;
  end
  
end

comp.compopt = compopt;

savename = [maskeddir '/' maskedfile];
save(savename,'comp','-v7.3');

disp(['Saved masked component file: ' savename]);

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [xx_a yy_a xx_b yy_b fita fitb] = shiftaxes(fita,fitb,cent,x_bin,y_bin)
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
fita(2) = fita(2) - xcent_a;
fitb(2) = fitb(2) - xcent_b;
fita(3) = fita(3) - ycent_a;
fitb(3) = fitb(3) - ycent_b;

return