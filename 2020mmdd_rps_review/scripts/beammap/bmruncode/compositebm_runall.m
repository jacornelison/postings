function compositebm_runall(expt,mapsize,year,makesingle,maskmaps, ...
                            rotmaps,compmaps,splitmaps,...
                            xyearrotmaps,xyearcompmaps,xyearsplitmaps);
% function compositebm_runall(expt,mapsize,year,makesingle,maskmaps, ...
%                             rotmaps,compmaps,splitmaps,...
%                             xyearrotmaps,xyearcompmaps,xyearplitmaps);
%
% General run script to go from individual-schedule windowed maps to
% composites, supplanting expt/year-specific run scripts.
%
%   1. ALL windowed maps in the bmrunlist should be have been made with
%      bm_winfitmap
%   2. ffbm_makecuts should have been run to generate the automated cuts
%   3. ffbm_makehandcuts should have been run to add visually bad maps to
%      the cut struct
% 
%   For each type of map: 1 = make the map, don't plot
%                         2 = make the map and plot
%                         3 = plot, don't make the map
%   
%   expt: 'bicep2','keck','bicep3'
%   year: 2012, etc.
%   mapsize: 1.2, 2, 4, ... (looks for these windowed maps)

compopt.expt = expt; 
compopt.year = year; 
if strcmp(expt,'bicep2') % BICEP2 has only one set, xyear not possible
  year = 2012;
  xyearrotmaps = 0;
  xyearcompmaps = 0;
  xyearsplitmaps = 0;
end
compopt.mapsize = mapsize; 
compopt.suffix = '';
sizestr = [strrep(num2str(mapsize),'.','p') 'deg'];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Gather all windowed map files from beammaps/maps into a single file
%  -> beammaps/maps_component/ffbm_year_allcomp_mapsizedeg.mat
compopt.onlychflag = 0;
if makesingle == 1 | makesingle == 2
  ffbm_makesinglecompfile(compopt);
end

% Plot unrotated maps - thumbnail is total number of maps for
% this detector
compopt.plottype = 'component';
compopt.plotdir = ['plots/' expt '_' num2str(year) ...
                   '_compmaps_unrotated_' sizestr];
if makesingle == 2 | makesingle == 3
  ffbm_plotCompMaps(compopt);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Mask the maps - take unrotated
%   beammaps/maps_component/ffbm_year_allcomp_mapsizedeg.mat ->
%   beammaps/maps_masked/ffbm_year_allcomp_masked_mapsizedeg.mat
if maskmaps == 1 | maskmaps == 2
  ffbm_maskground(compopt);
end

% Plot masked maps
compopt.plottype = 'masked';
compopt.plotdir = ['plots/' expt '_' num2str(year) ...
                   '_compmaps_masked_' sizestr];
compopt.componentdir = 'beammaps/maps_masked';
compopt.componentfile = ['ffbm_' num2str(year) '_allcomp_masked_' ...
                         sizestr];
compopt.plotfitxy = 0;
if maskmaps == 2 | maskmaps == 3
  ffbm_plotCompMaps(compopt);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Rotate maps - take masked
%   beammaps/maps_masked/ffbm_year_allcomp_masked_mapsizedeg.mat ->
%   beammaps/maps_rotated/ffbm_year_allcomp_rotated_mapsizedeg_coord.mat
if rotmaps == 1 | rotmaps == 2
  % First rotate to dk0
  compopt.coord = 'dk0';
  ffbm_rotatemaps(compopt);
  % Saves to ffbm_year_allcomp_rotated_mapsizedeg_dk0
  compopt.coord = 'xpyp';
  ffbm_rotatemaps(compopt);
  % Saves to ffbm_year_allcomp_rotated_mapsizedeg_xpyp
end

% Plot rotated maps
compopt.plottype = 'rotated';
compopt.componentdir = 'beammaps/maps_rotated';
if rotmaps == 2 | rotmaps == 3
  compopt.plotdir = ['plots/' expt '_' num2str(year) ...
                     '_compmaps_rotated_' sizestr '_dk0'];
  compopt.componentfile = ['ffbm_' num2str(year) ...
                      '_allcomp_rotated_' sizestr '_dk0']; 
  compopt.coord = 'dk0';
  ffbm_plotCompMaps(compopt);
  compopt.plotdir = ['plots/' expt '_' num2str(year) ...
                     '_compmaps_rotated_' sizestr '_xpyp'];
  compopt.componentfile = ['ffbm_' num2str(year) ...
                      '_allcomp_rotated_' sizestr '_xpyp']; 
  compopt.coord = 'xpyp';
  ffbm_plotCompMaps(compopt);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Composite maps - take rotated
%   beammaps/maps_rotated/ffbm_year_allcomp_rotated_mapsizedeg_coord.mat ->
%   beammaps/maps_composite/ffbm_year_all_mapsizedeg_coord.mat
compopt.onlycommon_ab = 1;
%compopt.radtozero = winsize;
if compmaps == 1 | compmaps == 2
  % First composite dk0
  compopt.coord = 'dk0';
  ffbm_compositemaps(compopt);
  % Now xpyp
  compopt.coord = 'xpyp';
  ffbm_compositemaps(compopt);
end

% Plot composite maps
% First plot dk0, then xpyp
compopt.plottype = 'composite';
if compmaps == 2 | compmaps == 3
  compopt.coord = 'xpyp';
  compopt.plotdir = ['plots/' expt '_' num2str(year) ...
	'_composites_' sizestr '_xpyp'];
  compopt.compositefile = ['ffbm_' num2str(year) '_all_' ...
	sizestr '_xpyp'];
  compopt.componentfile = ['ffbm_' num2str(year)  ...
	'_allcomp_rotated_' sizestr '_xpyp'];
  ffbm_plotCompMaps(compopt);
  compopt.coord = 'dk0';
  compopt.plotdir = ['plots/' expt '_' num2str(year) ...
                     '_composites_' sizestr '_dk0'];
  compopt.compositefile = ['ffbm_' num2str(year) '_all_' ...
                      sizestr '_dk0'];
  compopt.componentfile = ['ffbm_' num2str(year) ...
                      '_allcomp_rotated_' sizestr '_dk0']; 
  ffbm_plotCompMaps(compopt);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Make split/sum/diff maps
%   beammaps/maps_rotated/ffbm_year_allcomp_rotated_mapsizedeg_coord.mat ->
%   beammaps/maps_split/ffbm_year_half_mapsizedeg_coord_{1/2/halfsum/halfdiff}
compopt.makesplits = 1;
if splitmaps == 1 | splitmaps == 2
  compopt.coord = 'dk0';
  ffbm_compositemaps(compopt);
  compopt.coord = 'xpyp';
  ffbm_compositemaps(compopt);
end

% Plot split maps
if splitmaps == 2 | splitmaps == 3
end

compopt.makesplits = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Make x-year rotated maps (needs all relevant rotated years to exist)
% - take rotated
%   beammaps/maps_rotated/ffbm_year_allcomp_rotated_mapsizedeg_coord.mat
%   beammaps/maps_rotated/ffbm_year_allcomp_rotated_xyear_mapsizedeg_coord
if xyearrotmaps == 1 | xyearrotmaps == 2
  % First rotate to dk0
  compopt.coord = 'dk0';
  ffbm_xyearcombine(compopt);
  % Saves to ffbm_year_allcomp_rotated_xyear_mapsizedeg_dk0
  compopt.coord = 'xpyp';
  ffbm_xyearcombine(compopt);
  % Saves to ffbm_year_allcomp_rotated_xyear_mapsizedeg_xpyp
end

% Plot xyear maps
if xyearrotmaps == 2 | xyearrotmaps == 3
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% x-year Composite maps - take rotated
%   beammaps/maps_rotated/ffbm_year_allcomp_rotated_xyear_mapsizedeg_coord ->
%   beammaps/maps_composite/ffbm_year_all_xyear_mapsizedeg_coord.mat
compopt.onlycommon_ab = 1;
compopt.xyear = 1;
%compopt.radtozero = winsize;
if xyearcompmaps == 1 | xyearcompmaps == 2
  % First composite dk0
  compopt.coord = 'dk0';
  ffbm_compositemaps(compopt);
  % Now xpyp
  compopt.coord = 'xpyp';
  ffbm_compositemaps(compopt);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Make xyear split/sum/diff maps
%   beammaps/maps_rotated/ffbm_year_allcomp_rotated_mapsizedeg_coord.mat ->
%   beammaps/maps_split/ffbm_year_half_xyear_mapsizedeg_coord_ ... 
%                       {1/2/halfsum/halfdiff}
compopt.makesplits = 1;
if xyearsplitmaps == 1 | xyearsplitmaps == 2
  compopt.coord = 'dk0';
  ffbm_compositemaps(compopt);
  compopt.coord = 'xpyp';
  ffbm_compositemaps(compopt);
end

% Plot xyear split maps
if xyearsplitmaps == 2 | xyearsplitmaps == 3
end

return