function composite_runall_keck_2012(winsize,makesingle,unrotplots,maskmaps,...
    maskplots,rotmaps,rotplots,compmaps,compplots,splitmaps)
% composite_runall_keck_2012(winsize,makesingle,unrotplots,maskmaps,...
%                            maskplots,rotmaps,rotplots,compmaps,...
%                            compplots,splitmaps)
%
% Wrapper function to go from windowed/fitted maps and cuts to final
% composite beam maps suitable for beam map sims.  
%
% Before running this: 
%   1. ALL windowed maps in the bmrunlist should be have been made with
%      bm_makemap (the fits should all exist in the files)
%   2. ffbm_makecuts should have been run to generate the automated cuts
%   3. ffbm_makehandcuts should have been run to add visually bad maps to
%      the cut struct
%
% Note that here we're using all default options

compopt.expt = 'keck';
compopt.year = 2012;
switch winsize
  case 2
    compopt.suffix = '';
    compopt.componentfile = 'ffbm_2012_allcomp_2deg';
  case 8
    compopt.suffix = '_8deg';
end

% Gather all windowed map files from beammaps/maps into a single file
%  -> beammaps/maps_component/ffbm_2012_allcomp.mat
compopt.onlychflag = 0;
if makesingle
  ffbm_makesinglecompfile(compopt);
end

if winsize == 2
  compopt.suffix = '_2deg';
end

% Plot unrotated maps - thumbnail is total number of maps for
% this detector
compopt.plottype = 'component';
compopt.plotdir = ['plots/keck_2012_compmaps_unrotated' compopt.suffix];
if unrotplots
  ffbm_plotCompMaps(compopt);
end

% Mask the maps - take unrotated
%   beammaps/maps_component/ffbm_2012_allcomp.mat ->
%   beammaps/maps_masked/ffbm_2012_allcomp_masked.mat
if maskmaps
  ffbm_maskground(compopt);
end

% Plot masked maps
compopt.plottype = 'masked';
compopt.plotdir = ['plots/keck_2012_compmaps_masked' compopt.suffix];
compopt.componentdir = 'beammaps/maps_masked';
compopt.componentfile = ['ffbm_2012_allcomp_masked' compopt.suffix];
compopt.plotfitxy = 0;
if maskplots
  ffbm_plotCompMaps(compopt);
end

% Rotate maps - take masked
%   beammaps/maps_masked/ffbm_2012_allcomp_masked.mat ->
%   beammaps/maps_rotated/ffbm_2012_allcomp_rotated.mat
if rotmaps
  % First rotate to dk0
  compopt.coord = 'dk0';
  ffbm_rotatemaps(compopt);
  % Saves to ffbm_2012_allcomp_rotated_8deg_dk0
  compopt.coord = 'xpyp';
  ffbm_rotatemaps(compopt);
  % Saves to ffbm_2012_allcomp_rotated_8deg_xpyp
end

% Plot rotated maps
compopt.plottype = 'rotated';
compopt.componentdir = 'beammaps/maps_rotated';
if rotplots
  compopt.plotdir = ['plots/keck_2012_compmaps_rotated' ...
                     compopt.suffix '_dk0'];
  compopt.componentfile = ['ffbm_2012_allcomp_rotated' ...
                      compopt.suffix '_dk0'];
  compopt.coord = 'dk0';
  ffbm_plotCompMaps(compopt);
  compopt.plotdir = ['plots/keck_2012_compmaps_rotated' ...
                     compopt.suffix '_xpyp'];
  compopt.componentfile = ['ffbm_2012_allcomp_rotated' ...
                      compopt.suffix '_xpyp'];
  compopt.coord = 'xpyp';
  ffbm_plotCompMaps(compopt);
end

% Composite maps - take rotated
%   beammaps/maps_rotated/ffbm_2012_allcomp_rotated.mat ->
%   beammaps/maps_composite/ffbm_2012_all_onlycommonab.mat
compopt.onlycommon_ab = 1;
compopt.radtozero = winsize;
if compmaps
  % First composite dk0
  compopt.coord = 'dk0';
  compopt.rotatedfile = ['ffbm_2012_allcomp_rotated' ...
                      compopt.suffix '_dk0'];
  compopt.compositefile = ['ffbm_2012_all_onlycommonab' ...
                      compopt.suffix '_dk0'];
  ffbm_compositemaps(compopt);
  % Now xpyp
  compopt.coord = 'xpyp';
  compopt.rotatedfile = ['ffbm_2012_allcomp_rotated' ...
                      compopt.suffix '_xpyp'];
  compopt.compositefile = ['ffbm_2012_all_onlycommonab' ...
                      compopt.suffix '_xpyp'];
  ffbm_compositemaps(compopt);
end

% Plot composite maps
% First plot dk0, then xpyp
compopt.plottype = 'composite';
if compplots
  compopt.coord = 'dk0';
  compopt.plotdir = ['plots/keck_2012_composites' ...
                   compopt.suffix '_dk0'];
  compopt.compositefile = ['ffbm_2012_all_onlycommonab' ...
                    compopt.suffix '_dk0'];
  compopt.componentfile = ['ffbm_2012_allcomp_rotated' ...
                   compopt.suffix '_dk0'];
  ffbm_plotCompMaps(compopt);
  compopt.coord = 'xpyp';
  compopt.plotdir = ['plots/keck_2012_composites' ...
                     compopt.suffix '_xpyp'];
  compopt.compositefile = ['ffbm_2012_all_onlycommonab' ...
                      compopt.suffix '_xpyp'];
  compopt.componentfile = ['ffbm_2012_allcomp_rotated' ...
                      compopt.suffix '_xpyp'];  
  ffbm_plotCompMaps(compopt);
end

% Make split/sum/diff maps
compopt.makesplits = 1;
if splitmaps
  ffbm_compositemaps(compopt);
end

return




