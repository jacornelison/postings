function composite_runall_b3sched40_2015(makesingle,unrotplots,maskmaps,...
    maskplots,rotmaps,rotplots,compmaps,compplots,splitmaps)
% composite_runall_b3sched40_2015(makesingle,unrotplots,maskmaps,...
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

compopt.expt = 'b3';
compopt.year = 2015;
compopt.mapstoload = [63,64,65,66,68,69,70,71];

% Gather all windowed map files from beammaps/maps into a single file
%  -> beammaps/maps_component/ffbm_2015_sched40.mat
compopt.componentfile = 'ffbm_2015_sched40';
if makesingle
  ffbm_makesinglecompfile(compopt);
end

% Plot unrotated maps - thumbnail is total number of maps for
% this detector
compopt.hascomposite = 0;
compopt.plotdir = 'plots/b3_2015_compmaps_unrotated'; 
bmopt.expt = 'b3';
bmopt.rxNum = '0';
bmopt.filename = '20150201';
bmopt.plotdir = compopt.plotdir;
bmopt.t1 = '';
bmopt.t2 = '';
bmopt.author = 'KSK';
bmopt.bmnum = '';
bmopt.sched = '';
bmopt.run = 'b3_2015';

if unrotplots
  ffbm_plotCompMaps(compopt);
  make_bmhtml(bmopt);
end

% Mask the maps - take unrotated
%   beammaps/maps_component/ffbm_2015_sched40.mat ->
%   beammaps/maps_masked/ffbm_2015_sched40_masked.mat
compopt.maskedfile = 'ffbm_2015_sched40_masked';
if maskmaps
  ffbm_maskground(compopt);
end

% Plot masked maps
compopt.plotdir = 'plots/b3_2015_compmaps_masked';
compopt.componentdir = 'beammaps/maps_masked';
compopt.componentfile = 'ffbm_2015_sched40_masked';
bmopt.plotdir = compopt.plotdir;
if maskplots
  ffbm_plotCompMaps(compopt);
  make_bmhtml(bmopt);
end

% Rotate maps - take masked
%   beammaps/maps_masked/ffbm_2015_sched40_masked.mat ->
%   beammaps/maps_rotated/ffbm_2015_sched40_rotated.mat
compopt.rotatedfile = 'ffbm_2015_sched40_rotated';
if rotmaps
  ffbm_rotatemaps(compopt);
end

% Plot rotated maps
compopt.plotdir = 'plots/b3_2015_compmaps_rotated';
compopt.componentdir = 'beammaps/maps_rotated';
compopt.componentfile = 'ffbm_2015_sched40_rotated';
bmopt.plotdir = compopt.plotdir;
if rotplots
  ffbm_plotCompMaps(compopt);
  make_bmhtml(bmopt);
end

% Composite maps - take rotated
%   beammaps/maps_rotated/ffbm_2015_sched40_rotated.mat ->
%   beammaps/maps_composite/ffbm_2014_sched40_onlycommonab.mat
compopt.onlycommon_ab = 0;
compopt.compositefile = 'ffbm_2015_sched40';
if compmaps
  ffbm_compositemaps(compopt);
end

% Plot composite maps
compopt.plotdir = 'plots/b3_2015_composites';
compopt.hascomposite = 1;
bmopt.plotdir = compopt.plotdir;
if compplots
  ffbm_plotCompMaps(compopt);
  %make_bmhtml(bmopt);
end

% Make split/sum/diff maps
compopt.makesplits = 1;
if splitmaps
  ffbm_compositemaps(compopt);
end

return




