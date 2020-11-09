function run_b2bmsim(bmsimtype,bmsize,makepairmaps,coaddpairmaps,submit)
% function run_b2bmsim(bmsimtype,bmsize,makepairmaps,coaddpairmaps,submit)
% 
% Function to run all standard BICEP2 beam map sims for our "archival" beam
% map sim set.  Idea is to match real data maps as much as possible (in
% terms of sernums, daughters, etc)
% 
% bmsimtype: 'standard', 'split', 'floor'
% bmsize:      1.2, 2, 4, 6, 8

year = 2012;

[nbase daughter] = get_simnum(bmsimtype,bmsize);
rlz = 0;
type = 1; % All beam map sims are of the form sernum0001 
 
% Memory usage/time/queues
mempairmap = 40000; 
memcoadd = 30000;
maxtime = 400;
slurm_queue = 'general,serial_requeue,itc_cluster,regal';

% Global options
om = 1; % onlymissing
ntpg = 4; % ntags per job

% BICEP2 sim options
rlzc = 3; % realization chunk size
tagc = 1; % tag chunk size

% Coadd options
js = 0; % farm jacks separately for coadd
farmdps = 1;

% Deprojection options
deprojs = {[0,0,0,0],
           [1,0,0,0],
           [0,1,0,0],
           [1,1,0,0],
           [1,1,1,0],
           [1,1,0,1],
           [1,1,1,1],
           [1,1,0,2],
           [1,1,1,2]}; 

% Choose tags - should just be daughter-dependent
realfile = ['maps/1450/real_' daughter ...
            '_filtp3_weight3_gs_dp1100_jack0.mat'];
x = load(realfile);
tags = x.coaddopt.tags;
[tagsublist mtl] = get_tag_sublist(x.coaddopt);
clear x;

% Simopts
clear simopt
simopt.siginterp = 'linear'; % for beam map sims;
simopt.curveskyrotbeam = 1; 
simopt.curveskyrescale = 1; 
simopt.beamcen = 'obs';
simopt.diffpoint = 'ideal'; % Since the beam maps include this
simopt.chi = 'obs';
simopt.epsilon = 'obs';
simopt.noise = 'none';
simopt.sig = 'nopol';
simopt.sigmaptype = 'healmap';
simopt.maketod = false;
simopt.interpix = 0.1;
% Taken directly from reduc_coaddpairmaps.m, line 692
[p ind] = get_array_info([num2str(year) '0501'],...
                         'obs','obs','obs','obs',[],'obs','obs');
p.ukpv = (p.ukpv(ind.a)+p.ukpv(ind.b))/2;
simopt.ukpervolt = p.ukpv;
simopt.coord = 'C';
%simopt.force_ab_int_from_common_pix = false; % What is this?
simopt.update = 1;

% Mapopts
clear mapopt
mapopt.beamcen = 'obs';
mapopt.chi = 'obs';
mapopt.epsilon = 'obs';
mapopt.acpack = 0;
mapopt.gs = 1;
mapopt.filt = 'p3';
mapopt.deproj = 1:6; % true; % turn on 9,10 for residual beam subtraction
mapopt.update = true;
mapopt.curveskyrotbeam = 1; 
mapopt.curveskyrescale = 1; 
mapopt.realpairmapset = 'pairmaps/1450/real';

% Coaddopts
chflags = get_default_chflags([],num2str(year));

% Beam map flags - CLW-defined
% Doesn't seem to work with new get_array_info?  My maps probably don't
% include 'pair41' anyway (check this) so we don't need beam map flags
%chflags.filebase{9} = 'fp_data/fp_data_tocut';
%chflags.par_name{9} = 'pair41';
%chflags.low(9) = -1;
%chflags.high(9) = 0.5;

coaddopt.chflags = chflags;
coaddopt.daughter = daughter;
coaddopt.deproj_timescale = 'by_phase';
coaddopt.filt = 'p3';
coaddopt.gs = 1;
coaddopt.jacktype = '0123456789abcde';
coaddopt.temporaljack_splitdate = '20111101';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Make pairmaps

if makepairmaps
    
  [sigmapfilename mapopt.deproj_map] = get_inputdeprojmap(bmsimtype);
  simopt.mapopt = mapopt;
  
  % Choose your beam map!
  simopt.beammapfilename = get_beammap(bmsimtype,bmsize);
  
  if ~exist(['simrunfiles/' num2str(nbase) '_' daughter '_dithers.mat'],'file')
    par_make_simrunfiles(tags,simopt,nbase,daughter);
  end

  mkdir(['pairmaps/' num2str(nbase)]);
  system(['ln -s /n/panlfs2/bicep/bicep2/pipeline/pairmaps/1450/real pairmaps/' ...
	num2str(nbase) '/real'])
  
  % Check that this works with the latest version of runsim
  runsim(nbase,daughter,mapopt,'nopol','none',type,sigmapfilename,rlz,mtl,...
         true,rlzc,om,800,0,[],1,slurm_queue,0,mempairmap,tagc,0,maxtime,submit);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Coadd pairmaps

if coaddpairmaps
   
  coaddopt.sernum = [num2str(nbase) 'xxx1'];
  coaddopt.save_cuts = false;
  coaddopt.tagsublist = tagsublist;
  
  if farmdps
    for jj = 1:size(deprojs,1)
      coaddopt.deproj = deprojs{jj}; %coaddopt.deproj = deprojs(jj,:);
      farm_coaddpairmaps(tags,coaddopt,'OnlyMissing',om,'Realizations',rlz,...
                         'FarmJobs',1,'JobLimit',800,'FarmJacksSeparately',js,...
                         'Queue',slurm_queue,'MemRequire',memcoadd,...
                         'SplitSubmission',20,'UseCompiled',0,...
                         'maxtime',maxtime,'submit',submit);
    end
  else
    coaddopt.deproj = deprojs;
    farm_coaddpairmaps(tags,coaddopt,'OnlyMissing',om,'Realizations',rlz,...
                       'FarmJobs',1,'JobLimit',800,'FarmJacksSeparately',js,...
                       'Queue',slurm_queue,'MemRequire',memcoadd,...
                       'SplitSubmission',20,'UseCompiled',0,...
                       'maxtime',maxtime,'submit',submit);
  end

end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [nbase daughter]  = get_simnum(bmsimtype,bmsize)

% Choose sernum based on type/size
switch bmsimtype
  case 'standard'
    switch bmsize
      case 1.2
        nbase = 3625;
      case 2
        nbase = 3626;
      case 4
        nbase = 3627;
      case 6
        nbase = 3628;
      case 8
        nbase = 3629;
    end
  case 'split'
    switch bmsize
      case 1.2
        nbase = 3635;
      case 2
        nbase = 3636;
      case 4
        nbase = 3637;
      case 6
        nbase = 3638;
      case 8
        nbase = 3639;
    end
  case 'floor' % Constructed Gaussian maps from measured beamparams
    switch bmsize
      case 1.2
        nbase = 3645;
      case 2
        nbase = 3646;
      case 4
        nbase = 3647;
      case 6
        nbase = 3648;
      case 8
        nbase = 3649;
    end
end

% For BICEP2, daughter is just 'a' to match real
daughter = 'a';

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [sigmap deproj_map] = get_inputdeprojmap(bmsimtype)

sigmap = ...
    {'input_maps/planck/planck_derivs_nopix/HFI_SkyMap_100_2048_R1.10_nominal_nside0512_nlmax1280_bnone.fits'};

% Deprojection map -- for regular sims, this is the beam profile-smoothed
% Planck map.  For constructed Gaussian beam maps, we use
% Gaussian-smoothed Planck maps instead (so deprojection should work
% perfectly-ish)
switch bmsimtype
  case {'standard','split'}
    deproj_map = ...
        {'input_maps/planck/planck_derivs_nopix/synfast_deproj_143_nominal_B2.fits'};
  case 'floor'
    deproj_map = ...
        {'input_maps/planck/planck_derivs_nopix/HFI_SkyMap_143_2048_R1.10_nominal_nside0512_nlmax1280_b30.52.fits'};
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function filename = get_beammap(bmsimtype,bmsize)

sizestr = [strrep(num2str(bmsize),'.','p') 'deg'];

switch bmsimtype
  case 'standard' 
    filename = ['beammaps/maps_composite/ffbm_2012' ...
                '_all_' sizestr '_dk0.mat'];
  case 'split'
    filename = ['beammaps/maps_split/ffbm_2012' ...
                '_half_' sizestr '_dk0_diff.mat'];
  case 'floor'
    filename = ['beammaps/maps_constructed/bicep2_2012_obs_' sizestr '.mat'];
end

return

