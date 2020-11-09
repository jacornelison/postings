function run_keckbmsim(year,bmsimtype,bmsize,makepairmaps,coaddpairmaps,submit)
% function run_keckbmsim(year,bmsimtype,bmsize,makepairmaps,coaddpairmaps,submit)
% 
% Function to run all standard Keck beam map sims for our "archival" beam
% map sim set.  Idea is to match real data maps as much as possible (in
% terms of sernums, daughters, etc)
% 
% year:      2012, 2013, 2014, 2015
% bmsimtype: 'standard', 'split', 'floor'
% bmsize:      1.2, 2, 4, 6, 8

[nbase daughter] = get_simnum(year,bmsimtype,bmsize);
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

% Keck sim options
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
realfile = ['maps/1351/real_' daughter ...
            '_filtp3_weight3_gs_dp1100_jack01.mat'];
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
mapopt.realpairmapset = 'pairmaps/1351/real';

% Coaddopts
chflags = get_default_chflags([],num2str(year));

% Beam map flags.  For Keck 2012-2014, no need to implement additional
% beam map channel flags (see 20160824_k13bmsim_spectra and 
% 20160804_k2015bmsim_spectra).  However the 220s have some crappy beams
% which should really be taken out.
switch daughter
  case 'e'
    % 2015 format changed - 5xn_flags array, one perfreq
    chflags.filebase{10} = 'fp_data/fp_data_tocut';
    chflags.par_name{10} = 'nobeams_xyear';
    chflags.low(:,10) = -0.1;
    chflags.high(:,10) = 0.1;
end

coaddopt.chflags = chflags;
coaddopt.daughter = daughter;
coaddopt.deproj_timescale = 'by_phase';
coaddopt.filt = 'p3';
coaddopt.gs = 1;
coaddopt.jacktype = '0123456789abcde';
coaddopt.coaddtype = 1; % per-rx
switch daughter
  case 'a'
    coaddopt.temporaljack_splitdate = '20120719';
  case 'b'
    coaddopt.temporaljack_splitdate = '20130710';
  case 'd'
    coaddopt.temporaljack_splitdate = '20140717';
  case 'e'
    coaddopt.temporaljack_splitdate = '20150708';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Make pairmaps

if makepairmaps
    
  [sigmapfilename mapopt.deproj_map] = get_inputdeprojmap(year,bmsimtype);
  simopt.mapopt = mapopt;
  
  % Choose your beam map!
  simopt.beammapfilename = get_beammap(year,bmsimtype,bmsize);
  
  if ~exist(['simrunfiles/' num2str(nbase) '_' daughter '_dithers.mat'],'file')
    par_make_simrunfiles(tags,simopt,nbase,daughter);
  end

  mkdir(['pairmaps/' num2str(nbase)]);
  system(['ln -s /n/panlfs2/bicep/keck/pipeline/pairmaps/1351/real pairmaps/' ...
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
function [nbase daughter]  = get_simnum(year,bmsimtype,bmsize)

% Choose sernum based on type/size
switch bmsimtype
  case 'standard'
    switch bmsize
      case 1.2
        nbase = 3630;
      case 2
        nbase = 3631;
      case 4
        nbase = 3632;
      case 6
        nbase = 3633;
      case 8
        nbase = 3634;
    end
  case 'split'
    switch bmsize
      case 1.2
        nbase = 3640;
      case 2
        nbase = 3641;
      case 4
        nbase = 3642;
      case 6
        nbase = 3643;
      case 8
        nbase = 3644;
    end
  case 'floor' % Constructed Gaussian maps from measured beamparams
    switch bmsize
      case 1.2
        nbase = 3650;
      case 2
        nbase = 3651;
      case 4
        nbase = 3652;
      case 6
        nbase = 3653;
      case 8
        nbase = 3654;
    end
end

% Choose daughter based on year
switch year
  case 2012
    daughter = 'a';
  case 2013
    daughter = 'b';
  case 2014
    daughter = 'd';
  case 2015
    daughter = 'e';
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [sigmap deproj_map] = get_inputdeprojmap(year,bmsimtype)

% Input Planck map -- no beam smoothing applied
switch year
  case {2012,2013}
    sigmap = ...
        {'input_maps/planck/planck_derivs_nopix/HFI_SkyMap_100_2048_R1.10_nominal_nside0512_nlmax1280_bnone.fits'};
  case 2014
    sigmap = ...
        {'input_maps/planck/planck_derivs_nopix/HFI_SkyMap_100_2048_R1.10_nominal_nside0512_nlmax1280_bnone.fits',...
         'input_maps/planck/planck_derivs_nopix/HFI_SkyMap_143_2048_R1.10_nominal_nside0512_nlmax1280_bnone.fits'};
  case 2015
    sigmap = ...
        {'input_maps/planck/planck_derivs_nopix/HFI_SkyMap_100_2048_R1.10_nominal_nside0512_nlmax1280_bnone.fits',...
         'input_maps/planck/planck_derivs_nopix/HFI_SkyMap_143_2048_R1.10_nominal_nside0512_nlmax1280_bnone.fits',...
         'input_maps/planck/planck_derivs_nopix/HFI_SkyMap_217_2048_R1.10_nominal_nside0512_nlmax1280_bnone.fits'};
end

% Deprojection map -- for regular sims, this is the beam profile-smoothed
% Planck map.  For constructed Gaussian beam maps, we use
% Gaussian-smoothed Planck maps instead (so deprojection should work
% perfectly-ish)
switch year
  case {2012,2013}
    switch bmsimtype
      case {'standard','split'}
        deproj_map = ...
            {'input_maps/planck/planck_derivs_nopix/synfast_deproj_143_nominal_B2.fits'};
      case 'floor'
        deproj_map = ...
            {'input_maps/planck/planck_derivs_nopix/HFI_SkyMap_143_2048_R1.10_nominal_nside0512_nlmax1280_b30.52.fits'};
    end
  case 2014
    switch bmsimtype
      case {'standard','split'}
        deproj_map = ...
            {'input_maps/planck/planck_derivs_nopix/HFI_SkyMap_100_2048_R1.10_nominal_nside0512_nlmax1280_bKuber100.fits';...
             'input_maps/planck/planck_derivs_nopix/synfast_deproj_143_nominal_B2.fits'};
      case 'floor'
        deproj_map = ...
            {'input_maps/planck/planck_derivs_nopix/HFI_SkyMap_100_2048_R1.10_nominal_nside0512_nlmax1280_b43.23.fits';...
             'input_maps/planck/planck_derivs_nopix/HFI_SkyMap_143_2048_R1.10_nominal_nside0512_nlmax1280_b30.52.fits'};
    end
  case 2015
    switch bmsimtype
      case {'standard','split'}
        deproj_map = ...
            {'input_maps/planck/planck_derivs_nopix/HFI_SkyMap_100_2048_R1.10_nominal_nside0512_nlmax1280_bKuber100.fits';...
             'input_maps/planck/planck_derivs_nopix/synfast_deproj_143_nominal_B2.fits';...
             'input_maps/planck/planck_derivs_nopix/HFI_SkyMap_217_2048_R1.10_nominal_nside0512_nlmax1280_bKuber220.fits'};    
      case 'floor'
        deproj_map = ...
            {'input_maps/planck/planck_derivs_nopix/HFI_SkyMap_100_2048_R1.10_nominal_nside0512_nlmax1280_b43.23.fits';...
             'input_maps/planck/planck_derivs_nopix/HFI_SkyMap_143_2048_R1.10_nominal_nside0512_nlmax1280_b30.52.fits';...
             'input_maps/planck/planck_derivs_nopix/HFI_SkyMap_217_2048_R1.10_nominal_nside0512_nlmax1280_b20.06.fits'};    
     end
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function filename = get_beammap(year,bmsimtype,bmsize)

sizestr = [strrep(num2str(bmsize),'.','p') 'deg'];

switch bmsimtype
  case 'standard' % Use xyear composite maps
    filename = ['beammaps/maps_composite/ffbm_' num2str(year) ...
                '_all_xyear_' sizestr '_dk0.mat'];
  case 'split'
    filename = ['beammaps/maps_split/ffbm_' num2str(year) ...
                '_half_xyear_' sizestr '_dk0_diff.mat'];
  case 'floor'
    filename = ['beammaps/maps_constructed/k' num2str(year) ...
                '_obs_' sizestr '.mat'];
end

return

