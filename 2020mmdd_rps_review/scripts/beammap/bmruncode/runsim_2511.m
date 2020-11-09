function runsim_2511(doreal,dosim,dorealcoadd,dosimcoadd,doaps,...
    rlz,submit,daughter)
% runsim_2511(doreal,dosim,dorealcoadd,dosimcoadd,doaps,...
%    rlz,submit,daughter)
%
% Sim 2511: Beam map sim for Keck 2014
%   integral normalized, only common AB
% Following the format of runsim_0751 - much of this is shamelessly copied,
% so lots of unnecessary stuff
% Uses composite beam maps from 2014 beam mapping campaign
%
% One realization on the real sky: 25110001 
% 10 realizations on sim skies:    2511xxx2 (need different "types" because
% of how sigmapfilename is read)
%
% To make real sky: 
% >> doreal = 1; rlz = 0; 
% To make sims:
% >> dosim = 0; rlz = 1:10; 
%
% Daughter simsets
% 'a': coaddopt.coaddtype = 1, apsopt.coaddrx = 2;
% 'b': coaddopt.coaddytpe = 2; apsopt.coaddrx = 0;

nbase = 2511;

% Defaults: do nothing at all
if(~exist('doreal','var'))
  doreal = 0;
end
if(~exist('dosim','var'))
  dosim = 0;
end
if(~exist('dorealcoadd','var'))
  dorealcoadd = 0;
end
if(~exist('dosimcoadd','var'))
  dosimcoadd = 0;
end
if(~exist('doaps','var'))
  doaps = 0;
end
if(~exist('rlz','var'))
  rlz = 0;
end
if(~exist('submit','var'))
  submit = true;
end
if(~exist('daughter','var'))
  daughter = 'a';
end

% Memory usage
memreal = 15000;
% Extremely large amounts of memory for per-pair maps
switch daughter
  case {'a'}
    memcoadd = 30000;
    memcoaddcoadd = 12000;
    memcoaddcoadddp = 12000;
  case 'b'
    memcoadd = 120000;
    memcoaddcoadd = 120000;
    memcoaddcoadddp = 120000;
end
maxtime = 400;
maxtime_real = 400;
if dorealcoadd | dosimcoadd 
  maxtime = 12*60;
end
if doreal | dosim
  maxtime = 190;
end

% Pick appropriate queues
display('Using SLURM')
slurm_queue = 'general,serial_requeue,itc_cluster';
queue_sim = slurm_queue;
queue_coadd = slurm_queue;
queue_coadd_coadd = slurm_queue;
queue_other = slurm_queue;

% Global options
om = 0; % onlymissing option (0 = remake all)
ntpg = 4; % NTagsPerJob

% Sim options
if strcmp(get_experiment_name(),'keck')
  rlzc = 3; % rlzchunksize for runsim
end
tagc = 1; % tagchunksize for runsim

% Coadd options
js = 0; % option for FarmJacksSeparately for coadd
farmdps = 1;

% Which deprojection we want to do:
switch daughter
  case {'a'}
    deprojs = [0,0,0,0;  % No deprojection
      1,0,0,0;  % relgain
      0,1,0,0;  % A/B
      1,1,0,0;  % relgain + A/B offsets   
      1,1,1,0;  % beamwidth + A/B offsets + relgain
      1,1,0,1;  % ellipticity + A/B offsets + relgain
      1,1,1,1;  % all
      1,1,0,2]; % ellipticity subtraction
  case 'b' % per-pair maps take FOREVER so only do a few
    deprojs = [0,0,0,0;
      1,1,0,0;
      1,1,0,1;
      1,1,1,1;
      1,1,0,2];
end
rlz

% WHAT TAGS DO I WANT
realfile = 'maps/1351/real_d_filtp3_weight3_gs_dp1100_jack01.mat';
x = load(realfile);
tags = x.coaddopt.tags;
[tagsublist,mtl] = get_tag_sublist(x.coaddopt);
clear x;

clear simopt
simopt.sernum = [nbase 'xxx1'];
simopt.rlz = 0;
simopt.siginterp = 'linear'; % for bm sims
simopt.beamcen = 'obs';
simopt.diffpoint = 'ideal'; % Since the beam maps include this
simopt.chi = 'obs';
simopt.epsilon = 'obs';

% From my stuff:
simopt.noise = 'none';
simopt.sig = 'nopol';
simopt.sigmaptype = 'healmap';
simopt.maketod = false;
simopt.interpix = 0.1;

ukpervolt = get_ukpervolt('2014');

simopt.ukpervolt = ukpervolt;
simopt.coord = 'C';
%simopt.force_ab_int_from_common_pix = false; % What is this?
simopt.update = 0; % remake pairmaps which already exist (1 only missing)

clear mapopt

mapopt.beamcen = 'obs';
mapopt.chi = 'obs';
mapopt.epsilon = 'obs';
mapopt.acpack = 0;
mapopt.gs = 1;
mapopt.filt = 'p3';
mapopt.deproj = [1:6,9,10]; % 9/10 for residbeam
mapopt.update = true;

% My stuff
% Also made a symlink in pairmaps/2511/real 
mapopt.realpairmapset = 'pairmaps/1351/real';

clear coaddopt

chflags = get_default_chflags([],'2014');

coaddopt.chflags = chflags;
coaddopt.daughter = daughter;
coaddopt.deproj_timescale = 'by_phase';
coaddopt.filt = 'p3';
coaddopt.gs = 1;
coaddopt.jacktype = '0123456789abcde';

switch daughter
  case 'a'
    coaddopt.coaddtype = 1;
  case 'b'
    coaddopt.coaddtype = 2;
end

coaddopt.temporaljack_splitdate = '20140717';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Let's work (?!)

if doreal 
  % We call it "real" but need to treat it like a sim - only difference is input maps
  
  simopt.mapopt = mapopt;
  simopt.beammapfilename = ...
      'beammaps/maps_composite/ffbm_2014_all_onlycommonab.mat';
  if ~exist(['simrunfiles/' num2str(nbase) '_' daughter '_dithers.mat'],'file')
    par_make_simrunfiles(tags,simopt,nbase,daughter);
  end
  %{
  % Maps used for rev 1 (see posting 2015-06-09: Keck 2014 Beam Map Simulation)
  sigmapfilename = {...
      'input_maps/planck/planck_maps/HFI_SkyMap_100_2048_R1.10_nominal_l600_00p00_cel_uk.fits',...
  'input_maps/planck/planck_maps/HFI_SkyMap_143_2048_R1.10_nominal_l600_00p00_cel_uk.fits'};
  %}
  sigmapfilename = {...
      'input_maps/planck/planck_derivs_nopix/HFI_SkyMap_100_2048_R1.10_nominal_nside0512_nlmax1280_bnone.fits',...
  'input_maps/planck/planck_derivs_nopix/HFI_SkyMap_143_2048_R1.10_nominal_nside0512_nlmax1280_bnone.fits'};

  mapopt.deproj_map = {...
      'input_maps/planck/planck_derivs_nopix/HFI_SkyMap_100_2048_R1.10_nominal_nside0512_nlmax1280_bKuber100.fits';...
  'input_maps/planck/planck_derivs_nopix/synfast_deproj_143_nominal_B2.fits'};

  % Residbeam stuff
  mapopt.residbeam.beammap = ...
      'beammaps/maps_composite/ffbm_2014_all_onlycommonab.mat';
  mapopt.residbeam.skymap = ...
      {'input_maps/planck/planck_derivs_nopix/HFI_SkyMap_100_2048_R1.10_nominal_nside0512_nlmax1280_bnone.fits';...
  'input_maps/planck/planck_derivs_nopix/HFI_SkyMap_143_2048_R1.10_nominal_nside0512_nlmax1280_bnone.fits'};

  type = 1;
  runsim(nbase,daughter,mapopt,'nopol','none',type,sigmapfilename,rlz,mtl,...
      true,rlzc,om,800,0,[],1,queue_sim,0,[],tagc,0,maxtime,submit);
 
end

if dosim 

  simopt.mapopt = mapopt;
  simopt.beammapfilename = 'beammaps/maps_composite/ffbm_2014_all.mat';
  if ~exist(['simrunfiles/' num2str(nbase) '_' daughter '_dithers.mat'],'file')
    par_make_simrunfiles(tags,simopt,nbase,daughter);
  end
  
  mapopt.deproj_map = {...
      'input_maps/camb_planck2013_r0/map_unlens_n0512_rxxxx_sKuber100_cPl100_dPl100.fits';...
      'input_maps/camb_planck2013_r0/map_unlens_n0512_rxxxx_sB2bbns_cPl143_dPl143.fits'};

  % Signal sim
  type = 2;
  sigmapfilename = {...
      'input_maps/camb_planck2013_r0/map_unlens_n2048_rxxxx_sKuber100_cPl100_dNoNoi.fits',...
      'input_maps/camb_planck2013_r0/map_unlens_n2048_rxxxx_sB2bbns_cPl143_dNoNoi.fits'};      
  
  runsim(nbase,daughter,mapopt,'nopol','none',type,sigmapfilename,rlz,mtl,...
      true,rlzc,om,800,0,[],1,queue_sim,0,[],tagc,0,maxtime,submit);
  
end % dosim

if dorealcoadd
  
  for k = 1
    k
    
    coaddopt.sernum = sprintf('%04dxxx%01d',nbase,k);
    coaddopt.save_cuts = false;
    coaddopt.tagsublist = tagsublist;

    if farmdps
      for jj = 1:size(deprojs,1)
	coaddopt.deproj = deprojs(jj,:);
	farm_coaddpairmaps(tags,coaddopt,'OnlyMissing',om,'Realizations',rlz,...
	    'FarmJobs',1,'JobLimit',800,'FarmJacksSeparately',js,...
	    'Queue',queue_coadd_coadd,'MemRequire',memcoaddcoadddp,...
	    'SplitSubmission',20,'UseCompiled',0,...
	    'maxtime',maxtime,'submit',submit);
      end
    else
      coaddopt.deproj = deprojs;
      farm_coaddpairmaps(tags,coaddopt,'OnlyMissing',om,'Realizations',rlz,...
	  'FarmJobs',1,'JobLimit',800,'FarmJacksSeparately',js,...
	  'Queue',queue_coadd_coadd,'MemRequire',memcoaddcoadddp,...
	  'SplitSubmission',20,'UseCompiled',0,...
	  'maxtime',maxtime,'submit',submit);
    end
  end % simtypes
  
end % dorealcoadd

if dosimcoadd
  
  for k = 2
    k
    
    coaddopt.sernum = sprintf('%04dxxx%01d',nbase,k);
    coaddopt.save_cuts = false;
    coaddopt.tagsublist = tagsublist;

    if farmdps
      for jj = 1:size(deprojs,1)
	coaddopt.deproj = deprojs(jj,:);
	farm_coaddpairmaps(tags,coaddopt,'OnlyMissing',om,'Realizations',rlz,...
	    'FarmJobs',1,'JobLimit',800,'FarmJacksSeparately',js,...
	    'Queue',queue_coadd_coadd,'MemRequire',memcoaddcoadddp,...
	    'SplitSubmission',20,'UseCompiled',0,...
	    'maxtime',maxtime,'submit',submit);
      end
    else
      coaddopt.deproj = deprojs;
      farm_coaddpairmaps(tags,coaddopt,'OnlyMissing',om,'Realizations',rlz,...
	  'FarmJobs',1,'JobLimit',800,'FarmJacksSeparately',js,...
	  'Queue',queue_coadd_coadd,'MemRequire',memcoaddcoadddp,...
	  'SplitSubmission',20,'UseCompiled',0,...
	  'maxtime',maxtime,'submit',submit);
    end
  end % simtypes
  
end % dosimcoadd

if doaps
  
  if ~exist(['aps/' nbase],'dir')
    mkdir(['aps/' nbase]);
  end
  
  apsopt.save_coaddopts = 1;
  apsopt.pure_b = 'kendrick';
  apsopt.ukpervolt = simopt.ukpervolt;
  switch daughter
    case 'a'
      apsopt.coaddrx = 2; % 2: per-frequency
    case 'b'
      apsopt.coaddrx = 0;
  end
  
  d = dir(['maps/' num2str(nbase) '/*_' coaddopt.daughter '_*.mat']);

  for ii = 1:length(d)
    mapname = [num2str(nbase) '/' d(ii).name];
    reduc_makeaps(mapname,apsopt);
  end

  % If the regular sim, also make aps per-rx
  switch daughter
    case 'a'
      apsopt.coaddrx = 0;
      for ii = 1:length(d)
	mapname = [num2str(nbase) '/' d(ii).name];
	reduc_makeaps(mapname,apsopt);
      end
  end
  
end

return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Legacy sim code
    %d = dir('pairmaps/2510/0012/*.mat');
    %for ii = 1:length(d)
    %  tags2{ii} = d(ii).name(1:17);
    %end
    %tags = setdiff(mtl,tags2);
    %tags = mtl;
    %cmd = 'reduc_makesim(tags0,simopt)';
    %for i = 1:numel(tags)
    %  tags0 = tags(i);
    %  farmit('~/code/keckpipe/farmfiles/2510',cmd,'var',{'tags0','simopt'},...
    %     'queue','general,serial_requeue','mem',40000,'maxtime',maxtime);
    %end
    %subsim = true
