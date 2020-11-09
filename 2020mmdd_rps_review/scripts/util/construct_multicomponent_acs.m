function [ac,coaddopt,m]=construct_multicomponent_acs(expt,rlz,p,snum,planck_jacktype)
% [ac,coaddopt,m]=construct_multicomponent_acs(expt,rlz,p,snum,planck_jacktype)
%
% Outputs simulated Planck or B2 ac structure in uK_CMB containing lensed LCDM, CMB tensors,
% dust, and statistical noise:
% 
% Outputs:
%
% ac,coaddopt,m    - cell array of acs; coaddopt, m structure as usual
%
% Inputs:
%
% expt -  determines what noise realizations are used and frequency scaling of dust map
%         'b150'
%         'p100'
%         'p143'
%         'p217'
%         'p353'
% rlz - realization number, one for everything
% p   - model parameters as in Harvard multicomponent analysis
% snum- sim number to use for saving (if not specified, don't save)
% planck_jacktype - to select Planck det set, year, or half-ring split maps
%
%  The model parameters [p] are:
%
%    1. r, tensor-to-scalar ratio
%    2. A_L, lensing amplitude
%    3. P_sync, polarized synchrotron amplitude, in uK^2_CMB, at 150 GHz and ell=80.
%    4. P_dust, polarized dust amplitude, in uK^2_CMB, at 353 GHz and ell=80.
%    5. beta_sync, polarized synchrotron frequency spectral index
%    6. gamma_sync, polarized synchotron spatial spectral index
%    7. beta_dust, polarized dust frequency spectral index
%    8. gamma_dust, polarized dust spatial spectral index
%    9. E/B ratio for sync (value of 1 implies equal power in E and B)
%   10. E/B ratio for dust (value of 1 implies equal power in E and B)
%   11. epsilon, sync/dust spatial correlation parameter
%   12. T_greybody, dust greybody temperature (optional)
%   13. scaling for PSM foreground realizations from Jamie. T
%   14. scaling for old PSM foreground realizations (as in PRL) from Jamie T.
%

if ~exist('snum','var') || isempty(snum)
  snum=[];
end
if ~exist('planck_jacktype','var') || isempty(planck_jacktype)
  planck_jacktype=0;
end
if ~ischar(planck_jacktype)
  planck_jacktype=num2str(planck_jacktype);
end

% Maybe we've got a vector of rlz...
if length(rlz)>1
  for i=1:length(rlz)
    construct_multicomponent_acs(expt,rlz(i),p,snum,planck_jacktype);
  end
  return
end

% Maybe we're being asked for all Planck bands
if strcmpi(expt,'planck')
  planck_bands=[30 44 70 100 143 217 353];
  clear map
  for i=1:length(planck_bands)
    tmpexpt=['p' num2str(planck_bands(i),'%.3d')];
    [tmpac,coaddopt,m]=construct_multicomponent_acs(tmpexpt,rlz,p,[],planck_jacktype);
    if iscell(tmpac)
      for j=1:length(tmpac)
        ac{j}(i,:)=tmpac{j};
      end
    else
      ac(i,:)=tmpac;
    end
  end
  if ~isempty(snum)
    outpath=fullfile('maps',num2str(snum,'%.4d'));
    if ~exist(outpath,'dir')
      mkdir(outpath);
    end
    fname=[num2str(rlz,'%.3d') '8_a_filtp3_weight3_gs_dp1100_jack' planck_jacktype '.mat'];
    save(fullfile(outpath,fname),'ac','m','coaddopt');
  end
  return
end

% Scaling for PSM dust (if including at all)
use_old_psm=false;
psmamp=0;
if length(p)>=13  && isfinite(p(13)) && p(13)>0
  psmamp=p(13);
end
% Scaling for old PSM dust (if including at all)
if length(p)>=14  && isfinite(p(14)) && p(14)>0
  if psmamp~=0
    error(['Trying to specify both old and new PSM versions in same multicomponent map!']);
  end
  psmamp=p(14);
  use_old_psm=true;
end

rlzs=rlz;
rlzn=rlz;
rlzd=rlz;
rlzp=rlz;
rlzy=rlz;

obs_expt=get_experiment_name();

% Get maps
[ac_lcdm,coaddopt,m] = get_lcdm_ac(obs_expt,rlzs); % lensed LCDM
ac_r    = get_bmode_ac(obs_expt,rlzs); % r=0.1 tensors
if psmamp~=0
  if use_old_psm
    ac_d    = get_oldpsm_ac(obs_expt,expt,rlzd,planck_jacktype); % Old-PSM foreground realizations from Jamie T.
  else
    ac_d    = get_psm_ac(obs_expt,expt,rlzd,planck_jacktype); % PSM dust realizations from Jamie T
  end
end

% Power law maps: note Planck Paper XXX uses spatial index alpha=gamma-2
if isfinite(p(8))
  ac_p    = get_fgpowerlaw_ac(obs_expt,rlzp,'power-law dust',p(8)-2,p(10)); % Power-law dust foreground
else
  ac_p    = get_oldpsm_ac(obs_expt,'p353',rlzd,planck_jacktype);
end
ac_y    = get_fgpowerlaw_ac(obs_expt,rlzy,'power-law sync',p(6)-2,p(9)); % Power-law synch foreground

% Need fake coaddopt with weight=0 to pass to cal_coadd_ac
% This was we rescale the signal, but not any of the weights
co=[];
co.weight=0;

% Scale tensor map.  Scaling by r works in bandpower, i.e. map^2, so take sqrt.
ac_r = cal_coadd_ac(ac_r,sqrt(p(1)/0.1),co);

% If dust and sync maps are compatible, and correlated, correlate them now
if (p(11)~=0)
  if (p(9)~=p(10))
    error(['Incompatible parameters: dust and sync are to be correlated, but E/B ratio does not match.']);
  end
  if (p(6)~=p(8))
    error(['Incompatible parameters: dust and sync are to be correlated, but spatial spectral index does not match.']);
  end
  ac_y=dust_sync_mix_ac(ac_p,ac_y,p(11));
end

% Get bandpass for experiment 'expt'
% This uses Colin's code for loading bandpasses
bandpass=get_freq_band(expt);
bandpass=bandpass{1};

% For power-law sync and dust:
% Calculate frequency scaling using Colin's code
if length(p)>=12 && isfinite(p(12)) && p(12)>0
  tgraydust=p(12);
else
  tgraydust=[];
end
% fsc=freq_scaling(bandpass, beta, temp, nu0)
fsc_p=freq_scaling(bandpass,p(7),tgraydust,353);  % Dust scales by 353 intensity now!!!
fsc_y=freq_scaling(bandpass,p(5),[],150);

% Careful with amplitude scaling of foreground components...
% the amplitude parameters p(3), p(4) are defined in bandpower, so have a sqrt()
% the frequency scalings fsc_* are in map, so no sqrt().

% Dust amplitude defined for D_l^BB at *353* GHz and ell=100
ac_p = cal_coadd_ac(ac_p,fsc_p*sqrt(p(4)),co);
% Synch amplitude defined for D_l^BB at 150 GHz and ell=100
ac_y = cal_coadd_ac(ac_y,fsc_y*sqrt(p(3)),co);

% PSM foreground maps are already done per-frequency.  Don't apply any additional
% scaling with band for Planck bands.  Just turn on or off with psmamp, which is 1 or 0.
% For BICEP2/Keck, scale from 143 to real band.  This is a little ugly...
if lower(expt(1))=='p'
  fsc_d=1;
else
  fsc_d=freq_scaling(bandpass,p(7),tgraydust,143);
end
if psmamp~=0
  ac_d = cal_coadd_ac(ac_d,fsc_d*sqrt(psmamp),co);
end

% Get noise maps and add to signal map (add systematics map here if necessary)
ac_n = get_noise_ac(obs_expt,expt,rlzn,planck_jacktype);

% Put different component ac's into cells of a cell array
ac{1} = ac_n;
ac{2} = ac_lcdm;
ac{3} = ac_r;
ac{4} = add_acs(ac_p,ac_y);
if psmamp~=0
  ac{4} = add_acs(ac{4},ac_d);
end
% ac=stripacs(ac);

% Expand all to have same size
ac = expandac(ac);

if ~isempty(snum)
  outpath=fullfile('maps',num2str(snum,'%.4d'));
  if ~exist(outpath,'dir')
    mkdir(outpath);
  end
  fname=[num2str(rlz,'%.3d') '8_a_filtp3_weight3_gs_dp1100_jack' planck_jacktype '.mat'];
  save(fullfile(outpath,fname),'ac','m','coaddopt');
end

return


%%%%%%%%%%%%%%%%%%%%%%%%
% Expand all ac structures to have same size
function ac=expandac(ac)

if ~iscell(ac)
  return
end

sz=[0 0];
for i=1:length(ac)
  tmpsz=size(ac{i});
  sz(1)=max([tmpsz(1),sz(1)]);
  sz(2)=max([tmpsz(2),sz(2)]);
end
for i=1:length(ac)
  tmpsz=size(ac{i});
  if any(tmpsz<sz)
    ac{i}=repmat(ac{i},sz(1)/tmpsz(1),sz(2)/tmpsz(2));
  end
end

return

%%%%%%%%%%%%%%%%%%%%%%%%
% Note this is different from the standard pipeline function addac!
% We want to add multiple signal components, not coadd different observations!
function ac=add_acs(ac1,ac2,varargin)

if nargin > 2
  ac=add_acs(add_acs(ac1,ac2),varargin{:});
  return
end

% Expand as needed
n=max([size(ac1,2) size(ac2,2)]);
if n>1
  if size(ac1,2)==1
    ac1=repmat(ac1,1,n);
  end
  if size(ac2,2)==1
    ac2=repmat(ac2,1,n);
  end
end
% Assume weights, etc. are the same
ac=ac1;
sumfields={'wz','wcz','wsz'};
for i=1:size(ac,1)
  for j=1:size(ac,2)
    for k=1:length(sumfields)
      ac(i,j).(sumfields{k})=ac(i,j).(sumfields{k})+ac2(i,j).(sumfields{k});
    end
  end
end

return


%%%%%%%%%%%%%%%%%%%%%%%%%%
function ac=get_noise_ac(obs_expt,expt,rlz,planck_jacktype)

switch expt(1)
 
 case 'b'
  if ~strcmp(planck_jacktype,'0')
    error(['BICEP2 or Keck maps should be made with planck_jacktype=0.  Asked for planck_jacktype=' planck_jacktype]);
  end
  switch(obs_expt)
    case 'bicep2',
      load(sprintf('/n/bicepfs2/bicep2/pipeline/maps/0751/%03d6_a_filtp3_weight3_gs_dp1100_jack0.mat',rlz));
      ac=cal_ac_wrapper(ac,m,coaddopt);
    case 'keck',
      load(sprintf('/n/panlfs2/bicep/keck/pipeline/maps/1353/%03d6_ab_filtp3_weight3_gs_dp1100_jack0.mat',rlz));
      ac=cal_ac_wrapper(ac,m,coaddopt);
    otherwise, error(['Unknown experiment ' obs_expt']); 
  end
 case 'p'
  switch(obs_expt)
    case 'bicep2',
      load(sprintf('/n/bicepfs2/b2planck/pipeline/maps/1613/%03d6_a_filtp3_weight3_gs_dp1100_jack%s.mat',rlz,planck_jacktype));
      fac=[30,44,70,100,143,217,353];
      ind=find(fac==str2num(expt(2:end)));
      ac=cal_ac_wrapper(ac(ind,:),m,coaddopt);
    case 'keck',
      load(sprintf('/n/bicepfs2/b2planck/pipeline/maps/1614/%03d6_ab_filtp3_weight3_gs_dp1100_jack%s.mat',rlz,planck_jacktype));
      fac=[30,44,70,100,143,217,353];
      ind=find(fac==str2num(expt(2:end)));
      ac=cal_ac_wrapper(ac(ind,:),m,coaddopt);
  end
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function bandpass=get_freq_band(expt)

expt=upper(expt);
switch(expt)
  case 'B150', expt='B2_150';
end
lo=[];
lo.expt(1).name=expt;
lo = like_read_bandpass(lo);
bandpass = lo.bandpass;

return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ac,coaddopt,m]=get_lcdm_ac(expt,rlz)

switch(expt)
  case 'bicep2',
load(sprintf('/n/bicepfs2/bicep2/pipeline/maps/0751/%03d5_a_filtp3_weight3_gs_dp1100_jack0.mat',rlz));
[ac,coaddopt]=cal_ac_wrapper(ac,m,coaddopt);
  case 'keck',
load(sprintf('/n/panlfs2/bicep/keck/pipeline/maps/1353/%03d5_ab_filtp3_weight3_gs_dp1100_jack0.mat',rlz));
[ac,coaddopt]=cal_ac_wrapper(ac,m,coaddopt);
end

return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ac=get_bmode_ac(expt,rlz)

switch(expt)
  case 'bicep2',
load(sprintf('/n/bicepfs2/bicep2/pipeline/maps/0751/%03d4_a_filtp3_weight3_gs_dp1100_jack0.mat',rlz));
ac=cal_ac_wrapper(ac,m,coaddopt);
  case 'keck',
load(sprintf('/n/panlfs2/bicep/keck/pipeline/maps/1353/%03d4_ab_filtp3_weight3_gs_dp1100_jack0.mat',rlz));
ac=cal_ac_wrapper(ac,m,coaddopt);
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ac]=get_fgpowerlaw_ac(expt,rlz,seedstr,alpha_BB,BB_to_EE)

fname1=sprintf('%03d8_a_filtp3_weight3_gs_dp1100_jack0.mat',rlz);
fname2=sprintf('%03d2_a_filtp3_weight3_gs_dp1100_jack0.mat',rlz);
for snum=0560:0569
  if exist(fullfile('maps',num2str(snum,'%.4d'),fname1),'file')
    fname=fname1;
  else
    fname=fname2;
  end
  if exist(fullfile('maps',num2str(snum,'%.4d'),fname),'file')
    disp(['Found file ' fullfile('maps',num2str(snum,'%.4d'),fname)]);
    tmp=load(fullfile('maps',num2str(snum,'%.4d'),fname),'coaddopt');
    co=tmp.coaddopt;
    if ~isfield(co,'fgopt')
      continue
    end
    if strcmp(seedstr,co.fgopt.seedhash)
      disp(['In ' fname ' : found matching tag string ' seedstr]);
    else
      continue
    end
    if alpha_BB ~= co.fgopt.alpha_BB
      disp(['In ' fname ' : found alpha_BB=' num2str(co.fgopt.alpha_BB) ', wanted ' num2str(alpha_BB)]);
      continue
    end
    if BB_to_EE ~= co.fgopt.BB_to_EE
      disp(['In ' fname ' : found BB_to_EE=' num2str(co.fgopt.alpha_BB) ', wanted ' num2str(alpha_BB)]);
      continue
    end
    load(fullfile('maps',num2str(snum,'%.4d'),fname));
    % Foreground maps could be an ac structure that need make_map and abscal;
    % or a map structure that is already abscalled
    if exist('ac','var') || ~exist('map','var')
      if ~isfield(coaddopt,'coaddtype') || isempty(coaddopt.coaddtype)
        coaddopt.coaddtype=0;
      end
      ac=cal_ac_wrapper(ac,m,coaddopt);
    end
    return
  end
end

error(['Suitable power law foreground maps not found!']);

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ac,coaddopt,m]=get_psm_ac(obs_expt,expt,rlz,planck_jacktype)

switch expt(1)
 case 'b'
  freq=143;
  if ~strcmp(planck_jacktype,'0')
    error(['BICEP2 or Keck maps should be made with planck_jacktype=0.  Asked for planck_jacktype=' planck_jacktype]);
  end
  switch(obs_expt)
    case 'bicep2',
      load(sprintf('/n/panlfs2/bicep/jetolan/bicep2/maps/0704/xxx8_allcmb_ffp8_fg_bpm_%03d_0%03d_filtp3_weight3_gs_dp1100_jack0.mat',freq,rlz));
      ac=cal_ac_wrapper(ac,m,coaddopt);
    case 'keck',
      load(sprintf('/n/panlfs2/bicep/jetolan/keck/maps/0706/xxx8_cmb_subset2012_2013_rxall_ffp8_fg_bpm_%03d_0%03d_filtp3_weight3_gs_dp1100_jack0.mat',freq,rlz));
      ac=cal_ac_wrapper(ac,m,coaddopt);
    otherwise, error(['Unknown experiment ' obs_expt']);
  end
 case 'p'
  freq=str2num(expt(2:end));
  switch(obs_expt)
    case 'bicep2',
      load(sprintf('/n/panlfs2/bicep/jetolan/bicep2/maps/0704/xxx8_allcmb_ffp8_fg_bpm_%03d_0%03d_filtp3_weight3_gs_dp1100_jack0.mat',freq,rlz));
      ac=cal_ac_wrapper(ac,m,coaddopt);
    case 'keck',
      load(sprintf('/n/panlfs2/bicep/jetolan/keck/maps/0706/xxx8_cmb_subset2012_2013_rxall_ffp8_fg_bpm_%03d_0%03d_filtp3_weight3_gs_dp1100_jack0.mat',freq,rlz));
      ac=cal_ac_wrapper(ac,m,coaddopt);
  end
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ac,coaddopt,m]=get_oldpsm_ac(obs_expt,expt,rlz,planck_jacktype)

switch expt(1)
 case 'b'
  freq=143;
  if ~strcmp(planck_jacktype,'0')
    error(['BICEP2 or Keck maps should be made with planck_jacktype=0.  Asked for planck_jacktype=' planck_jacktype]);
  end
  switch(obs_expt)
    case 'bicep2',
      load(sprintf('/n/panlfs2/bicep/jetolan/bicep2/maps/0704_psm178/xxx8_allcmb_psm178_fg_bpm_%03d_0%03d_filtp3_weight3_gs_dp1100_jack0.mat',freq,rlz));
      ac=cal_ac_wrapper(ac,m,coaddopt);
    case 'keck',
      load(sprintf('/n/panlfs2/bicep/jetolan/keck/maps/0706_psm178/xxx8_cmb_subset2012_2013_rxall_psm178_fg_bpm_%03d_0%03d_filtp3_weight3_gs_dp1100_jack0.mat',freq,rlz));
      ac=cal_ac_wrapper(ac,m,coaddopt);
    otherwise, error(['Unknown experiment ' obs_expt']);
  end
 case 'p'
  freq=str2num(expt(2:end));
  switch(obs_expt)
    case 'bicep2',
      load(sprintf('/n/panlfs2/bicep/jetolan/bicep2/maps/0704_psm178/xxx8_allcmb_psm178_fg_bpm_%03d_0%03d_filtp3_weight3_gs_dp1100_jack0.mat',freq,rlz));
      ac=cal_ac_wrapper(ac,m,coaddopt);
    case 'keck',
      load(sprintf('/n/panlfs2/bicep/jetolan/keck/maps/0706_psm178/xxx8_cmb_subset2012_2013_rxall_psm178_fg_bpm_%03d_0%03d_filtp3_weight3_gs_dp1100_jack0.mat',freq,rlz));
      ac=cal_ac_wrapper(ac,m,coaddopt);
  end
end

return



%%%%%%%%%%%%%%%%%%%%%%%%%%5
function d=stripmap(b)

for k=1:numel(b)

  a=b(k);
  
  c.x_tic=a.x_tic;
  c.y_tic=a.y_tic;
  c.T=a.T;
  c.Tvar=a.Tvar;
  c.Titime=a.Titime;
  c.Q=a.Q;
  c.U=a.U;
  c.Qvar=a.Qvar;
  c.Uvar=a.Uvar;
  c.QUcovar=a.QUcovar;
  c.Pitime=a.Pitime;

  d(k)=c;

end

return


%%%%%%%%%%%%%%%%%%%%%%%%%%5
function ac_y=dust_sync_mix_ac(ac_p,ac_y,epsilon)

fn=fieldnames(ac_y);
for i=1:length(fn)
  if ismember(fn{i},{'x_tic','y_tic'})
    continue
  end
  ac_y.(fn{i})=ac_y.(fn{i})*(1-epsilon)+ac_p.(fn{i})*(epsilon); % correct scaling?
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%5
function [ac,coaddopt]=cal_ac_wrapper(ac,m,coaddopt)

calfac=[];
if isfield(coaddopt,'mapopt') && ~isempty(coaddopt.mapopt)
  if isfield(coaddopt.mapopt{1},'simopt')
    so=coaddopt.mapopt{1}.simopt;
    calfac=so.ukpervolt(1);
  end
end
if isempty(calfac)
  calfac=get_ukpervolt();
  calfac=calfac(1);
end
if isfield(coaddopt,'ukpv_applied') && ~isempty(coaddopt.ukpv_applied)
  % Trust that existence of coaddopt.ukpv_applied means map is truly
  % calibrated in CMB units
  % calfac=calfac/coaddopt.ukpv_applied(1);
  calfac=1;
else
  coaddopt.ukpv_applied=calfac;
end
co=[];
co.weight=3;
ac=cal_coadd_ac(ac,calfac,co);

return

