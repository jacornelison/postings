function [map,coaddopt,m,map_n,map_lcdm,map_r,map_d,map_s]=construct_multicomponent_maps(expt,rlz,p,snum,planck_jacktype)
% [map,coaddopt,m,map_n,map_lcdm,map_r,map_d]=construct_multicomponent_maps(expt,rlz,p,snum,planck_jacktype)
%
% Outputs simulated Planck or B2 map in uK_CMB containing lensed LCDM, CMB tensors,
% dust, and statistical noise:
% 
% Outputs:
%
% map,coaddopt,m   - all maps added together
% map_n            - noise map(s)
% map_lcdm         - lensed lcdm signal only map
% map_r            - tensors map
% map_d            - PSM dust map
% map_p            - power-law dust map
% map_y            - power-law synch map
% map_s            - total signal (map_lcdm+map_r+map_d)
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
    construct_multicomponent_maps(expt,rlz(i),p,snum,planck_jacktype);
  end
  return
end

% Maybe we're being asked for all Planck bands
if strcmpi(expt,'planck')
  planck_bands=[30 44 70 100 143 217 353];
  clear map
  for i=1:length(planck_bands)
    tmpexpt=['p' num2str(planck_bands(i),'%.3d')];
    [tmpmap,coaddopt,m,map_n,map_lcdm,map_r,map_d,map_s]=construct_multicomponent_maps(tmpexpt,rlz,p,[],planck_jacktype);
    map(i,:)=tmpmap;
  end
  if ~isempty(snum)
    outpath=fullfile('maps',num2str(snum,'%.4d'));
    if ~exist(outpath,'dir')
      mkdir(outpath);
    end
    fname=[num2str(rlz,'%.3d') '8_a_filtp3_weight3_gs_dp1100_jack' planck_jacktype '.mat'];
    save(fullfile(outpath,fname),'map','m','coaddopt');
  end
  return
end

rlzs=rlz;
rlzn=rlz;
rlzp=rlz;
rlzy=rlz;

obs_expt=get_experiment_name();

% Get maps
map_lcdm = get_lcdm_map(obs_expt,rlzs); % lensed LCDM
map_r    = get_bmode_map(obs_expt,rlzs); % r=0.1 tensors
[map_d,coaddopt,m] = get_dust_map(obs_expt); % PSM dust

% Power law maps: note Planck Paper XXX uses spatial index alpha=gamma-2
map_p    = get_fgpowerlaw_map(obs_expt,rlzp,'power-law dust',p(8)-2,p(10)); % Power-law dust foreground
map_y    = get_fgpowerlaw_map(obs_expt,rlzp,'power-law sync',p(6)-2,p(9)); % Power-law synch foreground

% Scale tensor map.  Scaling by r works in bandpower, i.e. map^2, so take sqrt.
map_r = cal_coadd_maps(map_r,sqrt(p(1)/0.1));

% If dust and sync maps are compatible, and correlated, correlate them now
if (p(11)~=0)
  if (p(9)~=p(10))
    error(['Incompatible parameters: dust and sync are to be correlated, but E/B ratio does not match.']);
  end
  if (p(6)~=p(8))
    error(['Incompatible parameters: dust and sync are to be correlated, but spatial spectral index does not match.']);
  end
  map_y=dust_sync_mix(map_p,map_y,p(11));
end

% Get bandpass for experiment 'expt'
% This uses Colin's code for loading bandpasses
bandpass=get_freq_band(expt);
bandpass=bandpass{1};

% For PSM dust only:
% Scale dust map using numbers from Ken Ganga
% http://b2p.planck.fr/index.php/AnalysisLogbook/2014-08-14DustConversion
f=str2num(expt(2:end));
fsc_d=1;
switch f
 case 100
  fsc_d=0.424885357233;
 case 143
  fsc_d=0.900287214698;
 case 150
  fsc_d=1;
 case 217
  fsc_d=3.12985190361;
 case 353
  fsc_d=25.1198188496;
end
d=0;
map_d = cal_coadd_maps(map_d,d*fsc_d);

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
map_p = cal_coadd_maps(map_p,fsc_p*sqrt(p(4)));
% Synch amplitude defined for D_l^BB at 150 GHz and ell=100
map_y = cal_coadd_maps(map_y,fsc_y*sqrt(p(3)));

% Add signal maps together
map_s = map_lcdm;
map_s = addmaps(map_s,map_r);
map_s = addmaps(map_s,map_d);
map_s = addmaps(map_s,map_p);
map_s = addmaps(map_s,map_y);

% Get noise maps and add to signal map (add systematics map here if necessary)
map_n = get_noise_map(obs_expt,expt,rlzn,planck_jacktype);
map = addmaps(map_s,map_n);

% Only return certain fields, necessary for array concatenation
map=stripmap(map);
map_n=stripmap(map_n);
map_lcdm=stripmap(map_lcdm);
map_r=stripmap(map_r);
map_d=stripmap(map_d);
map_s=stripmap(map_s);
map_p=stripmap(map_p);
map_y=stripmap(map_y);

if ~isempty(snum)
  outpath=fullfile('maps',num2str(snum,'%.4d'));
  if ~exist(outpath,'dir')
    mkdir(outpath);
  end
  fname=[num2str(rlz,'%.3d') '8_a_filtp3_weight3_gs_dp1100_jack' planck_jacktype '.mat'];
  save(fullfile(outpath,fname),'map','m','coaddopt');
end

return


%%%%%%%%%%%%%%%%%%%%%%%%
function map=addmaps(map1,map2)

% Expand as needed
n=max([size(map1,2) size(map2,2)]);
if n>1
  if size(map1,2)==1
    map1=repmat(map1,1,n);
  end
  if size(map2,2)==1
    map2=repmat(map2,1,n);
  end
end
map=map1;
for i=1:size(map,1)
  for j=1:size(map,2)
    map(i,j).T=map1(i,j).T+map2(i,j).T;
    map(i,j).Q=map1(i,j).Q+map2(i,j).Q;
    map(i,j).U=map1(i,j).U+map2(i,j).U;
  end
end

return


%%%%%%%%%%%%%%%%%%%%%%%%%%
function map=get_noise_map(obs_expt,expt,rlz,planck_jacktype)

switch expt(1)
 
 case 'b'
  if ~strcmp(planck_jacktype,'0')
    error(['BICEP2 or Keck maps should be made with planck_jacktype=0.  Asked for planck_jacktype=' planck_jacktype]);
  end
  switch(obs_expt)
    case 'bicep2',
      load(sprintf('/n/bicepfs2/bicep2/pipeline/maps/0751/%03d6_a_filtp3_weight3_gs_dp1100_jack0.mat',rlz));
      map=cal_map_wrapper(make_map(ac,m),m,coaddopt);
    case 'keck',
      load(sprintf('/n/panlfs2/bicep/keck/pipeline/maps/1353/%03d6_ab_filtp3_weight3_gs_dp1100_jack0.mat',rlz));
      map=cal_map_wrapper(make_map(ac,m),m,coaddopt);
    otherwise, error(['Unknown experiment ' obs_expt']); 
  end
 case 'p'
  switch(obs_expt)
    case 'bicep2',
      load(sprintf('/n/bicepfs2/b2planck/pipeline/maps/1613/%03d6_a_filtp3_weight3_gs_dp1100_jack%s.mat',rlz,planck_jacktype));
      fac=[30,44,70,100,143,217,353];
      ind=find(fac==str2num(expt(2:end)));
      map=cal_map_wrapper(make_map(ac(ind,:),m),m,coaddopt);
    case 'keck',
      load(sprintf('/n/bicepfs2/b2planck/pipeline/maps/1614/%03d6_ab_filtp3_weight3_gs_dp1100_jack%s.mat',rlz,planck_jacktype));
      fac=[30,44,70,100,143,217,353];
      ind=find(fac==str2num(expt(2:end)));
      map=cal_map_wrapper(make_map(ac(ind,:),m),m,coaddopt);
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
function map=get_lcdm_map(expt,rlz)

switch(expt)
  case 'bicep2',
load(sprintf('/n/bicepfs2/bicep2/pipeline/maps/0751/%03d5_a_filtp3_weight3_gs_dp1100_jack0.mat',rlz));
map=cal_map_wrapper(make_map(ac,m),m,coaddopt);
  case 'keck',
load(sprintf('/n/panlfs2/bicep/keck/pipeline/maps/1353/%03d5_ab_filtp3_weight3_gs_dp1100_jack0.mat',rlz));
map=cal_map_wrapper(make_map(ac,m),m,coaddopt);
end

return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function map=get_bmode_map(expt,rlz)

switch(expt)
  case 'bicep2',
load(sprintf('/n/bicepfs2/bicep2/pipeline/maps/0751/%03d4_a_filtp3_weight3_gs_dp1100_jack0.mat',rlz));
map=cal_map_wrapper(make_map(ac,m),m,coaddopt);
  case 'keck',
load(sprintf('/n/panlfs2/bicep/keck/pipeline/maps/1353/%03d4_ab_filtp3_weight3_gs_dp1100_jack0.mat',rlz));
map=cal_map_wrapper(make_map(ac,m),m,coaddopt);
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [map]=get_fgpowerlaw_map(expt,rlz,seedstr,alpha_BB,BB_to_EE)

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
      map=make_map(ac,m,coaddopt);
      map=cal_map_wrapper(map,m,coaddopt);
    end
    return
  end
end

error(['Suitable power law foreground maps not found!']);

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [map,coaddopt,m]=get_dust_map(expt)

switch(expt)
  case 'bicep2',
load('/n/bicepfs2/b2planck/pipeline/maps/psmdust/0012_a_filtp3_weight3_gs_dp1100_jack0.mat');
map=cal_map_wrapper(make_map(ac,m),m,coaddopt);
  case 'keck',
warning('No Keck reobserved PSM map available!');
load('/n/bicepfs2/b2planck/pipeline/maps/psmdust/0012_a_filtp3_weight3_gs_dp1100_jack0.mat');
map=cal_map_wrapper(make_map(ac,m),m,coaddopt);
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
function map_y=dust_sync_mix(map_p,map_y,epsilon)

fn=fieldnames(map_y);
for i=1:length(fn)
  if ismember(fn{i},{'x_tic','y_tic'})
    continue
  end
  map_y.(fn{i})=map_y.(fn{i})*(1-epsilon)+map_p.(fn{i})*(epsilon); % correct scaling?
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%5
function map=cal_map_wrapper(map,m,coaddopt)

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
end
map=cal_coadd_maps(map,calfac);

return

