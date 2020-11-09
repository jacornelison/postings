function fix_ukpervolt(map,ukpervolt,pre_ukpervolt)
% function fix_ukpervolt(map,ukpervolt,pre_ukpervolt)
% written to fix simulation ukpervolt
% 
% map = string specifying the maps to fix.  takes wildcards
% ukpervolt = new calibration
%  
% if previously a wrong rescaling attempt with cal_coadd_ac(...,xx,cc,1,1)
% that scaled variance and signal was done, pre_ukpervolt allows to roll this 
% back and then apply the fix
%
% NOTE NOTE NOTE
% This is all awfull! Please be very sure you know what you are doing. This function 
%  make an in place replacement of potentially very valuable files.
%  
% example: 
% ukpervolt=get_default_ukpervolt('2014');
% fix_ukpervolt('1351/001[245]_c_filtp3_weight3_gs_dp1100_jack?1.mat',ukpervolt);
%

error('Are you sure you want to use this function?')

if ~exist('map','var') || isempty(map)
  error('Must input map!')
end
if ~exist('ukpervolt','var') || isempty(ukpervolt)
  ukpervolt=get_ukpervolt;
  warning('Using default ukpervolt from 2012');
end
  
% find the requested maps
[s d]=system_safe(['/bin/ls maps/' map]);
d=strip_nonsense(d);
maps1=strread(d,'%s');
nmaps=length(maps1);

for ii=1:nmaps

  mapn=maps1{ii};
  % check type - this should only be applied to signal-only sims
  type=str2num(mapn(14));
  if ~(type>=2 & type<=5)
    warning(['written only to fix the ukpervolt of signal-only sim! Skipping ' mapn])
  end

  disp(['Loading ' mapn]);
  load(mapn);
  
  ukpv_old=coaddopt.mapopt{1}.simopt.ukpervolt
  
  % check if it needs updating
  % assumes that ukpv_old is either the same size as ukpervolt
  % or that it is a singular value.  
  % now that simopt.ukpervolt can take nelements, this is not 
  % necessarily true....
  if isfield(coaddopt.mapopt{1}.simopt,'simukpervolt') && strcmp(coaddopt.mapopt{1}.simopt.simukpervolt,'fixed_20150506')
    disp(['Skip since ',coaddopt.mapopt{1}.simopt.simukpervolt])
    clear ac;
    clear m;
    clear coaddopt;
    continue
  end
  
  % wrong cal_coadd_ac including variance rescaling was previously applied,
  % fix this first by going back into the old ukpervolt:
  if exist('pre_ukpervolt','var') && ~isempty(pre_ukpervolt)
    ac=cal_coadd_ac(ac,ukpv_old/pre_ukpervolt,coaddopt,1,1);
    ukpv_old = pre_ukpervolt;
  end
  
  % re calibrate, only apply to signal portion:
  ac=cal_coadd_ac(ac,ukpv_old./ukpervolt,coaddopt,0,1);

  % update the simopt for future reference
  coaddopt.mapopt{1}.simopt.ukpervolt=ukpervolt;
  
  coaddopt.mapopt{1}.simopt.simukpervolt = 'fixed_20150506';  

  % resave
  saveandtest(mapn,'ac','m','coaddopt','-v7.3');
  setpermissions(mapn);

  % update user
  disp(['Fixed ' mapn]);
  
  
  clear ac;
  clear m;
  clear coaddopt;

end



