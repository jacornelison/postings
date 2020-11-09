function cut=get_default_round2_cuts(expt,year)
% cut=get_default_round2_cuts(expt,year)
%
% get a set of vanilla cuts
% - all available should be listed here - although some are set to "inf" etc
%
% added year option in for keck as well
%

if(~exist('expt','var') || isempty(expt))
  expt=get_experiment_name;
end
if(~exist('year','var'))
  year='2012';
end

cut.elnod_fracdel=0.3;

switch expt
 case 'bicep2'
  cut.elnod_ab_ba=0.03;
 case 'keck'
  cut.elnod_ab_ba=0.04;
 case 'bicep3'
  cut.elnod_ab_ba=0.04;
end

cut.elnod_nancount=1;

% these have to differ B2/Keck
switch expt
 case 'bicep2'
  cut.elnod_mean=[2000,8000];
 case 'keck'
  cut.elnod_mean=[1000,30000];
 case 'bicep3'
  cut.elnod_mean=[90,10000];
end
switch expt
 case 'bicep2'
   cut.elnod_median=[3000,6000];
 case 'keck'
  switch year
   case '2012'
    cut.elnod_median=[2200,7000];
   case '2013'
    cut.elnod_median=[2200,7000];
   case '2014'
    % rows in frequency ascending order (100GHz then 150GHz)
    cut.elnod_median=[4500,10000; 3000,9000];
   case '2015'
    cut.elnod_median=[4500,10000; 3000,9000; 2000,8000];
   case '2016'
    % Preliminary
    % 150GHz 210GHz and 220GHz
    cut.elnod_median=[3000,9000; 3500,13000; 2000,8000];
   case '2017'
    % Preliminary for 210, 220, 270 GHz
    cut.elnod_median=[3500,13000; 2000,8000; 2000,15000];
   case '2018'
    % Preliminary for 210, 220, 270 GHz
    cut.elnod_median=[3500,13000; 2000,8000; 2000,15000];
   case '2019'
    % Preliminary for 150 (copied from 2016), 210, 220, 270 GHz
    cut.elnod_median=[3000,9000; 3500,13000; 2000,8000; 2000,15000];
  end
 case 'bicep3'
  cut.elnod_median=[2000,5000];
end

switch expt
 case 'bicep2'
  cut.elnod_gof=100;
 case 'keck'
  cut.elnod_gof=250;
 case 'bicep3'
  cut.elnod_gof=75;
end

switch expt
 case 'keck'
  switch year
   case {'2017','2018','2019'}
    % Preliminary
    cut.elnod_chisq_dif=10;
    cut.elnod_altminnoise=0.15;    
   case '2016'
    % Preliminary
    cut.elnod_chisq_dif=10;
    cut.elnod_altminnoise=0.15;
   case '2015'
    % Preliminary
    cut.elnod_chisq_dif=10;
    cut.elnod_altminnoise=0.15;
   case '2014'
    cut.elnod_chisq_dif=10;
    cut.elnod_altminnoise=0.15;
   case '2013'
    cut.elnod_chisq_dif=10;
    cut.elnod_altminnoise=0.15;
   case '2012'
    cut.elnod_chisq_dif=8;
    cut.elnod_altminnoise=0;
  end
 case 'bicep2'
  cut.elnod_chisq_dif=8;
  cut.elnod_altminnoise=0;
 case 'bicep3'
  cut.elnod_chisq_dif=10;
  cut.elnod_altminnoise=0;
end

cut.rtes_frac=[0.1,0.95];
cut.rnorm=[0,Inf];
cut.pjoule=[0,Inf];

% abandon these as extra complexity which adds little
%cut.elnod_90=[1.1,1.5];
%cut.elnod_10=[0.7,0.9];

cut.fp_cor=1;
cut.skewness_dif=0.2;

cut.scanset_std=2.5;
if strcmp(expt, 'keck') && ismember(year, {'2018'})
  % Make exception for Keck 270 bands
  % In ascending order: 210, 220 and 270 GHz
  cut.scanset_std = [0,2.5;0,2.5;0,5.0];
elseif strcmp(expt, 'keck') && ismember(year, {'2019'})
  % In ascending order: 150, 210, 220, and 270 GHz
  cut.scanset_std = [0,2.5;0,2.5;0,2.5;0,5.0];
end

cut.fb_wn_sd_p0=Inf;
cut.fb_1f_sd_p0=Inf;

% cuts related to destepping
cut.num_fj=5;
cut.num_destep=5;
cut.max_fj_gap=1000;

cut.stationarity_ab=[0,0.7];
cut.stationarity_dif=[0,0.2];

switch expt
 case {'bicep2','keck'}
  switch year
   case {'2017','2018','2019'}
    cut.tfpu_mean=[0.2,0.35];
   otherwise
    cut.tfpu_mean=[0.2,0.3];
  end
  case 'bicep3'
    cut.tfpu_mean=[0.25,0.35];
end
cut.tfpu_std=5d-5;
cut.enc_az_diff=3e4;
cut.az_range=100;

cut.passfrac_halfscan=0.9;

switch expt
 case 'bicep2'
  cut.passfrac_scanset=0.5;
 otherwise
  cut.passfrac_scanset=0.3;
end

switch expt
 case 'keck'
  switch year
   case {'2017','2018','2019'}
    % Preliminary
    cut.satcom=6;
   case '2016'
    % Preliminary
    cut.satcom=6;
   case '2015'
    cut.satcom=6;
   case '2014'
    cut.satcom=6; % GPT 20140702 - so far define only for Keck 2014
  end
 case 'bicep3'
  cut.satcom=6;
end

return
