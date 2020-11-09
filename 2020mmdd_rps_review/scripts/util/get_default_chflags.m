function chflags=get_default_chflags(expt,year)
% chflags=get_default_chflags(expt,year)
%
% get the default channel flags for Keck and B2 and BICEP3

if(~exist('expt','var') || isempty(expt))
  expt=get_experiment_name;
end
if(~exist('year','var'))
  year=2012;
end
if ~isnumeric(year)
  year=str2num(year);
end

switch expt
 case 'bicep2'
  chflags.filebase={'fp_data/fp_data_ukpvflags', 'fp_data/fp_data_ukpvflags', 'fp_data/fp_data_aboffsetflags', 'fp_data/fp_data_beamshapeflags',...
   'fp_data/fp_data_beamshapeflags', 'fp_data/fp_data_beamshapeflags', 'fp_data/fp_data_beamshapeflags', ...
    'fp_data/fp_data_contaminatedpairs'};
  chflags.par_name={'ukpv', 'ukpvaoverb', 'aboffset_pct', 'sigma_pct', 'ellip', 'sigma_diff_pct', 'pc_diff', 'pairdiff_contam'};
  chflags.low=[3000 0.9550 0 0.9000 0 -0.0700 0 -1];
  chflags.high=[3750 1.0450 0.2000 1.1000 0.1500 0.0700 0.0900 0.27911];

 case 'keck'
   
   chflags.filebase={'fp_data/fp_data_radioptflags','fp_data/fp_data_ukpvflags', 'fp_data/fp_data_ukpvflags','fp_data/fp_data_aboffsetflags', 'fp_data/fp_data_beamshapeflags','fp_data/fp_data_beamshapeflags', 'fp_data/fp_data_beamshapeflags', 'fp_data/fp_data_beamshapeflags','fp_data/fp_data_contaminatedpairs'};
   
   switch year
   case 2012
    chflags.par_name={'pos_err','ukpv_pct', 'ukpvaoverb','aboffset_pct', 'sigma_pct', 'ellip', 'sigma_diff_pct', 'pc_diff', 'pairdiff_contam'};
    chflags.low=[0 0.9 .93 0 0.9 0 -0.07 0 -1];
    chflags.high=[0.1 1.25 1.07 0.3 1.1 0.15 0.07 0.09 0.14043];
   case 2013
    chflags.par_name={'pos_err','ukpv_pct', 'ukpvaoverb','aboffset_pct', 'sigma_pct', 'ellip', 'sigma_diff_pct', 'pc_diff', 'pairdiff_contam'};
    chflags.low=[0 0.9 .93 0 0.9 0 -0.07 0 -1];
    chflags.high=[0.1 1.25 1.07 0.3 1.1 0.15 0.07 0.09 0.15105];
   case 2014
    chflags.par_name={'pos_err','ukpv_pct', 'ukpvaoverb','aboffset_pct', 'sigma_pct', 'ellip', 'sigma_diff_pct', 'pc_diff', 'pairdiff_contam'};
    chflags.low=[0 0.9 .93 0 0.9 0 -0.07 0 -1];
    chflags.high=[0.1 1.1 1.07 0.3 1.1 0.15 0.07 0.09 0.14474];
   case 2015
    chflags.par_name={'pos_err','ukpv_pct', 'ukpvaoverb','aboffset_pct', 'sigma_pct', 'ellip', 'sigma_diff_pct', 'pc_diff', 'pairdiff_contam'};
    % get_array_info applies channel flags while building [p,ind], so limits
    % can be specified at per-rx granularity since each rx loads independent
    % fp_data file. We emulate per-frequency limits by setting each per-rx
    % limit appropriately.
    chflags.low=[
      0.00 0.90 0.93 0.0 0.9 0.00 -0.07 0.00 0.0000 %  95GHz
      0.00 0.65 0.90 0.0 0.9 0.00 -0.07 0.00 0.0000 % 220GHz
      0.00 0.90 0.93 0.0 0.9 0.00 -0.07 0.00 0.0000 %  95GHz
      0.00 0.65 0.90 0.0 0.9 0.00 -0.07 0.00 0.0000 % 220GHz
      0.00 0.90 0.93 0.0 0.9 0.00 -0.07 0.00 0.0000 % 150GHz
      ];
    chflags.high=[
      0.10 1.10 1.07 0.3 1.1 0.15  0.07 0.09 0.13358 %  95GHz
      0.15 1.35 1.10 0.3 1.1 0.15  0.07 0.09 0.13358 % 220GHz
      0.10 1.10 1.07 0.3 1.1 0.15  0.07 0.09 0.13358 %  95GHz
      0.15 1.35 1.10 0.3 1.1 0.15  0.07 0.09 0.13358 % 220GHz
      0.10 1.10 1.07 0.3 1.1 0.15  0.07 0.09 0.13358 % 150GHz
      ];
   case 2016
    chflags.par_name={'pos_err','ukpv_pct', 'ukpvaoverb','aboffset_pct', 'sigma_pct', 'ellip', 'sigma_diff_pct', 'pc_diff', 'pairdiff_contam'};
    % Adopt the same approach applied in K2015 - frequency-dependent thresholds
    chflags.low=[
      0.00 0.65 0.90 0.0 0.9 0.00 -0.07 0.00 0.0000 % 210GHz
      0.00 0.65 0.90 0.0 0.9 0.00 -0.07 0.00 0.0000 % 220GHz
      0.00 0.65 0.90 0.0 0.9 0.00 -0.07 0.00 0.0000 % 210GHz
      0.00 0.65 0.90 0.0 0.9 0.00 -0.07 0.00 0.0000 % 220GHz
      0.00 0.90 0.93 0.0 0.9 0.00 -0.07 0.00 0.0000 % 150GHz
      ];
    chflags.high=[
      0.15 1.35 1.10 0.3 1.1 0.15  0.07 0.09 0.13332 % 210GHz
      0.15 1.35 1.10 0.3 1.1 0.15  0.07 0.09 0.13332 % 220GHz
      0.15 1.35 1.10 0.3 1.1 0.15  0.07 0.09 0.13332 % 210GHz
      0.15 1.35 1.10 0.3 1.1 0.15  0.07 0.09 0.13332 % 220GHz
      0.10 1.10 1.07 0.3 1.1 0.15  0.07 0.09 0.13332 % 150GHz
      ];
   case 2017
    chflags.par_name={'pos_err','ukpv_pct', 'ukpvaoverb','aboffset_pct', 'sigma_pct', 'ellip', 'sigma_diff_pct', 'pc_diff', 'pairdiff_contam'};
    % Adopt the same approach applied in K2015 - frequency-dependent thresholds
    chflags.low=[
      0.00 0.65 0.90 0.0 0.9 0.00 -0.07 0.00 0.0000 % 210GHz
      0.00 0.65 0.90 0.0 0.9 0.00 -0.07 0.00 0.0000 % 220GHz
      0.00 0.65 0.90 0.0 0.9 0.00 -0.07 0.00 0.0000 % 210GHz
      0.00 0.65 0.90 0.0 0.9 0.00 -0.07 0.00 0.0000 % 220GHz
      0.00 0.65 0.90 0.0 0.9 0.00 -0.07 0.00 0.0000 % 270GHz - copy from 220 GHz 
      ];
    chflags.high=[
      0.15 1.35 1.10 0.3 1.1 0.15  0.07 0.09 0.1299 % 210GHz
      0.15 1.35 1.10 0.3 1.1 0.15  0.07 0.09 0.1299 % 220GHz
      0.15 1.35 1.10 0.3 1.1 0.15  0.07 0.09 0.1299 % 210GHz
      0.15 1.35 1.10 0.3 1.1 0.15  0.07 0.09 0.1299 % 220GHz
      0.15 1.35 1.10 0.3 1.1 0.15  0.07 0.09 0.1299 % 270GHz - copy from 220 GHz 
      ];
   case 2018
    chflags.par_name={'pos_err','ukpv_pct', 'ukpvaoverb','aboffset_pct', 'sigma_pct', 'ellip', 'sigma_diff_pct', 'pc_diff', 'pairdiff_contam'};
    % Adopt the same approach applied in K2015 - frequency-dependent thresholds
    % 20190121: temporarily loosen cuts of chflags from FFBM. 
    % It will be updated once we got products from K2018 FFBM analysis.
    chflags.low=[
      0.00 0.65 0.90 -10 -10 -10 -10 -10 0.0000 % 210GHz
      0.00 0.65 0.90 -10 -10 -10 -10 -10 0.0000 % 220GHz
      0.00 0.65 0.90 -10 -10 -10 -10 -10 0.0000 % 210GHz
      0.00 0.65 0.90 -10 -10 -10 -10 -10 0.0000 % 220GHz
      0.00 0.65 0.80 -10 -10 -10 -10 -10 0.0000 % 270GHz
      ];
    chflags.high=[
      0.15 1.35 1.10 10 10 10 10 10 0.12735 % 210GHz
      0.15 1.35 1.10 10 10 10 10 10 0.12735 % 220GHz
      0.15 1.35 1.10 10 10 10 10 10 0.12735 % 210GHz
      0.15 1.35 1.10 10 10 10 10 10 0.12735 % 220GHz
      0.60 1.35 1.20 10 10 10 10 10 0.12735 % 270GHz
      ];
   otherwise
    error('Channel flags requested for year that is not implemented')
   end
   
 case 'bicep3'
  chflags.filebase={'fp_data/fp_data_radioptflags','fp_data/fp_data_ukpvflags', 'fp_data/fp_data_ukpvflags','fp_data/fp_data_aboffsetflags', 'fp_data/fp_data_beamshapeflags','fp_data/fp_data_beamshapeflags', 'fp_data/fp_data_beamshapeflags', 'fp_data/fp_data_beamshapeflags','fp_data/fp_data_contaminatedpairs'};

  switch year
   case 2016
    chflags.par_name={'pos_err','ukpv_pct', 'ukpvaoverb','aboffset_pct', 'sigma_pct', 'ellip', 'sigma_diff_pct', 'pc_diff', 'pairdiff_contam'};
    % Adopt the same similar number in Keck 95 for now
    % changed ukpv_pct, ukpvaoverb, and pairdiff_contam
    chflags.low=[
      0.00 0.85 0.92 0.0 0.9 0.00 -0.07 0.00 0.0000  
      ];
    chflags.high=[
      0.10 1.15 1.08 0.3 1.1 0.15  0.07 0.09 0.1557
      ]; 
   case 2017
    chflags.par_name={'pos_err','ukpv_pct', 'ukpvaoverb','aboffset_pct', 'sigma_pct', 'ellip', 'sigma_diff_pct', 'pc_diff', 'pairdiff_contam'};
    % use the same as B2016
    chflags.low=[
      0.00 0.85 0.92 0.0 0.9 0.00 -0.07 0.00 0.0000  
      ];
    chflags.high=[
      0.10 1.15 1.08 0.3 1.1 0.15  0.07 0.09 0.1961
      ]; 
   case 2018
    % copy 2018
    chflags.par_name={'pos_err','ukpv_pct', 'ukpvaoverb','aboffset_pct', 'sigma_pct', 'ellip', 'sigma_diff_pct', 'pc_diff', 'pairdiff_contam'};
    chflags.low=[
      0.00 0.85 0.92 0.0 0.9 0.00 -0.07 0.00 0.0000  
      ];
    chflags.high=[
      0.10 1.15 1.08 0.3 1.1 0.15  0.07 0.09 0.1961
      ]; 
   otherwise
    error('Channel flags requested for year that is not implemented')
   end 

end
