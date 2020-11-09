 function dPdT = calc_power_resp(band)
%
%    Input:
%        - band structure containing either:
%                - band.name: One of the following: 
%    'B2_150','Bicep2','B2','K','B2K','BK','BK14_150','BK15_150' for  150 band
%     'K95', 'BK14_95', 'BK15_95'  for 100 band 
%      'K220','BK15_220' for 220 band
%
%                - band.v_cen and band.frac_bw to use a tophat bandpass
%   Output:
%         - dPdT_rj  (in the Raleigh Jeans regime) 
%         - dPdT_cmb (in the CMB regime temperature regime)
%         - and their ratio to convert K_rj temps into K_cmb temps.

[h,k,c]= constants() ;
Tcmb = 2.72548;
% define the frequency axis from 0 to 500GHz
dv = 0.1e9 ;
v = 0.1e9:dv:500e9;
 
if(~exist('band'))
  band=[];
end
%% define bandpass
if(isfield(band,'name'))
  l.opt.expt = {band.name};
  l  = like_read_bandpass(l) ;
  if isempty(l.bandpass{1})
     disp('Band.name option must be one of the pre-defined expt name. See this doc to check')
     return
  else
    % interp onto common frequency axis
    bandpass = interp1(l.bandpass{1}(:,1)*1e9,l.bandpass{1}(:,2),v);
  end
  
elseif (isfield(band,'v_cen')) & (isfield(band,'frac_bw'))
  disp('Using a top hand bandpass')
  bandpass = ones([1,numel(v)]);
  v_lo = band.v_cen-band.v_cen*band.frac_bw/2;
  v_hi = band.v_cen+band.v_cen*band.frac_bw/2;
  bandpass(find(v<v_lo)) = 0;
  bandpass(find(v>v_hi)) = 0;

else
  disp('band.name or band.v_cen/frac_bw must be define to obtain a responsivity')
  return
end

% March 27th: dB replaces the numerical derivative 
%with the planck_dIdT.m function which provides the  analytical derivative
% Below is the code for the analyitical result
% the factor of 0.5* c^2/v^2 is to convert to single moded power per unit freq.
dPdTv = @(v,T) bandpass .* planck_dIdT(v,T) * 0.5 * c^2 ./v.^2 ; 

dPdT.cmb = nansum(dPdTv(v,Tcmb))*dv;
dPdT.rj = nansum(dPdTv(v,1000))*dv;
dPdT.rj2cmb = dPdT.cmb/dPdT.rj;
dPdT

% Below is the code for the numerical derivative which gives the same results
%%  derive dP/dTcmb power responsivity to a small CMB temperature variation
%dT = 0.04;
%Qv_cmb =         planck_v(v,Tcmb,1) .* bandpass;
%Qv_cmb_plus_dT = planck_v(v,Tcmb+dT,1) .* bandpass ;
%dPdT.cmb = (nansum(Qv_cmb_plus_dT).*dv - nansum(Qv_cmb).*dv)./dT ; % W/K_cmb
%
%%  derive dP/dTrj power responsivity in K_rj regime.
%dT = 0.1;
%Qv_cmb_plus_dT = planck_v(v,1000+dT,1) .* bandpass ;
%Qv_cmb =         planck_v(v,1000,1).* bandpass;
%dPdT.rj = (nansum(Qv_cmb_plus_dT)* dv - nansum(Qv_cmb)* dv)./dT ; % W/K_cmb
%dPdT.rj2cmb = dPdT.cmb/dPdT.rj;
%dPdT
