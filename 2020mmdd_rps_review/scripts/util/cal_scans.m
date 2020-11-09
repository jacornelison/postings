function [v,en]=cal_scans(v,fs,en,indg,indrg,yr,lc,bias0,p,ind,tesbias)
% v=cal_scans(v,fs,en,indg,indrg,yr,lc,bias0,p,ind,tesbias)
%
% Relative cal the channels
% fs   =  scan indices
% en   =  structure of the elnod calibration data
% indg =  the set of channels to apply cal to (normally ind.gl)
% indrg = the set of channels to take median gain over (normally ind.rgl)
% yr = the year, for removing small values in data from 2014 onward
% lc   =  structure of the load curve data, used for 2014 and later data
% bias0 = a frequency dependent rescaling when calculating median in power units
%         indexed in frequency-ascending order and
%         numel(bias0) = numel(unique(p.band(ind.la)))
% p and ind = structures returned from get_array_info
% tesbias = the biases of the detectors during the scanset for use when converting to power units
%
% This function will delibaretly fail if we're in 2014 or later but the new input
% variables (lc,bias0,p,ind,tesbias) aren't passed in
%
% No test for before/after gain stability is done here - this is
% deferred to the mapping stage
%
% By normalizing by the median this function removes common mode elnod
% response changes which are due to changes in atmospheric opacity (or
% overall instrument gain).
% Notice that in BICEP2 data trends among groups of detectors remain
% so this may not be OK.

disp('cal_scans...')

% get the elnod gains
g=en.g(:,:,2);
% replicate to a new variable so that data are not arbitrarily NaN'd during relgain
% calibration
gv=g;

% -ve values make no sense
g(g<0)=NaN;

% identical zero values are channels which were masked in MCE software
g(g==0)=NaN;

% GPT - also get rid of small values to exclude whole rx that will get cut
if ~exist('yr','var') || isempty(yr)
  yr=2012;
end
if ischar(yr)
  yr=str2num(yr);
end
expt=get_experiment_name();
if strcmp(expt,'keck') && yr>=2014
  g(g<100)=NaN;
end

% Convert to power units if 2014 or beyond
if yr>=2014 && exist('lc','var') && isfield(lc,'g') && ~isempty(lc.s)
  pu=true;
  % need another copy of the gains with the nans for median calculation
  gains=g;
else
  pu=false;
end

if pu
  % We want to convert fb units to nominal power units with the dP/dI correction,
  % according to Walt's Post "elnod response in power units" in Keck logbook.
  % http://bicep0.caltech.edu/~spuder/keck_analysis_logbook/analysis/20141106_elnod_power/
  % See also "RelGain Correction with Maps Pager"
  % http://bicep0.caltech.edu/~spuder/keck_analysis_logbook/analysis/20141125_RelGain_Correction/
  % Here are the values we need and where they come from:
  % r_shunt: Take mean value from calib file, or assume 0.003 Ohm if not available.
  % r_sensor: Comes from lc.g(:,:,4).  Note this value is not saved in the cut_params
  % file. From that file, we would have to use rtes_frac and rnorm to recover r_sensor.
  % tes_bias: want the bias as it was recorded at the start of the first field scan, so
  % need d.mce0.tes.bias(fs.s(1),:)
  % The factor we need to multiply the fb units by is
  % tes_bias*(1-r_shunt/r_sensor)/(1+r_shunt/r_sensor) 
  % elnod fits also need to be converted to power units. Currently the goodness of fit
  % parameter is still in the old units.  We'd probably have to re-run elnod to get the
  % changed units.  This shouldn't make much difference.  Cut thresholds can be
  % changed, if necessary.
  disp('converting to nominal units of power to compute median of elnod gains')
  bias=double(tesbias(fs.s(1),:));
  % get values from partial load curve analysis lc.g
  if isfield(lc,'r_sh')
    rsh=nanmean(lc.r_sh); %r_shunt
  else
    rsh=0.003; %r_shunt
  end
  rs=nanmean(lc.g(:,:,4),1); % r_sensor
  rn=nanmean(lc.g(:,:,3),1); % r_normal
  % p.mce_col runs between 0 and 15 (0 to 30 for B3).
  % Convert this to run between 1 and 80 (diff max. for B3)so that
  % channels can be matched with their bias value
  indtes=(max(p.mce_col+1)*p.mce+p.mce_col+1)';
    
  recalfactor=(1-rsh./rs)./(1+rsh./rs); % This is the dP/dI correction
  rcf=nan(size(recalfactor)); % The value the elnod gains will be multiplied by
  rfrac=rs./rn;  rfrac(rfrac<0.0 | rfrac>1.0)=NaN;

  % apply recalfactor and a bias factor to change the values of these units of power
  % to similar values they would have been if left in feedback units.  It should be
  % about 2600 for 150 GHz and 1000 for 100 GHz (to the nearest 100).  Only do this
  % for rgl, since those are the only values used to calculate the median. Also
  % because dark detectors have band=0 for 150 detectors in 2014.

  freq=unique(p.band(indg)); % should produce error if indices with multiple
                             % frequencies are incorrectly passed to this function
   
  % which index of bias0 to use for this frequency
  bind=find(unique(p.band(ind.la))==freq);
  
  for k=indrg

    rcf(k)=bias(indtes(k))*recalfactor(k)/bias0(bind);
    % make sure we're not multiplying by a negative number. Maybe will want to come
    % back and make sure we're not multiplying by an absurdly large number, but cuts
    % might take care of that.
    % Also, don't allow unreliable values in the median calculation
    if rs(k)<= rsh | isnan(rfrac(k)) ; rcf(k)=NaN; end
    gains(:,k)=gains(:,k)*rcf(k);
  end
  
  %set up field names 
  rglf=sprintf('rgl%03d', freq);
  gmedf=sprintf('med%03d', freq);
  en.gmed.(gmedf)=nanmedian(gains(:,ind.(rglf)),2);
  
end

% Take ratio of each channel's gain to median of good channels at each
% time step - we use the median to give some immunity to crazy
% channels and/or channels which come and go.

if ~pu
  rg=g./repmat(nanmedian(g(:,indrg),2),[1,size(g,2)]);
else
  % gmedf should be correct fieldname for the frequency.  We should only
  % be looking at one frequency at a time in this function  
  rg=gv./repmat(en.gmed.(gmedf),[1,size(g,2)]);
end

for j=1:length(fs.t)
  
  % find the preceding elnod
  x=fs.t(j)-en.t; x(x<0)=NaN;
  [dummy,b]=nanmin(x); if(isnan(dummy)) b=[]; end
  % find the subsequent elnod
  x=fs.t(j)-en.t; x(x>=0)=NaN;
  [dummy,a]=nanmax(x); if(isnan(dummy)) a=[]; end
  [b a]
  
  % For each detector scale to match median
  for i=indg
    s=fs.sf(j);  e=fs.ef(j);
    
    % take the mean before/after gain
    rg1=rg(b,i); rg2=rg(a,i);
    m=mean([rg1,rg2]);
    
    % apply the gain cal regardless of stability
    v(s:e,i)=v(s:e,i)/m;
    
  end
end

return
