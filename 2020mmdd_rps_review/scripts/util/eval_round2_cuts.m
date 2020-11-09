function [c,cp]=eval_round2_cuts(cp,c1,cut,p,ind)
% [c,cp]=eval_round2_cuts(cp,c1,cut,p,ind)
%
% Evaluate cuts which are applied at reduc_coaddpairmaps stage.
% The decision being made is what scanset/pairs to include in the
% coadded map and many of these are simple cuts.
%
% The overall philosophy is to delay calculating things until the last point
% at which it is possible to do so - to maximize the ability to make
% changes without re-running.
%
% Using the cut parameters in cp and the values and options in
% structure "cut" generate a series of cut masks in c - while the cut-params can
% be any size the output cut size must be something that expand_cuts
% knows how to make into nhalf-scan x nchannels.
% Creation of each cut mask is triggered by the existence of a field
% in structure "cut" which specifies a range/threshold etc. Each cut
% can be evaluated using one or more cut parameters according to
% arbitrary rules as coded below.
% 
% Once the cuts have been evaluated we expand them to full size and
% then combine to an overall mask. From this additional cut-on-cut
% parameters are generated and cut upon.
%
% There is a special complication with the passfrac_halfscan and
% scanset cuts - they depend on the overall round1 cut mask which is
% therefore passed in as an additional argument c1
%
% p is used to tell which channels are in which rx

disp('evaluating round2 cuts...');

% copy "utility" parameters
c.nhs=1; % cause expand cuts to do the right thing
c.nch=cp.nch;
c.nrx=cp.nrx;
% nmce added in BICEP3-era, so ensure compatibility with older data products
if isfield(cp,'nmce')
  c.nmce=cp.nmce;
else
  c.nmce=cp.nrx; % True for BICEP2 and Keck
end
%c.ntiles=cp.ntiles;

% manual cut from aux_data file
utc=cp.fs_t(1)+datenum(2005,5,17)-53507; % time of first halfscan
mc=ParameterRead('aux_data/manual_cuts.csv');
st=datenum(mc.starttime,'yyyymmdd HH:MM:SS');
et=datenum(mc.endtime,'yyyymmdd HH:MM:SS');
c.manual=true(1,c.nrx);
for i=1:length(st)
  c.manual(st(i)<utc&et(i)>utc,mc.rx(i)+1)=false;
end

% note that elnod_mean, fracdel and ab_ba could all be derived cut
% parameters added at this stage but they are not for historical
% reasons only

% mean of before/after elnods must be in specified range
if(isfield(cut,'elnod_mean'))
  c.elnod_mean=cp.elnod_mean>cut.elnod_mean(:,1)&...
               cp.elnod_mean<cut.elnod_mean(:,2);
  % exempt the darks (which will normally fail)
  c.elnod_mean(:,ind.d)=true;
  % both A and B of each pair must pass
  x=c.elnod_mean(ind.a)&c.elnod_mean(ind.b);
  c.elnod_mean(ind.a)=x; c.elnod_mean(ind.b)=x;
end

% fractional delta of before/after elnods must be below threshold
if(isfield(cut,'elnod_fracdel'))
  c.elnod_fracdel=cp.elnod_fracdel<cut.elnod_fracdel;
  % exempt the darks (which will normally fail)
  c.elnod_fracdel(:,ind.d)=true;
  % both A and B of each pair must pass
  x=c.elnod_fracdel(ind.a)&c.elnod_fracdel(ind.b);
  c.elnod_fracdel(ind.a)=x; c.elnod_fracdel(ind.b)=x;  
end

% ratio-of-the-ratio of a/b before/after values must be within
% threshold of unity
if(isfield(cut,'elnod_ab_ba'))
  c.elnod_ab_ba=cp.elnod_ab_ba<cut.elnod_ab_ba;
  % exempt the darks (which will normally fail)
  c.elnod_ab_ba(:,ind.d)=true;
  % both A and B of each pair must pass
  x=c.elnod_ab_ba(ind.a)&c.elnod_ab_ba(ind.b);
  c.elnod_ab_ba(ind.a)=x; c.elnod_ab_ba(ind.b)=x;  
end

% want no NaNs within el nods
if(isfield(cut,'elnod_nancount'))
  c.elnod_nancount=cp.elnod_nancount<cut.elnod_nancount;
  % both A and B of each pair must pass
  x=c.elnod_nancount(ind.a)&c.elnod_nancount(ind.b);
  c.elnod_nancount(ind.a)=x; c.elnod_nancount(ind.b)=x;
end

% cut on elnod goodness of fit 
if(isfield(cut,'elnod_gof') && isfield(cp,'elnod_gof')) 
  c.elnod_gof=cp.elnod_gof<cut.elnod_gof;
  % exempt the darks (which will normally fail)
  c.elnod_gof(:,ind.d)=true;
  % both A and B of each pair must pass
  x=c.elnod_gof(ind.a)&c.elnod_gof(ind.b);
  c.elnod_gof(ind.a)=x;c.elnod_gof(ind.b)=x;
end

% cut on elnod pair diff chisq, we don't want pair mismatched elnods
if(isfield(cut,'elnod_chisq_dif') && isfield(cp,'elnod_chisq'))
  c.elnod_chisq_dif=true(1, c.nch);
  val=cp.elnod_chisq(ind.b)<cut.elnod_chisq_dif;
  %else if the alternative minimum noise is requested
  if(isfield(cut,'elnod_altminnoise') && isfield(cp,'elnod_var'))
    val_alt=cp.elnod_var(ind.b)./cut.elnod_altminnoise.^2<cut.elnod_chisq_dif;
    val=val|val_alt;
  end 
  c.elnod_chisq_dif(ind.a)=val;c.elnod_chisq_dif(ind.b)=val;
  %exempt the darks
  c.elnod_chisq_dif(:,ind.d)=true;
end


% cut on detectors fractional resistance
% set bad values (from bad fits) to NaN and write cut in such a way that NaN passes 
if(isfield(cut,'rtes_frac'))
  cp.rtes_frac(cp.rtes_frac>1.0)=NaN;
  cp.rtes_frac(cp.rtes_frac<0.0)=NaN;
  c.rtes_frac=~(cp.rtes_frac<cut.rtes_frac(1) | cp.rtes_frac>cut.rtes_frac(2));
end

% cut on detectors normal resistance
% set bad values (from bad fits) to NaN and write cut in such a way that NaN passes 
if(isfield(cut,'rnorm'))
  cp.rnorm(cp.rnorm<0.0)=NaN;
  c.rnorm=~(cp.rnorm<cut.rnorm(1) | cp.rnorm>cut.rnorm(2));
end
% cut on detectors joule power
% set bad values (from bad fits) to NaN and write cut in such a way that NaN passes 
if(isfield(cut,'pjoule'))
  cp.pjoule(cp.pjoule<0.0)=NaN;
  c.pjoule=~(cp.pjoule<cut.pjoule(1) | cp.pjoule>cut.pjoule(2));
end


if isfield(cut,'elnod_median')
  freqs = unique(p.band(ind.la)');
  rxs = unique(p.rx');
  % Map a given receiver to a particular frequency.
  rxfreqs = arrayfun(@(rx) unique(p.band(p.rx'==rx & p.band'~=0 ...
      & ~isnan(p.band'))), rxs);

  % If elnod_median doesn't already exist, the "traditional"/historic
  % on-the-fly calculation from elnod_mean must be done first.
  if ~isfield(cp,'elnod_median')
    % initialize arrays
    cp.elnod_median=zeros(1,cp.nrx);

    for freq=freqs
      indrglf = sprintf('rgl%03d', freq);
      en_mean = cp.elnod_mean(ind.(indrglf));

      % Keck 2014 and later is recognizable by multiple frequencies, so
      % filter too small elnods.
      if length(freqs)>1
        en_mean(en_mean<100) = NaN;
      end

      cp.elnod_median(rxfreqs==freq) = nanmedian(en_mean);
    end
  end

  for rx=rxs
    freqidx = find(freqs == rxfreqs(rx+1));
    c.elnod_median(rx+1) = ...
        cp.elnod_median(rx+1)>cut.elnod_median(freqidx,1) ...
      & cp.elnod_median(rx+1)<cut.elnod_median(freqidx,2);
  end
end

% The 10 and 90th percentiles of the elnod values are constant
% fraction of the median unless weather is bad and channels start
% to saturate.
cp.elnod_90=prctile(cp.elnod_mean(ind.rgl),90)./cp.elnod_median;
if(isfield(cut,'elnod_90'))
  c.elnod_90=cp.elnod_90>cut.elnod_90(1)&...
             cp.elnod_90<cut.elnod_90(2);  
end
cp.elnod_10=prctile(cp.elnod_mean(ind.rgl),10)./cp.elnod_median;
if(isfield(cut,'elnod_10'))
  c.elnod_10=cp.elnod_10>cut.elnod_10(1)&...
             cp.elnod_10<cut.elnod_10(2);  
end

% white-noise and 1/f metrics
if(isfield(cut,'fb_wn_sd_p0'))
  % this cut applies only to paired channels - true by default
  c.fb_wn_sd_p0=true(1,c.nch);
  % take cut from diff
  val=cp.fb_wn_sd_p0(ind.b)<cut.fb_wn_sd_p0;
  % and apply to sum and diff
  c.fb_wn_sd_p0(ind.a)=val; c.fb_wn_sd_p0(ind.b)=val;
end
if(isfield(cut,'fb_1f_sd_p0'))
  % this cut applies only to paired channels - true by default
  c.fb_1f_sd_p0=true(1,c.nch);
  % take cut from diff
  val=cp.fb_1f_sd_p0(ind.b)<cut.fb_1f_sd_p0;
  % and apply to sum and diff
  c.fb_1f_sd_p0(ind.a)=val; c.fb_1f_sd_p0(ind.b)=val;
end

% Skewness:
if(isfield(cut,'skewness_dif'))
  % this cut applies only to paired channels - true by default
  c.skewness_dif=true(1,c.nch);
  % take cut from diff
  val=abs(cp.skewness(ind.b))<cut.skewness_dif;
  % and apply to sum and diff
  c.skewness_dif(ind.a)=val; c.skewness_dif(ind.b)=val;
end

% Satcom:
if(isfield(cut,'satcom'))
  % this cut applies only to paired channels - true by default
  c.satcom=true(1,c.nch);
  % take cut from diff
  val=cp.satcom(ind.b)<cut.satcom;
  % and apply to sum and diff
  c.satcom(ind.a)=val; c.satcom(ind.b)=val;
end

% std over full scanset
if(isfield(cut,'scanset_std'))
  % this cut applies only to paired channels - true by default
  c.scanset_std=true(1,c.nch);
  % Expand cut over bands if it is needed
  if numel(cut.scanset_std) == 1
    cut.scanset_std=[0 cut.scanset_std];
    cut.scanset_std=repmat(cut.scanset_std,length(freqs),1);
  end
  % take cut from diff
  val=[];
  for ii=1:cp.nrx
    % Identify the band of Rx, in ascending order of frequency
    band_order=find(unique(rxfreqs) == rxfreqs(ii));
    % Use "ismember" instead of "intersect" to keep order
    rx_indb=ind.b(ismember(ind.b,find(p.rx == ii-1)));
    val_tmp=cp.scanset_std(rx_indb) >= cut.scanset_std(band_order,1) & ...
            cp.scanset_std(rx_indb) <  cut.scanset_std(band_order,2);
    val=[val val_tmp];
  end
  % and apply to sum and diff
  c.scanset_std(ind.a)=val; c.scanset_std(ind.b)=val;
end

% Correlation across focal plane post p3 filter
if(isfield(cut,'fp_cor'))
  % don't cut if the entire focal plane's correlation was calculated as NaN
  c.fp_cor=cp.fp_cor<cut.fp_cor | isnan(cp.fp_cor);
end

% take the std of the std of the halfscans as a measure of noise
% stationarity - this is a derived cut parameter which can be added
% at this stage.
% For Gaussian noise this would be a constant number determined by the
% effective number of points in each half-scan. Pair diff is close to
% this while pair sum number is bigger and more variable due to 1/f
% atmospheric noise.
cp.stationarity_ab=nanstd(cp.fb_std_p3,1)./nanmean(cp.fb_std_p3,1);
cp.stationarity=nanstd(cp.fb_std_sd_p3,1)./nanmean(cp.fb_std_sd_p3,1);

if(isfield(cut,'stationarity_ab'))
  c.stationarity_ab=cp.stationarity_ab>cut.stationarity_ab(1)&...
                    cp.stationarity_ab<cut.stationarity_ab(2);
  % both A and B of each pair must pass
  x=c.stationarity_ab(ind.a)&c.stationarity_ab(ind.b);
  c.stationarity_ab(ind.a)=x; c.stationarity_ab(ind.b)=x;  

end

if(isfield(cut,'stationarity_sum'))
  % this cut applies only to paired channels - true by default
  c.stationarity_sum=true(1,c.nch);
  % take cut from sum
  val=cp.stationarity(ind.a)>cut.stationarity_sum(1)&...
      cp.stationarity(ind.a)<cut.stationarity_sum(2);
  % and apply to sum and diff
  c.stationarity_sum(ind.a)=val; c.stationarity_sum(ind.b)=val;  
end

if(isfield(cut,'stationarity_dif'))
  % this cut applies only to paired channels - true by default
  c.stationarity_dif=true(1,c.nch);
  % take cut from diff
  val=cp.stationarity(ind.b)>cut.stationarity_dif(1)&...
      cp.stationarity(ind.b)<cut.stationarity_dif(2);
  % and apply to sum and diff
  c.stationarity_dif(ind.a)=val; c.stationarity_dif(ind.b)=val;  
end

% focal plane temperature(s) must be in range
if(isfield(cut,'tfpu_mean'))
  c.tfpu_mean=cp.tfpu_mean>cut.tfpu_mean(1)&...
              cp.tfpu_mean<cut.tfpu_mean(2);
end

% focal plane temperature fluctuations must be less than threshold
if(isfield(cut,'tfpu_std'))
  c.tfpu_std=cp.tfpu_std<cut.tfpu_std;
end

% az encoder differences must be less than threshold
% (not present for B2 so don't try to calc even if asked for)
if(isfield(cut,'enc_az_diff')&&isfield(cp,'enc_az_diff'))
  c.enc_az_diff=abs(cp.enc_az_diff(:,5))<cut.enc_az_diff;
end

%az range cut
if(isfield(cut,'az_range') && isfield(cp,'az_range'))
  c.az_range=abs(cp.az_range)<cut.az_range;
end


% Evaluate cuts which are the analog of the per halfscan cuts but now
% working on per scanset basis - i.e. take the max cut parameter value
% over halfscans or similar - this is a bad way to do such cuts -
% better to go back to reduc_makepairmaps and do it explicitly per
% halfscan but maybe we want to anyway to avoid recompute

% total number of nans over all half-scan must be less than threshold
if(isfield(cut,'fb_nancount'))
  c.fb_nancount=sum(cp.fb_nancount)<cut.fb_nancount;
  % both A and B of each pair must pass
  x=c.fb_nancount(ind.a)&c.fb_nancount(ind.b);
  c.fb_nancount(ind.a)=x; c.fb_nancount(ind.b)=x;
end

% Cuts related to destepping
cp.num_fj=reshape(cp.num_fj,1,[]);
cp.num_destep=reshape(cp.num_destep,1,[]);

% Too many flux jumps in this channel
if(isfield(cut,'num_fj'))
  c.num_fj=cp.num_fj<=cut.num_fj;
  % both A and B of each pair must pass
  x=c.num_fj(ind.a)&c.num_fj(ind.b);
  c.num_fj(ind.a)=x; c.num_fj(ind.b)=x;
end

% Too many desteppings in this channel
if(isfield(cut,'num_destep'))
  c.num_destep=cp.num_destep<=cut.num_destep;
  % both A and B of each pair must pass
  x=c.num_destep(ind.a)&c.num_destep(ind.b);
  c.num_destep(ind.a)=x; c.num_destep(ind.b)=x;
end

% Neighbor channel behaving badly (but not consistently raging)
if(isfield(cut,'max_fj_gap'))
  % identify inconsistently raging channels
  x=(cp.num_fj>=10)&(cp.max_fj_gap<=cut.max_fj_gap);
  % all channels good by default
  y=true(size(cp.max_fj_gap));
  c.max_fj_gap=true(size(cp.max_fj_gap));
  for rx=unique(p.rx)'
    thisrx = p.rx==rx;
    for col=unique(p.mce_col(~isnan(p.mce_col)))'
      rp=(p.mce_col==col) & thisrx;
      y(rp)=y(rp)&~any(x(rp));
    end
    for row=unique(p.mce_row(~isnan(p.mce_row)))'
      % Pull out of loop to minimize repeated calculations.
      thismce_row = p.mce_row==row;
      coldiv8 = floor(p.mce_col/8);
      % Is this correct for BICEP3?  Do we still have 8 cols per RC?
      for rc=unique(coldiv8)'
        rp=thismce_row & thisrx & (coldiv8==rc);
        y(rp)=y(rp)&~any(x(rp));
      end
    end
  end
  % In principle, we should deal with "oddball" crosstalk examples here.
  % Need a mask of culprit/victim channels.

  % both A and B of each pair must pass
  x=y(ind.a)&y(ind.b);
  c.max_fj_gap(ind.a)=x; c.max_fj_gap(ind.b)=x;
end


% COPY ABOVE TYPE CUTS HERE FROM eval_round1_cuts as needed - and be
% careful to make sure A|B bad causes sum&diff bad


% do the expand here because we need to do a combine below
c=expand_cuts(c,p);


% DO PASSFRAC_HALFSCAN CUT

% find fraction of halfscans passing round 1 overall cut
cp.passfrac_halfscan=sum(c1.overall,1)/size(c1.overall,1);

% apply halfscan fraction cut
if(isfield(cut,'passfrac_halfscan'))
  c.passfrac_halfscan=cp.passfrac_halfscan>cut.passfrac_halfscan;
  % both A and B of each pair must pass
  x=c.passfrac_halfscan(ind.a)&c.passfrac_halfscan(ind.b);
  c.passfrac_halfscan(ind.a)=x; c.passfrac_halfscan(ind.b)=x;
end


% DO PASSFRAC_SCANSET CUT

% calculate overall second round mask
cm=combine_cuts(c);
% combine with first round to form nhalfscan*nchannel mask
cme=c1.overall&repmat(cm,[c1.nhs,1]);

% determine the fraction passing within each rx
for r=unique(p.rx)'
  % find the mask values for the good light channels in this rx
  x=cme(:,intersect(find((p.rx)==r),ind.gl));    
  % find the fraction passing
  cp.passfrac_scanset(r+1)=sum(x(:))/numel(x);
end

% apply scanset fraction cut
if(isfield(cut,'passfrac_scanset'))
  c.passfrac_scanset=true(size(cm));
  for r=unique(p.rx)'
    % test if above threshold for this rx
    y=cp.passfrac_scanset(r+1)>=cut.passfrac_scanset;
    % find all channels of this rx
    rp=(p.rx)==r;
    % write the logical mask value
    c.passfrac_scanset(:,rp)=y;
  end
end

return
