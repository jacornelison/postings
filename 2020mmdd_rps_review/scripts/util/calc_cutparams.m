function cp=calc_cutparams(d,fs,en,lc,p,ind,dg)
% cp=calc_cutparams(d,fs,en,lc,p,ind,dg)
%
% Calculates cut parameter values to be used in
% accepting and rejecting individual half-scans
% for each detector.
% Called by reduc_plotcuts (for investigation) and reduc_makepairmaps
% (for real).
%
% Cut parameters are derived from timestream data and are typically
% floating point values. They can be of varying size depending on the
% quantity. They may correspond directly to a cut or be used in
% combination to derive cuts. They may also not end up being used for
% a cut at all, but just be a quantity which it is interesting to plot
% over time.
% 
% Most cut parameters are calculated from non-pair-diff data. However
% some are calculated from tod with has been sum/diff'ed. In this case
% the user of the quantity (normally eval_cuts) needs to know this.
%
% p is used to determine which channels are in which rx
% ind is used at the moment only to determine the sets of a/b channels -
% please try very hard to avoid adding anything which depends on the sets
% of "good" channels in ind as we don't want to commit to which these are
% until reduc_coaddpairmaps.
%
% dg is the deglitching structure saved by reduc_initial.

expt=get_experiment_name();
disp(['calculating cut parameters for ' expt '...'])

% keep a record of relevant sizes
cp.nhs=length(fs.s);
cp.nch=size(d.mce0.data.fb,2);
cp.nrx=length(unique(p.rx));
cp.nmce=length(unique(p.mce));

% useful receiver <-> frequency mappings
freqs = unique(p.band(ind.la)');
rxs = unique(p.rx');
% Map a given receiver to a particular frequency.
rxfreqs = arrayfun(@(rx) unique(p.band(p.rx'==rx & p.band'~=0 & ...
    ~isnan(p.band'))), rxs);

% precompute the p0/p3 before/after sum/diff -
% this uses loads of extra memory but we need in various places and
% slow to recompute
d_p0=filter_scans(d,fs,'p0');
d_p3=filter_scans(d,fs,'p3');

% do the sum diff over whole scanset - relevant to see if steps
% corrupting noise model, as well as over the elnods.

fsp.sf(1)=en.sf(1); fsp.ef(1)=en.ef(1);
fsp.sf(2)=fs.sf(1); fsp.ef(2)=fs.ef(end);
if (length(en.sf) > 1)
  fsp.sf(3)=en.sf(2); fsp.ef(3)=en.ef(2);
end

sd=sumdiff_pairs(d,p,fsp,ind.a,ind.b);

d_sd_p0=filter_scans(sd,fs,'p0');
d_sd_p0=filter_scans(d_sd_p0,en,'p0');

d_sd_p3=filter_scans(sd,fs,'p3');
d_sd_p3=filter_scans(d_sd_p3,en,'p3');


% CUTS WITH SIZE 1 x NRECEIVERS

% Correlation coeff over the focal plane post p3 filter

for i=1:cp.nrx
  
  % find the channels in this rx
  chs=intersect(find((p.rx)==i-1),ind.gl);

  % Exclude half-scans for which more than 10% of channels contain NaN's
  mapind = [];
  for k=1:length(fs.sf)
    fracnans=sum(isnan(sum(d_p3.mce0.data.fb(fs.sf(k):fs.ef(k),chs))))/length(chs);
    if fracnans<0.10
      mapind=[mapind fs.sf(k):fs.ef(k)];
    end
  end
  
  % find the corr coeff matrix for this rx
  cor=corrcoef(d_p3.mce0.data.fb(mapind,chs));
  
  % NaN out the diagonal
  for j=1:length(cor)
    cor(j,j)=NaN;
  end
  
  % reduce to per rx cut parameter
  cp.fp_cor(i)=nanmean(nanmean(cor));
end

% Focal plane temperature cut -
% which register to use depends on which experiment we are
% use presence of antenna.hk0/1/2 vs just .hk as test if we are Keck or B2
mapind=make_mapind(d,fs);

switch(expt)
  case 'keck'
    % we are Keck - one temp val per focal plane
    for r=unique(p.rx)'
        if ismember(r,[0,1,3])
          % this register selected on advice of Darren in email 14 May 2012
          tfpu(:,r+1)=d.antenna0.(['hk' num2str(r)]).ntd_temp(mapind,7);
        else
          % use other FPU NTD if rx2 or rx4 as per control NTD in schedlib
          tfpu(:,r+1)=d.antenna0.(['hk' num2str(r)]).ntd_temp(mapind,8);
        end
    end
  case 'bicep2'
    % we are B2
    tfpu=d.antenna0.hk.fast_temp(mapind,8);
  case 'bicep3'
    % We are BICEP3.  Many MCEs and one rx to rule them all.  Use the NTD on UC plate
    % that's not the one we're controlling on.
    tfpu=d.antenna0.hk0.ntd_temp(mapind,2);
end
tfpu=double(tfpu);

% The focal plane temperatures "ring" any time a frame is dropped by the
% mediator (presumably because of a GCP readout filter...??), but these are
% not physical temperature fluctuations. Therefore, remove a small chunk of
% data to keep them from affecting the tfpu cuts.
%
% Get a list of the skips spanning the scans block (which mapind spans)
b.sf = fs.sf(1);
b.ef = fs.ef(end);
skips = identify_skips(d, b);
% Then filter to keep just the filter errors and mediator dropped frames.
tmp = ismember(skips.name, {'filter failures','mediator dropped frames'});
skips = structcut(skips, tmp);
% If there are identified frame failures, NaN out segments of the data stream.
if numel(skips.sf) > 0
  % Increase the span by half the transfer function size on each size. This is
  % actually probably a bit of an overestimate (since tf.deconvkern includes a
  % low-pass filter and the MCE readout transfer function), but it's the best
  % we have available (and it's also only ~33 samples in Keck).
  kernlen = max(arrayfun(@(tf) length(tf.deconvkern), d.tf));
  skips.sf = skips.sf - floor(kernlen/2);
  skips.ef = skips.ef + floor(kernlen/2);
  % Then build an appropriate mask
  skipmask = false(length(d.t), 1);
  for ii=1:length(skips.sf)
    skipmask(skips.sf(ii):skips.ef(ii)) = true;
  end
  % Then apply NaNs to effectively disable the frame drop's effect on tfpu_std.
  tfpu(skipmask(mapind),:) = NaN;
  % Keep track if numel(skips.sf) > 0
  cp.skipmask = true;
end
if ~isfield(cp,'skipmask')
  cp.skipmask = false;
end

% take the mean and std of the fpu temps
cp.tfpu_mean=nanmean(tfpu);
cp.tfpu_std=nanstd(tfpu);
clear tfpu

% For Keck 2014+, the relative calibration (for atmospheric tau) has changed
% and we record a new quantity: the median that all channels are scaled to,
% derived from elnod analysis in power units. Report this as the "elnod_median"
% even though that isn't technically correct because it replaces the historic
% use for elnod_median.
if isfield(en,'gmed')
  cp.elnod_median = zeros(size(rxs));
  % The medians are stored by per-frequency naming in a structure fields, so
  % pull those out by dynamic field. Assign all receivers of the given
  % frequency.
  for freq=freqs
    fname = sprintf('med%d', freq);
    cp.elnod_median(rxfreqs==freq) = mean(en.gmed.(fname));
  end
% Improve beyond the historic case which defined
%    elnod_median = median(elnod_mean) = median(mean(elnod1,elnod2))
% but the relgain factors in cal_scans() are computed for each channel using a
% median which takes the form
%    elnod_median = mean(median(elnod1),median(elnod2))
% The historic case can be differentiated from this "fixed" case by absence or
% presence of the cp.elnod_median field, respectively.
else
  % Grab the gains and filter as is done in cal_scans().
  gains = en.g(:,:,2);

  % Get the year from the data itself, and like cal_scans(), only filter for
  % years >= 2014.
  ddate = datevec(utc2datenum([fs.t(1),fs.s(1)]));
  if ddate(1)>=2014
    gains(gains<100) = NaN;
  end

  for freq=freqs
    fname = sprintf('rgl%03d', freq);
    en_med = mean(nanmedian(gains(:,ind.(fname)), 2));
    cp.elnod_median(rxfreqs==freq) = en_med;
  end
end


% CUTS WITH SIZE NHALF-SCAN x NCHANNELS

% take scan std with p0 filtering
cp.fb_std_p0=scan_std(d_p0,fs);

% take scan std with p3 filtering
cp.fb_std_p3=scan_std(d_p3,fs);

% take un-relgain calibrated scan std
cp.fb_std_uncal=fs.std;

% take post sum diff scan std with p0 filtering
cp.fb_std_sd_p0=scan_std(d_sd_p0,fs);

% take scan std with p3 filtering
cp.fb_std_sd_p3=scan_std(d_sd_p3,fs);

% count number of nans in the half-scan
for i=1:length(fs.s)
  s=fs.sf(i); e=fs.ef(i);
  cp.fb_nancount(i,:)=sum(isnan(d.mce0.data.fb(s:e,:)));  
end

% following two are aberations - all other cut params have their
% "natural" number of columns, but these two are pre-expanded to be
% per channel

% flux jumps - if they occur, flag entire col
cp.is_fj_col=false(size(cp.fb_std_p0));
for mce=unique(p.mce)'
  for col=unique(p.mce_col)'
    rp=find((p.mce_col==col) & (p.mce==mce));
    for i=1:length(fs.s)
      s=fs.sf(i); e=fs.ef(i);
      if isfield(d.mce0.data,'numfj') && ~isempty(d.mce0.data.numfj) && nanmedian(nanmedian(d.mce0.rc1.data_mode))~=7
        isfj=logical(nansum(nansum(d.mce0.data.numfj(s:e,rp))));
      else
        isfj=0;
      end
      ismax=(nanmax(cp.fb_std_p0(i,rp),[],2))>1000;
      if isfield(dg,'fluxjumps')
        isfjcomp=logical(nansum(nansum(dg.fluxjumps(s:e,rp))));
      else
        isfjcomp=false;
        for j=1:length(rp)
          isfjcomp=isfjcomp | any(dg.fjlist(rp(j)).sf<=e & dg.fjlist(rp(j)).ef>=s);
        end
      end
      cp.is_fj_col(i,rp)=isfj | ismax | isfjcomp;
    end
  end
end

% flux jumps - if they occur, flag entire row
cp.is_fj_row=false(size(cp.fb_std_p0));
for mce=unique(p.mce)'
  for row=unique(p.mce_row)'
    rp=find((p.mce_row==row) & (p.mce==mce));
    for i=1:length(fs.s)
      s=fs.sf(i); e=fs.ef(i);
      if isfield(d.mce0.data,'numfj') && ~isempty(d.mce0.data.numfj) && nanmedian(nanmedian(d.mce0.rc1.data_mode))~=7
        isfj=logical(nansum(nansum(d.mce0.data.numfj(s:e,rp))));
      else
        isfj=0;
      end
      ismax=(nanmax(cp.fb_std_p0(i,rp),[],2))>1000;
      if isfield(dg,'fluxjumps')
        isfjcomp=logical(nansum(nansum(dg.fluxjumps(s:e,rp))));
      else
        isfjcomp=false;
        for j=1:length(rp)
          isfjcomp=isfjcomp | any(dg.fjlist(rp(j)).sf<=e & dg.fjlist(rp(j)).ef>=s);
        end
      end
      cp.is_fj_row(i,rp)=isfj | ismax | isfjcomp;
    end
  end
end



% CUTS WITH SIZE 1 x NCHANNELS

% mean of before/after elnod gains
cp.elnod_mean=nanmean(en.g(:,:,2),1);

% fractional delta of before/after elnod gains
% - as opacity changes with time this can differ without being a sign
% of pathological behavior
if(size(en.g,1)==1)
  % if only 1 elnod present we can't take the delta...
  cp.elnod_fracdel=nan(1,size(en.g,2));
else
  cp.elnod_fracdel=abs(en.g(2,:,2)-en.g(1,:,2))./cp.elnod_mean;
end

% take ratio-of-the-ratio of elnod pair gains before/after - this
% ratio better be close to one or else we are in danger of injecting
% false polarization signal into the maps
if(size(en.g,1)==1)
  cp.elnod_ab_ba=nan(1,size(en.g,2));
else
  cp.elnod_ab_ba=zeros(1,cp.nch);
  x=(en.g(1,ind.a,2)./en.g(1,ind.b,2))./(en.g(2,ind.a,2)./en.g(2,ind.b,2));
  x=abs(1-x);
  cp.elnod_ab_ba(1,ind.a)=x; cp.elnod_ab_ba(1,ind.b)=x;
end

% number of NaNs during el nods -- don't want any!
cp.elnod_nancount=squeeze(sum(en.g(:,:,4),1));


% elnod goodness of fit from elnod fitting routine
% take max of lead/trail
cp.elnod_gof=max(en.g(:, :, 3),[],1);

% elnod pair diff chi-square statistic
% take max of lead/trail
tmp=nan(2, cp.nch);
for i=1:length(en.s)
  npts=length(en.sf(i):en.ef(i));  
  tmp(i,ind.b)=nansum(d_sd_p0.mce0.data.fb(en.sf(i):en.ef(i),ind.b).^2) ...
    ./(nanstd(diff(d_sd_p0.mce0.data.fb(en.sf(i):en.ef(i),ind.b)))).^2/npts;  
  tmp(i,ind.a)=nansum(d_sd_p0.mce0.data.fb(en.sf(i):en.ef(i),ind.a).^2) ...
    ./(nanstd(diff(d_sd_p0.mce0.data.fb(en.sf(i):en.ef(i),ind.a)))).^2/npts;
end  
cp.elnod_chisq=max(tmp,[],1);

% elnod pair diff var statistic
tmp=nan(2, cp.nch);
for i=1:length(en.s)
  npts=length(en.sf(i):en.ef(i));
  tmp(i,ind.b)=nansum(d_sd_p0.mce0.data.fb(en.sf(i):en.ef(i),ind.b).^2)/npts;
  tmp(i,ind.a)=nansum(d_sd_p0.mce0.data.fb(en.sf(i):en.ef(i),ind.a).^2)/npts;
end
cp.elnod_var=max(tmp,[],1);

% In BICEP3 first-light data, no partial load curves.  Fill in NaN.
if ~isfield(lc,'g')
  lc.g=NaN*ones(2,cp.nch,4);
end

% detectors fractional resistance, rtes_frac=rtes/rnorm
% take mean before/after
x=lc.g(:,:,4)./lc.g(:,:,3);
cp.rtes_frac=nanmean(x,1);

% save rnorm,pjoule(end) for monitoring
% take mean before/after

cp.rnorm=nanmean(lc.g(:,:,3),1);
cp.pjoule=nanmean(lc.g(:,:,2),1);

% take noise level of average half-scan psd's in freq range 1.5-2Hz
% (avoids Keck PT line)
[psd,f]=mean_psd(d_sd_p0,fs);
fp=f>1.5&f<2;
cp.fb_wn_sd_p0=sqrt(mean(psd(fp,:)));
% take noise level of average half-scan psd's in freq range .1-.3Hz
fp=f>.1&f<.3;
cp.fb_1f_sd_p0=sqrt(mean(psd(fp,:)));
% power at the Pulse tube frequency for Keck
%if strcmp('keck',expt)
  % pulse frequency is 1.2 Hz.  Find the nearest PSD frequency
%  ptfreq=1.2;
%  freqind=find(abs(ptfreq-f)==min(abs(ptfreq-f)));
%  cp.pt_power=psd(freqind,:);
%end

%Calculate the NET as well
[psd,f]=mean_psd(d_sd_p3,fs);
fp=f>.1&f<1;
%NET=noise*(zenith temp)/(elnod gain)*(convert to uK)
cp.net=sqrt(mean(psd(fp,:)))./nanmean(nanmedian(en.g(:,ind.rgl,2),2),1)*12.7*1e6;
%copy pair diff into a.  this is per-pair NET
cp.net(ind.a)=cp.net(ind.b); 
%change 0 net to nan
cp.net(cp.net==0)=NaN;

% Calculate post-p3 skewness
mapind=make_mapind(d,fs);
cp.skewness=skewness(d_sd_p3.mce0.data.fb(mapind,:));

% Noise model depends on entire scanset being reasonable - no steps
% and nothing crazy inbetween scans
s=fs.sf(1); e=fs.ef(end);
cp.scanset_std=nanstd(sd.mce0.data.fb(s:e,:));

% The next three relate to destepping.  They are designed to be used in
% second-round cuts and will be propagated to neighbor channels according
% to logic in eval_round2_cuts.  In general, we want to cut channels that
% have too many flux jumps; channels that haven't been fully destepped;
% and neighbors of channels that have too many (and unevenly distrib-
% uted) jumps.

if isfield(dg,'njump')
  % Number of flux-jump steps found in each channel.
  % (Is it more than we can handle?)
  cp.num_fj=dg.njump(:,1)';
  % cp.num_fj(cp.num_fj>dg.maxstep)=Inf;

  % Number of times this channel should have been destepped.
  % (Was the destepping actually applied?)
  cp.num_destep=dg.njump(:,2)';
  cp.num_destep(dg.njump(:,2)>0 & dg.njump(:,3)==0)=Inf;

  % Check for constant raging.  If the channel is raging / flux jumping,
  % is it uniform throughout the scanset?  Or turning on and off?
  for i=1:length(dg.fjlist)
    clear tmp
    tmp(:,1)=[en.sf(1); dg.fjlist(i).ef];
    tmp(:,2)=[dg.fjlist(i).sf; max(en.ef(end),fsp.ef(2))];
    cp.max_fj_gap(i)=max(tmp(:,2)-tmp(:,1));
  end
else
  cp.num_fj=zeros(size(cp.scanset_std));
  cp.num_destep=zeros(size(cp.scanset_std));
  cp.max_fj_gap=zeros(size(cp.scanset_std));
end

% see 20130917_azcut post by JET

%construct az fixed signal
az=ground_signal(d_sd_p3,fs,ind.e);

%find places where mean az signal changes
for ch=ind.rglb
  
  %smooth with a rect. tophat
  xx=az.v(:,:,ch)';
  ker=[1,1,1;1,1,1;1,1,1];
  xx=conv2(xx,ker);
  
  %create a running mean
  for i=3:size(xx,1)-3   %we don't get to use the first and last halfscans
    mean_before_x=nanmean(xx(1:i-1,:),1);
    mean_after_x=nanmean(xx(i:end,:),1);

    %use var before and after as estimate of noise
    std_before_x=nanstd(xx(1:i-1,:),1);
    std_after_x=nanstd(xx(i:end,:),1);
    ee=sum((std_after_x-std_before_x).^2);
    
    %chisqaure
    x2=sum((mean_after_x-mean_before_x).^2)./ee;
   
    %store
    t4(i,ch)=x2;
  end
 
end

cp.satcom=zeros(1,length(ind.e));
val=nanmax(t4,[],1);
cp.satcom(ind.rglb)=val(ind.rglb);
cp.satcom(ind.rgla)=val(ind.rglb);

% CUTS WITH SIZE 1 x SOMETHING ELSE

% record TES bias values for each tile/col as they were at start of
% first field scan
cp.tes_bias=d.mce0.tes.bias(fs.s(1),:);

% az encoder differences
if(isfield(d.antenna0.pmac,'enc_az_diff'))
  cp.enc_az_diff=nanmedian(d.antenna0.pmac.enc_az_diff);
end

%cut on crazy az ranges
% Dropped frames also disturb the pointing information, so reuse the skip
% mask created to avoid cutting unphysical transients in tfpu_std.
az = d.pointing.hor.az;
if numel(skips.sf) > 0
  az(skipmask) = NaN;
end
cp.az_range=range(az(fs.sf(1):fs.ef(end)));
clear az

%NET per column
seq_col_list=unique(p.sequential_mce_col(isfinite(p.sequential_mce_col)));
for seq_col=seq_col_list'
    rp=find(p.sequential_mce_col==seq_col & strcmp(p.pol,'B'));
    if ~isempty(rp)
      cp.net_per_col(seq_col)=1/sqrt(nansum(1./cp.net(rp).^2));
    else
      cp.net_per_col(seq_col)=NaN;
    end
end

%NET by tile
seq_tile_list=unique(p.sequential_tile(isfinite(p.sequential_tile)));
for seq_tile=seq_tile_list'
  rp=find(p.sequential_tile==seq_tile);
  if ~isempty(rp)
    cp.net_per_tile(seq_tile)=1/sqrt(nansum(1./cp.net(rp).^2));
  else
    cp.net_per_tile(seq_tile)=NaN;
  end
end

%NET by Rx
for rx=unique(p.rx)'
  rp=find(p.rx==rx);
  cp.net_per_rx(rx+1)=1/sqrt(nansum(1./cp.net(rp).^2));
end


% CUTS WITH SIZE NHALF-SCAN x 1

% check fs indices make sense
cp.fs_inrange=fs.sf>0 & fs.ef>0 & fs.sf<fs.ef ...
    & fs.sf<=size(d.t,1) & fs.ef<=size(d.t,1);

% keep a record of the time of each half-scan
cp.fs_t=fs.t;

% record number of data points in each half-scans
cp.fs_npts=fs.ef-fs.sf+1;


% CUTS WITH SIZE NHALF-SCAN x NMCE

% check mce sync matching/ordering
for i=1:length(fs.s)
  s=fs.sf(i); e=fs.ef(i);
  x=double(d.mce0.syncBox.sampleNumber(s:e,:));
  % diff over time dimension - less than one indicates out of order
  % frames
  cp.syncsampnum_diff1(i,:)=sum(diff(x,[],1)<1);
  % diff over rx dimension - anything but zero indicates mismatched
  % frames
  cp.syncsampnum_diff2(i,:)=sum(bsxfun(@minus,x,mode(x,2))~=0);
end



return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

function az=ground_signal(d,fs,ind)
% az=ground_signal(d,fs,ind,psm)
% 
% extracted from ground_subtract
%
% ex.
% az=ground_signal(d,fs,ind.e);
% 


disp('ground_signal...');


nscanset=size(unique(fs.set),1);
nch=length(ind);

% for each scanset
for i=1:nscanset

  %find halfscans in that scanset
  q=find(fs.set == i);
  nhalfscans=length(q);
  npoints=fs.ef(q(1))-fs.sf(q(1))+1;
  
  %cut fs structure to only include this scanset.
  fs_tmp=structcut(fs,q);
  mapind=make_mapind(d,fs_tmp); 
  
  %take light bolo only in these halfscans
  if isfield(d.antenna0,'bolo')
    v=double(d.antenna0.bolo.mag(mapind,ind));
  else
    v=double(d.mce0.data.fb(mapind,ind));
  end
  v0=v;

  %mask out Nan in timestream
  numNaN=zeros(size(v));
  numNaN(isnan(v))=1;

  %reshape
  v=reshape(v,npoints,nhalfscans,nch);
  v0=reshape(v0,npoints,nhalfscans,nch);  

  %split left and right going
  l=find(fs_tmp.inc == 1);
  r=find(fs_tmp.inc == 0);
  vl=v(:,l,:);
  vr=v(:,r,:);
  vl0=v0(:,l,:);
  vr0=v0(:,r,:);
  
  
  %actual az coords
  az.az=reshape(d.pointing.hor.az(mapind),npoints,nhalfscans);
  az.azl=mean(az.az(:,l),2);
  az.azr=mean(az.az(:,r),2);


  %flip right going scans 
  vr0=flipdim(vr0,1);
  
  %put back in full strucutre
  v0(:,l,:)=vl0;
  v0(:,r,:)=vr0;
end

 %%%%%%%%%%%%%%%%%%%%%%%%
  %store info
  az.az=az.azl; %we use left going as the ordering in az
  az.v=v0;
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%
  
return
