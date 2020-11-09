function [tod,mapopt]=prepare_pairmap_tod(tag,mapopt)
% [tod,mapopt]=prepare_pairmap_tods(tag,mapopt)
%
% Prepares detector timestreams for map binning by pair summing/differencing,
% applying time-stream filters (poly/scan-sync), computes weights, etc. Also
% generates various diagnostic and accounting metadata and stores to mapopt.
%
% Extracted out of reduc_makepairmaps and reduc_matrix and merged to share a
% single implementation.
%
% INPUTS
%
%   tag       Tag name or struct of tag data.
%
%   mapopt    Map making option structure, as from get_default_mapopt().
%
% OUTPUTS
%
%   tod       Struct of timestream data (.d) along with other descriptor
%             data structures (i.e. .en, .fs, .p, .ind, .simopt, etc).
%
%   mapopt    As input with additional diagnostic fields added.
%
% EXAMPLE
%
%   tag = '20170507E08_dk023';
%   mapopt.sernum = '0000real';
%   mapopt = get_default_mapopt(mapopt);
%   [tod,mapopt] = prepare_pairmap_tod(tag, mapopt);
%   ...
%

tic();

if isstruct(tag)
  % data was passed directly (e.g. from reduc_makesim)
  todstruct = tag; clear tag
  liftvars(todstruct); clear todstruct
  mapopt.filename = 'none';
else
  % read in tod
  if logical(strfind(mapopt.sernum,'real'))
    if(~isfield(mapopt,'altcal'))
      filename=['data/real/' tag(1:6) '/' tag '_tod.mat'];
    else
      filename=['data/real/' tag(1:6) '/' tag '_tod_altcal.mat'];
    end
    if(isfield(mapopt,'powerunits'))
      if mapopt.powerunits
        filename=['data/real/' tag(1:6) '/' tag '_tod_pu.mat'];
      end
    end
  else
    filename=['scratch/' mapopt.sernum(1:4) '/' mapopt.sernum(5:end) '/' tag '_tod.mat'];
  end
  fprintf(1,'\nloading tod from %s...\n',filename);
  load(filename);
  mapopt.filename=filename;
end

mapopt.tag=tag;

% Flag of whether making un-summed/diffed raw detector maps
isabmap = isfield(mapopt, 'abmap') && mapopt.abmap;
% Flag of whether making a simulation
issim = exist('simopt','var');
% Flag of whether making a matrix
isobsmat = isfield(mapopt, 'sigmaptype') && isfield(mapopt, 'siginterp');
% Flag of whether making a reweight sim/matmrix
issubsim = issim && isfield(simopt,'subsim') && simopt.subsim;
issubsim = issubsim | (isobsmat && isfield(mapopt, 'flavor_subset') ...
    && mapopt.flavor_subset);

% remove large register blocks not needed here to save memory
d.mce0.data.err=[];
d.mce0.data.raw=[];
d.mce0.data.numfj=[];
d.antenna0.hk.fast_voltage=[];

% read psm - point source mask
psmfilename=sprintf('pntsrcmask/%s_pntsrcmask',tag);
if exist([psmfilename,'.mat'],'file') && mapopt.usepsm==1
  fprintf(1, 'using point source mask %s\n', psmfilename);
  load(psmfilename); mapopt.psmfilename=psmfilename;
end
% if no psm exists or we have been told not to use it set it null
if ~exist('psm','var') || mapopt.usepsm==0
  psm=false(size(d.mce0.data.fb));
end

toc() % time to read tag

tic()

% if requested, perform additional lowpass filtering on the time streams
if isfield(mapopt,'extralpf')
  tic()
  % determine if this is a tag that should be lowpass filtered
  ri=get_run_info(tag);
  % if this tag comes after tagstart, proceed with the additional low pass
  % filtering
  if datenum(ri.tstart{1},'dd-mmm-yyyy:HH:MM:SS') > mapopt.extralpf.dtnumstart
    % This attempts to follow how the low pass filtering is done during
    % reduc_initial Except we want to know the sample rate so we can scale the
    % cutoff frequency correctly
    [sampratio,samprate]=get_sampratio_rate(d);
    N=32;
    nyq=samprate/2;
    fpnt=mapopt.extralpf.cofreq(1)/nyq;
    spnt=mapopt.extralpf.cofreq(2)/nyq;
    hlpf=firpm(N,[0 fpnt spnt 1],[1 1 0 0])';
    % record the extra low pass filter that's applied
    mapopt.extralpf.lpf=hlpf;
    mapopt.extralpf.info=sprintf(['low pass filter firpm, fpnt = %0.4f, ' ...
        'spnt = %0.4f with %d taps = [ %s ]'], fpnt, spnt, length(hlpf), ...
        sprintf('%0.3f, ', hlpf));

    % convolve timestream with new low pass filter
    s=fs.sf(1)-N/2-1; e=fs.ef(end)+N/2+1;
    disp('Performing Additional Low Pass Filter')
    d.mce0.data.fb(s:e,:)=convn(d.mce0.data.fb(s:e,:),hlpf,'same');
    mapopt.extralpf.applied=true;
  else
    mapopt.extralpf.applied=false;
  end
  toc()
end

tic()

% Get array info and index arrays
[p,ind]=get_array_info(tag,mapopt.beamcen,mapopt.chi,[],[],[],mapopt.epsilon);

% Modify detector pol angle due to waveplates for this date/time
p=do_hwp_rot(d.t(1),p);

% substitute dark squids for light channels in some fashion
if any(mapopt.mapdarksquids)
  if isobsmat
    error('mapping dark squids is not supported for matrices')
  end
  [d,simopt]=map_squids(d,p,ind,mapopt.mapdarksquids);
  % need to switch to simulation-like behavior
  issim = true;
end

if any(mapopt.mapfpntds)
  if isobsmat
    error('mapping focal plane ntds is not supported for matrices')
  end
  [d,simopt]=map_fpntds(d,fs,p,ind,mapopt.mapfpntds);
  % need to switch to simulation-like behavior
  issim = true;
end

% if we are real record p info in mapopt for the record - this might
% be useful for diagnostics on the pointing/polpar etc if we ever
% decide to vary that as a function of time - on the hand it is quite big
if ~issim
  mapopt.p=p;
end

switch mapopt.chi
  case 'assim'
    p.chi=simopt.p.chi;
    p.chi_thetaref=simopt.p.chi_thetaref;
end

switch mapopt.epsilon
  case 'assim'
    p.epsilon=simopt.p.epsilon;
end

switch mapopt.beamcen
  case 'assim'
    p.r=simopt.p.r;
    p.theta=simopt.p.theta;
end

% calculate the cuts, etc

% For reweight simulations, do *not* load real pairmaps or apply any cuts
if issubsim
  disp('subsim in effect - forcing weight to 0');
  mapopt.weight=0;
else
  % If we are real data and not a reweight matrix, calculate the cuts
  if ~issim
    % calculate the cutparams
    cp=calc_cutparams(d,fs,en,lc,p,ind,dg);
    % evaluate the cuts
    [c,cp]=eval_round1_cuts(cp,mapopt.cut,p,ind);
    % combine to final mask array (with stat printout)
    c.overall=combine_cuts(c,ind);
    % record this info for re-use when pairmapping sims and when
    % coadding real and sim
    mapopt.c.cp=cp; % cut parameters
    mapopt.c.c=c;   % cuts (logical masks)

  % For simulations or a non-reweight matrix, load cuts from real pairmaps
  else
    % load the cuts, weights and var from real pairmaps file
    % - warning - this means we can't run sims until real data has been
    % processed with the exact same set of options
    if isfield(mapopt,'realpairmapopt')
      realpairmapopt=mapopt;
      fld=fieldnames(mapopt.realpairmapopt);
      for k=1:numel(fld)
        realpairmapopt.(fld{k})=mapopt.realpairmapopt.(fld{k});
      end
    else
      realpairmapopt=mapopt;
    end
    [pathstr,filename]=fileparts(get_pairmap_filename(tag,realpairmapopt));
    filename=fullfile(mapopt.realpairmapset,filename);
    fprintf(1,'loading pre-calc cuts and var/weight from file %s\n',filename);
    x=load(filename,'mapopt');
    c=x.mapopt.c.c;
    w=x.mapopt.hs.w;
    v=x.mapopt.hs.s.^2;
    clear x;
  end

  % apply the cuts - nan out the bad half-scan/channels
  d=apply_cuts(d,fs,c.overall);
end

% Can weight by anything we want. This rarely used option weights by
% reciprocal of the pre-sum/diff, pre-filter half scan variance.
if mapopt.weight==1 && ~issim && ~issubsim
  % calc var from d as normal
  v=scan_std(d,fs).^2;
  % weight sum/diff equally using sum of un-diff vars
  v(:,ind.la)=v(:,ind.la)+v(:,ind.lb); % sum A and B var
  v(:,ind.lb)=v(:,ind.la);             % copy sum to B
  w=1./v;                              % conv to weight
end

if ~isabmap
  % take sum/diff of A/B pairs
  % store sum in A, diff in B
  % do for matrix to get NaNs in correct place
  [d,p]=sumdiff_pairs(d,p,fs,ind.a,ind.b);
end

% do secondary residual relgain correction using atmosphere
% (doesn't this just mean we should have used the atmosphere in the
% first place?)
if mapopt.resrelgain == 1
  if isobsmat
    error('preforming resrelgain is not supported for matrices')
  end
  d=remove_relgain2(d,fs,ind.la);
end

p.chi_mean=zeros(size(p.chi));
p.delalpha=zeros(size(p.chi));
p.rr=zeros(size(p.chi));
p.beta=zeros(size(p.chi));
a=ind.a; b=ind.b;

% calc alpha_horn and delta alpha_horn values per pair as specified by Ken Ganga
p.chi_mean(a)=mean([p.chi(a),p.chi(b)],2)-45;
p.delalpha(a)=-(p.chi(b)-p.chi(a)-90);
p.chi_mean(b)=p.chi_mean(a);
p.delalpha(b)=p.delalpha(a);

% calc r and beta parameters for diff data as specified by Ken Ganga
p.gamma=(1-p.epsilon)./(1+p.epsilon);
p.rr(b)=0.5*sqrt((p.gamma(a)+p.gamma(b)).^2.*cosd(p.delalpha(a)).^2 + ...
    (p.gamma(a)-p.gamma(b)).^2.*sind(p.delalpha(a)).^2);
p.beta(b)=0.5*atan2((p.gamma(a)-p.gamma(b)).*sind(p.delalpha(a)), ...
    (p.gamma(a)+p.gamma(b)).*cosd(p.delalpha(a)))*180/pi;

% not sure how to handle r' and beta' in original Ganga doc so ignored...

% take a copy tod pre-filter
d0 = d;
% filter TODs iff making maps or non-reweight matrices
if ~isobsmat || (isobsmat && ~issubsim)
  % filter half scans
  if ~iscell(mapopt.filt)
    % filter sum/diff the same
    [d,fc]=filter_scans(d,fs,mapopt.filt,ind.e,psm,mapopt.filttype);
    mapopt.filtc = fc;
    clear fc
  else
    % filter sum/diff differently
    [d,fc{1}]=filter_scans(d,fs,mapopt.filt{1},ind.a,psm,mapopt.filttype);
    [d,fc{2}]=filter_scans(d,fs,mapopt.filt{2},ind.b,psm,mapopt.filttype);
    if isempty(fc{1})
      mapopt.filtc = [];
    else
      sz = size(fc{1});
      sz(2) = length(ind.e);
      mapopt.filtc = zeros(sz);
      mapopt.filtc(:,ind.a,:) = fc{1};
      mapopt.filtc(:,ind.b,:) = fc{2};
      clear sz
    end
    clear fc
  end
  % keep a record of what was subtracted over and above the mean
  e=filter_scans(d0,fs,'p0',ind.e,psm,mapopt.filttype);
  d.psub=e.mce0.data.fb-d.mce0.data.fb;
  clear e
  % only keep the original for matrix-applied filtering later
  if ~isobsmat
    clear d0
  end

  % take a copy of pre-gsub tod
  fb_ngsub = d.mce0.data.fb;
  if mapopt.gs==1
    % do the scan sync removal if requested
    d=ground_subtract(d,fs,ind.l,psm);
  end
  % keep a record of what was subtracted
  d.gsub=fb_ngsub-d.mce0.data.fb;
  clear fb_ngsub

  % apply correction for imperfect polarization efficiency
  % this is in effect an additional cal applied only to difference data
  if ~isabmap
    for k=ind.lb
      d.mce0.data.fb(:,k)=d.mce0.data.fb(:,k)/p.rr(k);
    end
  end
end

% take variance post filter, sum/diff, efficiency correction -
% resulting variance maps should reflect pixel scatter in
% corresponding signal maps

% for real data we calc var from d - var is used to construct
% expected noise maps regardless of how weighting is being done
if ~issim
  v=scan_std(d,fs).^2;

else % issim == true
  % when doing substitution sims we set variance uniform
  if issubsim
    v=ones(length(fs.s),length(ind.e));

  % for noise only sim we also calc var from d - the
  % resulting variance maps need to be indistinguishable from the
  % real ones or the noise model is not working.
  elseif ~strcmp(simopt.noise,'none')
    v=scan_std(d,fs).^2;
  end
  % (for signal only sims v remains the one from real data read in
  % above.)
end

% option 0 - weight uniformly
if mapopt.weight==0
  w=ones(size(v));
end

% option 2 - weight by post filter sum/diff half scan 1/var
if mapopt.weight==2 && ~issim
  w=1./v;
end

% option 3 - same as 2 except weight is calculated over all
% half-scans. This much reduces the sample variance issue which JK
% is obsessed with but places stringent requirements on noise
% stationarity which only a small fraction of pair-sum data probably
% meets. However since the noise model will not work for such
% non-stationary data possibly that objection is moot.
if mapopt.weight==3 && ~issim
  mapind=make_mapind(d,fs);
  w=1./nanvar(d.mce0.data.fb(mapind,:));
  % expand to usual size
  w=repmat(w,[length(fs.sf),1]);
end

% Catch infinite weights - this might legitimately happen although
% the associated data had better also be NaN...
w(~isfinite(w))=NaN;

% Record half-scan stats of the data as actually binned into map for
% diagnostics and for use when mapping sims. Note that this has cut
% masking applied as opposed to the similar values in cut-parameters structure cp.
mapopt.hs.s=sqrt(v);
mapopt.hs.w=w;
% also calculate and record the max value occuring in each scan
% - values large compared to std should not occur
mapopt.hs.m=scan_max(d,fs);

% expand out var and weight values over each half scan
d.v=exp_scan_val(d,fs,v,~isobsmat);
d.w=exp_scan_val(d,fs,w,~isobsmat);

% histogram the "deviation" (tod/std) for each channel
mapind=make_mapind(d,fs);
for j=ind.e
  [bc,n]=hfill(d.mce0.data.fb(mapind,j)./sqrt(d.v(mapind,j)),100,-10,10);
  mapopt.devhist.bc=bc;
  mapopt.devhist.n(j,:)=n;
end

% keep some additional records to help determine which part of the
% sky was being scanned during this scanset
mapopt.traj.source=d.antenna0.tracker.source(fs.s(1),:);
mapopt.traj.scan=d.antenna0.tracker.scan(fs.s(1),:);
mapopt.traj.eloff=d.antenna0.tracker.horiz_off(fs.s(1),2)*180/pi;
dk=d.antenna0.tracker.horiz_off(fs.s(1),3)*180/pi;
dk(dk<0)=dk(dk<0)+360; mapopt.traj.dk=dk;
% Add information about the RA center to uniquely identify tags after starting
% the double-length scanset schedules.
%
% Same as calculated in reduc_initial() to create d.pointing.
lat = median(double(d.antenna0.tracker.siteActual(:,2))/3.6e6);
lon = median(double(d.antenna0.tracker.siteActual(:,1))/3.6e6);
% Translate az/el field center to ra/dec.
elcen = d.pointing.hor.el(fs.sf(1));
azcen = d.pointing.hor.az(fs.sf(1)) - d.azoff(fs.sf(1));
[mapopt.traj.racen,dum] = azel2radec(azcen, elcen, d.t(fs.sf(1)), lat, lon);

toc() % time to setup data

% Restore unfiltered data when making matrices. Filtering matrix will be
% applied in reduc_matrix
if isobsmat
  d.mce0.data.fb = d0.mce0.data.fb;
  clear d0
end

% always available:
tod.tag = tag;
tod.p = p;
tod.ind = ind;
tod.d = d;
tod.en = en;
tod.fsb = fsb;
tod.fs = fs;
tod.psm = psm;

% conditionally available:
if ~issim
  % only in real data
  tod.lc = lc;
  tod.dg = dg;
end
if issim
  % only for simulations
  tod.simopt = simopt;
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%

function [d,simopt]=map_squids(d,p,ind,g)
% Assign dark squids to light channels. Return a simopt to trick accumulate into using
% real map cuts and weights, even though this is not a sim.
simopt.noise='none';
simopt.rlz=0;

% Assign timestream of each column's dark squid to the A channel and zero the B
% channels. Columns 8 and 14 (zero indexed) appear to have bad squids
% (see June 3, 2010 posting from RWO) so just use squid before.
for k=1:numel(ind.o)
  col=p.mce_col(ind.o(k));
  detind=find(p.mce_col==col);

  if col==8 | col==14
    colind=ind.o(k-1);
  else
    colind=ind.o(k);
  end

  d.mce0.data.fb(:,detind)=repmat(d.mce0.data.fb(:,colind),[1,numel(detind)]);
end

% Zero all B channels
d.mce0.data.fb(:,ind.b)=0;

% Multiply by gains
if numel(g)>1
  for k=ind.l
    d.mce0.data.fb(:,k)=d.mce0.data.fb(:,k)*g(k);
  end
end

end % function map_squids

%%%%%%%%%%%%%%%%%%%%%%%%%

function [d,simopt]=map_fpntds(d,fs,p,ind,g)
% Assign focal plane temps to all channels and mutiply each channel by a number meant
% to represent thermal response. Return a simopt to trick accumulate into using
% real map cuts and weights, even though this is not a sim.
simopt.noise='none';
simopt.rlz=0;

% Average tiles 2 and 3 NTDs
t_nK = 1e9*nanmean(d.antenna0.hk.fast_temp(:,[1,2]),2);

% Remove mean
mapind=make_mapind(d,fs);
t_nK = t_nK - mean(t_nK(mapind));

% Place into all light channels
d.mce0.data.fb(:)=NaN;
d.mce0.data.fb(mapind,ind.l)=repmat(t_nK(mapind),[1,numel(ind.l)]);

% Multiply by thermal responses (units per nK)
if numel(g)>1
  for k=ind.l
    d.mce0.data.fb(:,k)=d.mce0.data.fb(:,k)*g(k);
  end
end

end % function map_fpntds

