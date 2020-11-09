function d=read_run(tag,tend,ch)
% d=read_run(tag)
%
% Read a chunk of data from arc files, massage a bit and return
% accepts either "tag" or time range in any format load_arc accepts
%
% e.g.
% d=read_run('20100103a')
% d=read_run('110119 06:38:53','110119 06:44:45')

disp('read_run...');

if exist('tend','var') && ~isempty(tend)
  tstart=tag;
  tend=tend;
else
  % find start and end times for this tag
  r=get_run_info;
  i=strmatch(tag,r.tag,'exact');
  if isempty(i)
    error(['ERROR in read_run: tag ' tag ' not found in list returned by get_run_info']);
  end
  tstart=r.tstart{i};
  tend=r.tend{i};
end

chanstr={};
% If the user specified a channel list, break the channels out
% by MCE and GCP index.
if exist('ch','var') && ~isempty(ch)
  dt=utc2datenum(tstart);
  [p,ind]=get_array_info(datestr(dt,'yyyymmdd'));
  if ~isfield(p,'mce')
    if isfield(p,'rx')
      p.mce=p.rx;
    else
      p.mce=zeros(size(p.r));
    end
  end
  for i=0:max(p.mce)
    chpermce=reshape(p.gcp(ch(p.mce(ch)==i)),1,[]);
    tmpfb=['mce' num2str(i) '.data.fb'];
    tmpfj=['mce' num2str(i) '.data.numfj'];
    tmp='(';
    if ~isempty(chpermce)
      tmp=[tmp num2str(chpermce+1,'%d,')];
      if tmp(end)==','
        tmp=tmp(1:(end-1));
      end
    end
    tmp=[tmp ')'];
    chanstr=[chanstr,{[tmpfb tmp]},{[tmpfj tmp]}];
  end
else
  ch=[];
end

% read from the archive using Walt's load_arc command.
% currently select all registers ('*'), but first
% specify empty channel lists for registers we don't
% want.  These guys will not be loaded into RAM.
d=load_arc('arc',tstart,tend,{'mce*.data.err()',chanstr{:},'*'});

% utc comes back from reader as two column mjd/secinday - make
% single col
d=make_utc_single_col(d);
  
% make string fields into strings
d.antenna0.tracker.source=cast(d.antenna0.tracker.source,'char');
d.antenna0.tracker.scan=cast(d.antenna0.tracker.scan,'char');

% Convenient to shorten name of slow and fast time
d.ts=d.antenna0.time.utcslow;
d.t=d.antenna0.time.utcfast;

% Find the date from the system time (because one of the kludges is
% to fix the GPS time)
indx=find(d.antenna0.frame.utc~=0);
if ~isempty(indx)
  [year,month,day]=datevec(d.antenna0.frame.utc(indx(1))+678942);
else
  error('Returned no data or all 0 in d.antenna0.frame.utc');
end

% Bad GPS time kludge
% Manufacturer says went bad on "Feb 13"
% B3 was fixed 20160406. Keck was fixed around 20160408:04:00
% fixing it when not broken shouldn't do harm if a fix is in fact
% possible
sched=year*1e4+month*1e2+day;
if(sched>=20160212&sched<20160409)
  % The GPS cards in both Keck and B3 were on the fritz during this
  % period and antenna0.time.utc is wrong. This means
  % reconstructed pointing will be wrong. Attempt to kludge on the
  % assumption that system time is still OK (which it ought to be
  % since the control computers do ntp sync to station timeserver).
  
  % test the offset between the last GPS timestamp in each frame
  % and the frame.utc time
  % Testing some random times on Feb-10 and Apr-08 gives offsets in
  % the region of 0.020 seconds
  % plot((d.antenna0.frame.utc-d.antenna0.time.utcfast(20:20:end))*86400)
  % see http://bicep.rc.fas.harvard.edu/bkcmb/analysis_logbook/analysis/20160412_gps_problem
  
  % take the delta between the samples within each frame and the end
  % of that frame
  sr=length(d.t)/length(d.ts);
  rel=cvec(repmat(d.antenna0.time.utcfast(sr:sr:end)',sr,1));
  del=d.antenna0.time.utcfast-rel;
  % expand the system clock utc
  sys=cvec(repmat(d.antenna0.frame.utc',sr,1));
  % add together and also apply an empirical correction based on
  % the offset observed when GPS is working
  d.t=sys+del-0.020/86400;
  % account for dropped samples where time goes to zero; preserve them
  d.t(d.antenna0.time.utcfast==0 | rel==0) = 0;
  % take the slow time from the system time also (not sure what
  % this is used for)
  d.ts=d.antenna0.frame.utc;
  % sometimes the system time has large 1 or 2 frame glitches. NaN them.
  % find glitches as anything more than 2x away from expected max time diff from median
  tsmed=nanmedian(d.ts); tmed=nanmedian(d.t);
  tsord=sort(d.ts(d.ts>0)); tord=sort(d.t(d.t>0));
  tsch=tsord(round(length(tsord)/2))-tsord(1);
  tch=tord(round(length(tord)/2))-tord(1);
  d.ts((d.ts>(tsmed+2*tsch))|(d.ts<(tsmed-2*tsch)))=NaN;
  d.t((d.t>(tmed+2*tch))|(d.t<(tmed-2*tch)))=NaN;
  
  warning('Attempting correction for bad GPS time in early 2016');
end

% Experiment-specific kludges
switch get_experiment_name()
  case 'bicep2'

    % Early 2011 BICEP2 data has no online NTD calibrations, so we apply them
    % here instead.  Use GCP thermal calibration source file for the curves.
    if((sched>=20110202)&&(sched<=20110315))
      warning('Correcting early 2011 B2 data to apply ntdcal');
      d=ntdcal(d,'aux_data/ntd_cal/TempConvert.cc');
    end

  case 'keck'

    % the sign of the data is reversed in Keck cf B2
    signflip=true;

    % another lame sign flip... one of the columns in K4 is wired reversely
    % during the 2012 season
    if(isfield(d,'mce3') && year==2012)
      k4col1signflip=true;
    end

    % another sign flip for 2013 through the end of 2014 when K4 contains D2,
    % which shouldn't have the sign flip
    if(isfield(d,'mce3') && year>=2013 && year<2015)
      k4d2signflip=true;
    end
end

% for Keck and BICEP3 there are multiple mce's - concatenate these
for m={'mce1','mce2','mce3','mce4','mce5'}
  if(isfield(d,m))
    d.mce0=structcat(2,[d.mce0,getfield(d,m{1})]);
    d=rmfield(d,m);
  end
end
if isempty(ch)
  ch=1:size(d.mce0.data.fb,2);
end

% flip the sign - see above
if exist('signflip','var') && signflip
  d.mce0.data.fb=-d.mce0.data.fb;
end

if exist('k4col1signflip','var') && k4col1signflip
  k4col1=1618:1650;  %hardwired in since we haven't read in array info
  ck4col1=ismember(ch,k4col1);
  d.mce0.data.fb(:,ck4col1)=-d.mce0.data.fb(:,ck4col1);
end

if exist('k4d2signflip','var') && k4d2signflip
  k4=1585:2112;  %hardwired in since we haven't read in array info
  ck4=ismember(ch,k4);
  d.mce0.data.fb(:,ck4)=-d.mce0.data.fb(:,ck4);
end

% shrink mce0.data.numfj
d.mce0.data.numfj=int8(d.mce0.data.numfj);

return

%%%%%%%%%%%%%%%%%%%%%%%%
%  function d=make_utc_single_col(d)
%  this function migrated to make_utc_single_col.m
