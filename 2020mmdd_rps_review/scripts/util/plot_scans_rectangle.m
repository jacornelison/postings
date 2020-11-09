% plot_scans_rectangle - plot timestreams from a
% scanset in rectangular form, az x half-scans.
%   plot_scans_rectangle(tag,ch)
%
% Can pass options as parameter-value pairs:
%
%   plot_scans_rectangle(...,'sumdiff',true)
%      do pair sum/diff (default=true)
%   plot_scans_rectangle(...,'raw',true)
%      work from arc files (default=false)
%   plot_scans_rectangle(...,'time',true)
%      use UTC on x axis (default=false)
%   plot_scans_rectangle(...,'feed_az',true)
%      reckon per-feed az on y axis (default=false)
%   plot_scans_rectangle(...,'filt','p2')
%      filter to use (default='p2')
%   plot_scans_rectangle(...,'caxis',[-0.5 0.5])
%      color stretch (default=[-0.5 0.5])
%   plot_scans_rectangle(...,'ukpv',3250)
%      abscal number (default none)
%   plot_scans_rectangle(...,'savedir','~/posting/')
%      directory to save plots
%   plot_scans_rectangle(...,'newfig',true)
%      new figure window for each channel (default false)      
function h=plot_scans_rectangle(tag,ch,varargin)

if nargin<2 || isempty(ch)
  ch=[];
end

S=[];
S.sumdiff=true;
S.raw=false;
S.time=false;
S.feed_az=false;
S.filt='p2';
S.caxis=[-0.5 0.5];
S.ukpv=[];
S.savedir='';
S.newfig=false;
S.hsmean=false;

if length(varargin)==1 && isstruct(varargin{1})
  S=varargin{1};
else
  for i=1:2:length(varargin)
    S.(varargin{i})=varargin{i+1};
  end
end

if iscell(tag)
  for i=1:length(tag)
    plot_scans_rectangle(tag{i},ch,S);
  end
  return
end

[p,ind]=get_array_info(tag);
if isempty(ch)
  ch=ind.e;
end

if S.raw
  d=read_run(tag);
  [sampratio,samprate]=get_sampratio_rate(d);

  % calc these after lowpass...
  % determine nominal az/el and az/el offset of center pixel at each
  % time step (adds d.az_off,el_off and rawpoint.hor)
  d=arcvar_to_azel(d);

  % fetch the pointing model parameters
  % dropped samples throw off mean(d.t)
  pm=get_pointing_model(mean(d.t(d.t~=0)));
  % apply the inverse pointing model to generate d.pointing.hor
  d=invpointing_model(d,pm);

  % do the transform from az/el to ra/dec 
  %  - method and functions stolen from Cynthia bcn code point.h 
  % (lat/lon stored in milli arcsec)
  lat=median(double(d.antenna0.tracker.siteActual(:,2))/3.6e6);
  lon=median(double(d.antenna0.tracker.siteActual(:,1))/3.6e6);
  [d.pointing.cel.ra,d.pointing.cel.dec]=azel2radec(...
      d.pointing.hor.az,d.pointing.hor.el,d.t,lat,lon);
  d.pointing.cel.dk=d.pointing.hor.dk-parallactic(...
      d.pointing.hor.az,d.pointing.hor.el,d.t,lat,lon);

  % find start,end of field scan half-scans
  fs=find_scans(d.antenna0.tracker.scan_off(:,1),bitand(d.array.frame.features,2^1),sampratio);
  fs.t=d.t(fs.sf);
  if(~isempty(fs.s))
    % get scan info
    fs=get_scan_info(d.antenna0.tracker.scan(fs.s,:),fs);

    % tweak sf and ef such that half-scan cover only specifed scan throw
    % and all have exact same length
    % the 1.13 is turnfrac parameter used by Cynthia in BICEP1 analysis
    fs=tweak_field_scans(d,fs,1.13);
  end
  % make fsb (field-scan-blocks) structure
  fsb=find_blk(bitand(d.array.frame.features,2^1),sampratio);
  fsb.t=d.t(fsb.sf);

  % use fsb to fill in fs.set field Denis introduced
  for i=1:length(fsb.s)
    fs.set(fs.t>d.t(fsb.sf(i))&fs.t<d.t(fsb.ef(i)),1)=i;
  end

  % d.array.frame.features goes away before the last scan of a block is
  % actually finished. Fix this up by forcing block end to match the end
  % of the last scan
  for i=1:length(fsb.s)
    q=find(fs.set==i);
    fsb.e(i)=fs.e(q(end));
    fsb.ef(i)=fs.ef(q(end));
  end
else
  load(fullfile('data','real',tag(1:6),[tag '_tod.mat']));
end

if S.sumdiff
  % Poor man's pair sum/diff (no relgain unless it's already done in TOD file)
  d=sumdiff_pairs(d,p,fs,ind.a,ind.b);
end

% p2 subtract field scans
d=filter_scans(d,fs,S.filt,ind.e);

cc=false(size(d.t));
for i=1:length(fs.sf)
  cc(fs.sf(i):fs.ef(i))=true;
end

for jj=1:length(ch)

  chnum=ch(jj);

  xx=zeros(fs.ef(1)-fs.sf(1)+1,length(fs.sf));
  for i=1:length(fs.sf)
    xx(:,i)=d.mce0.data.fb(fs.sf(i):fs.ef(i),chnum);
    if mod(i,2)==0
      xx(:,i)=xx(end:-1:1,i);
    end
  end
  if ~isempty(S.ukpv)
    xx=xx*S.ukpv;
  end
  if jj==1 || S.newfig
    h=figure;
    setwinsize(gcf,700,400);
  end
  if S.feed_az
    % Reckon per-feed pointing
    % see http://bicep0.caltech.edu/~spuder/analysis_logbook/analysis/20130430_reckon/reckon.html
   [el,az] = reckon(d.pointing.hor.el,d.pointing.hor.az,p.r(chnum),90-p.theta(chnum)+d.pointing.hor.dk);
  else
    az=d.pointing.hor.az;
  end
  while max(az)<0
    az=az+360;
  end
  if ~S.time
    imagesc(1:length(fs.sf),az(fs.sf(1):fs.ef(1)),xx);
    xlabel('half-scan');
  else
    [yr mo dy hr mn sc]=mjd2date(fs.t);
    tt=datenum(yr,mo,dy,hr,mn,sc);
    imagesc(tt,az(fs.sf(1):fs.ef(1)),xx);
    xlabel('time');
    datetick('x');
    xlim([tt(1) tt(end)]);
  end 
  if S.feed_az
    ylabel('feed az / deg');
  else
    ylabel('boresight az / deg');
  end
  set(gca,'ydir','normal');
  colormap gray;
  if length(S.caxis)==1
    caxis([-1 1]*S.caxis);
  else
    caxis(S.caxis);
  end
  colorbar;
  title([strrep(tag,'_','\_') ' chan ' num2str(chnum)]);
  COLCOL={'red','green','cyan','magenta','yellow'};
  if exist('utc','var') && ~isempty(utc)
    if ~iscell(utc)
      utc={utc};
    end
    for i=1:length(utc)
      dn=utc2datenum(utc{i});
      [yr mo dy hr mn sc]=datevec(dn);
      mjd=date2mjd(yr,mo,dy,hr,mn,sc);
      hsnum=interp1(fs.t,1:length(fs.t),mjd);
      hold on;
      line(hsnum*[1 1],ylim,'color',COLCOL{i},'linestyle','--','linewidth',2);
      % text(hsnum,mean(ylim),utc{i},'color',COLCOL{i},'rotation',90,'horizontalalignment','center','verticalalignment','bottom','fontsize',14,'fontweight','bold');
    end
  end
  if ~isempty(S.savedir)
    if S.sumdiff
      fname=fullfile(S.savedir,[tag '-diff' num2str(chnum) '.png']);
    else
      fname=fullfile(S.savedir,[tag '-chan' num2str(chnum) '.png']);
    end
    mkpng(fname);
  elseif length(ch)>1 
    pause
  end
end
return
