function [d,ds]=deskip_scans(d,fsb)
% [d,ds]=deskip_scans(d,fsb)
%
% Deal with BICEP2/Keck skipped samples.
%
% These affect the antenna layer, and are usually single
% missing samples in antenna.*.  In some cases, there may
% be strings of skipped samples.
%
% Currently, we deal with skips by interpolating fast
% registers in antenna0.  There must be appropriate cuts
% to remove severe skips, or merging failures, when
% making maps!
%
% The 'ds' structure contains information about samples
% that have been interpolated, or larger skips that have
% not been interpolated.

disp('deskip... deal with BBCPCI missed samples for BICEP2/Keck');

% Get sync box number (data sample identifier)
% as recorded by antenna layer and MCE
syncant=d.antenna0.syncBox.sampleNumber;
syncmce=d.mce0.syncBox.sampleNumber(:,1);

% Identify a skip as any instance when they don't match
skipmask=syncant~=syncmce;

% Sample index
sampnum=1:length(skipmask);
ds.samp=sampnum(skipmask);
ds.interp_samp=[];

% Select registers to be deskipped.
% Currently: all fast regs under antenna0 except sync
fastregs=[];
brds=fieldnames(d.antenna0);
for i=1:length(brds)
  if strcmp(brds{i},'syncBox')
    continue
  end
  regs=fieldnames(d.antenna0.(brds{i}));
  for j=1:length(regs)
    if size(d.antenna0.(brds{i}).(regs{j}),1)==length(syncant)
      k=length(fastregs)+1;
      fastregs(k).map='antenna0';
      fastregs(k).brd=brds{i};
      fastregs(k).reg=regs{j};
      fastregs(k).nch=size(d.antenna0.(brds{i}).(regs{j}),2);
    end
  end
end

% for each blk
for i=1:length(fsb.sf)

  % find the scan block
  s=fsb.sf(i); e=fsb.ef(i);

  % find the skipped samples within the block
  skiplist=find(skipmask(s:e));

  % loop over skipped samples
  for j=1:length(skiplist)
    idx=(s-1)+skiplist(j);

    % Find nearby non-skip samples to use for interpolation
    % Start with all samples within +-10 ticks
    dist=10;
    goodmask=false(size(sampnum));
    goodmask(max(idx-dist,1) : min(idx+dist,sampnum(end)))=true;
    % and remove skips
    goodmask(skipmask)=false;

    % If not enough good neighbors, give up
    if sum(goodmask)<dist
      continue
    end

    ds.interp_samp=[ds.interp_samp; idx];

    % Loop over selected fast regs, interpolate each.
    for k=1:length(fastregs)
      d.(fastregs(k).map).(fastregs(k).brd).(fastregs(k).reg)(idx,:) ...
        =interp1(sampnum(goodmask), ...
        double(d.(fastregs(k).map).(fastregs(k).brd).(fastregs(k).reg)(goodmask,:)), ...
        idx,'linear','extrap');
    end 

    % Also interpolate d.t in the same way.
    d.t(idx,:) = interp1(sampnum(goodmask),double(d.t(goodmask,:)),idx,'linear','extrap');

  end
end

ds.bad_skip_samp=setdiff(ds.samp,ds.interp_samp);

return

