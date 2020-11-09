function [d,dg]=deglitch(d,sc,ch,p)
% [d,dg]=deglitch(d,sc,ch,p)
%
% This is a somewhat complicated deglitch evolved for use in BICEP2.
% It identifies four types of glitch:
%   1. Simple transients, using an 8-"grass unit" criterion
%   2. Large steps in feedback
%   3. Steps in the numfj register
%   4. Common-mode polarized small steps.
%
% if d.tf exists, then use the deconv kernel to determine the number of
% of points around glitches to cut, otherwise cut 1 sec of data
%
% The simple transients are replaced with NaN.  The large steps,
% numfj steps, and common-mode polarized steps have their offsets
% removed, unless a channel has > 10 such steps per scanset.  The
% channels in the same column and same-row/same-readout-card are
% also destepped.
%
% It's possible to find out what parts of the data were removed
% by the deglitcher by looking for NaNs in the TES time streams.
% The dg output contains information about the large steps and
% destepping in each channel.

disp('deglitch... (Spuder style - replace glitches with NaN, destep jumps)')

% Keep a record of glitches / jumps / desteps
%   1. njump(:,1) : number of large steps in this channel
%   2. njump(:,2) : number of (effective) times this channel should be destepped
%   3. njump(:,3) : number of times this channel actually was destepped
dg.njump=zeros(length(p.gcp),3);
%   4. fj(ch).sf/ef : list of large steps by channel
dg.fjlist(length(p.gcp)).sf=[];
dg.fjlist(length(p.gcp)).ef=[];
dg.maxstep=10;

% for each blk
tmp_ch=ch; clear ch

for i=1:length(sc.sf)
  % find the scan
  s=sc.sf(i); e=sc.ef(i);
  % find sample rate
  samprate=round(1/nanmedian(diff(d.t(s:e))*3600*24));

  % for each mce (same as for each rx in Keck)
  for imce=1+(0:max(p.mce))  
    mcech=tmp_ch(1+p.mce==imce);

    % get the tod for required channels
    x=d.mce0.data.fb(s:e,mcech);
    xj=d.mce0.data.numfj(s:e,mcech);

    % Set aside unmodified version
    z=x;

    % If the entire deconv period is ridiculously short, NaN it out!
    if size(x,1)<5
      d.mce0.data.fb(s:e,mcech)=NaN;
      continue
    end

    % Identify possible glitches of several types.
    %   1. Glitches of any kind, including transient spikes;
    %   2. Large steps or flux jumps
    %   3. Small, common-mode, polarized glitches
    % (see note /~spuder/analysis_logbook/analysis/20120522_polarjump/)
    mask=make_glitch_masks(x,xj,p,mcech);

    % Grow forward/backward in time around the "trigger points"
    % - for single sample glitch this will NaN 1s of time centered on
    % the glitch, unless d.tf exists.  in that case glitch around 
    % central sample of deconv kernel
    [mask,maskreg]=grow_glitch_masks(mask,d.tf(imce,:),samprate,s,e);

    % Keep track of large steps in this deconv block
    dg.njump(mcech,1)=dg.njump(mcech,1)+sum(mask{2},1)'/sum(maskreg);
    for ich=1:length(mcech)
      gl=find_blk([false;mask{2}(:,ich);false]);
      dg.fjlist(mcech(ich)).sf=[dg.fjlist(mcech(ich)).sf; gl.s+s-2];
      dg.fjlist(mcech(ich)).ef=[dg.fjlist(mcech(ich)).ef; gl.e+s-2];
    end

    % Make mask of channels to destep.  This will include
    %   1. Channels in maskjf
    %   2. Channels in same column
    %   3. Channels in same row and same readout card
    % But only if total # steps < 10.
    maskds=make_destep_mask(mask{2},sum(maskreg),p,mcech,dg.maxstep);

    % Also always destep small-polarized-glitches
    % Not doing for now -- one thing at a time!
    % maskds=maskds|mask{3};

    % Destep every glitch in its own channel AND its pol partner.
    maskds=maskds|make_partner_mask(mask{1},p,mcech);

    % Let window size be 1/2 second
    winsize=ceil((samprate+1)/2);

    % sample index, for convenience
    idx=1:size(maskds,1);

    for ich=1:length(mcech)

      % If the channel is hopeless, leave it alone.
      dg.njump(mcech(ich),2)=dg.njump(mcech(ich),2)+sum(maskds(:,ich))/sum(maskreg);
      if sum(maskds(:,ich))>=(dg.maxstep*sum(maskreg))
        disp(['Failing to deglitch channel ' num2str(ich) ', too many glitches.']);
        continue
      end

      % mask the glitches
      x(mask{1}(:,ich)|maskds(:,ich),ich)=NaN;

      % Make a list of the starts and ends of NaN regions around the
      % flux jumps.  Due to ringing of the filter, regions around the elnod
      % may have been masked out above
      gl=find_blk(maskds(:,ich));

      % Ignore NaN periods at start or end.
      cc=(gl.s>1) & (gl.e<idx(end)) & isfinite(gl.e);
      gl.s=gl.s(cc);
      gl.e=gl.e(cc);

      % Keep track of number of times destepped in this deconv block
      dg.njump(mcech(ich),3)=dg.njump(mcech(ich),3)+length(gl.s);

      if length(gl.s)>0
        disp(['Destepping channel ' num2str(ich) ', ' num2str(length(gl.s)) ' times']);
      end

      for j=1:length(gl.s)

        % skip if already merged into prev destep
        if(j>1 && gl.s(j)<=gl.e(j-1))
          gl.e(j)=max(gl.e(j-1),gl.e(j));
          continue
        end

        % now expand out to the first/last NaN in the fj region
        gl.s(j)=find(~isnan([0;x(1:gl.s(j),ich)]),1,'last');
        gl.e(j)=gl.e(j)-1+find(~isnan([x(gl.e(j):end,ich);0]),1,'first');

        % Remove step associated with large glitch
        x=destep(x,ich,idx,gl.s(j),gl.e(j),winsize);
      end
    end 

    % put the data back
    d.mce0.data.fb(s:e,mcech)=x;

  end % loop over rx

end % loop over deconv blocks

return

%%
function mask=make_glitch_masks(x,xj,p,ch)

  % First: find large excursions in "grass units"

  % take the abs of the diff to look for sharp spike in the tod
  y=abs(diff(x)); y=[y;zeros(size(y(1,:)))];

  % normalize by the median - this gives a measure of spike height
  % versus the "grass" which doesn't care about the presence of the
  % spikes
  y1=bsxfun(@rdivide,y,2*nanmedian(y));
  % take max over 2 samples to symmetrize zero-padding.
  % this avoids failure of the AND in mask{2}.
  y1=max(y1,[zeros(size(y(1,:))); y1(1:(end-1),:)]);

  % Second: find large deltas in feedback units
  y2=find_steps(x);

  % Third: find flux jumps recorded in numfj register
  y3=abs(diff(xj)); y3=[y3;zeros(size(y3(1,:)))];

  % Fourth: construct common-mode step parameter.
  ind=make_ind(p);
  cha=ismember(ch,ind.rgla); chb=ismember(ch,ind.rglb);
  y4=abs(nanmedian(diff(x(:,cha)-x(:,chb),1),2));
  y4=y4/nanmedian(y4)/1.5;
  y4=[y4;zeros(size(y4(1,:)))];

  % And construct output masks using these four parameters.
  mask{1}=(y1>8);
  mask{2}=(mask{1} & y2>=40) | y3>0;
  mask{3}=repmat(y4>5,1,length(ch));

  return

%%
function [mask,maskreg]=grow_glitch_masks(mask_in,tf,samprate,s,e)

  if isfield(tf,'deconvkern')
    % take only current receiver
    % determine which deconv kernel was applied
    test=[tf(:).dc];
    tfindx=find([mean([s e])>=[test.sf] & mean([s e])<=[test.ef]]);
    dckern=tf(tfindx).deconvkern;
    if isempty(dckern)
      maskreg=ones(samprate+1,1);
    else
      [a,maxind]=max(dckern);
      kernlen=length(dckern);
      masklen=kernlen*2;
      % make mask odd and place central sample at center
      if ~mod(masklen,2); 
        masklen=masklen+1;
      end
      maskreg=zeros(masklen,1,'single');
      maskreg(1:kernlen)=1;
      maskreg=circshift(maskreg,ceil(masklen/2)-maxind);
    end
  else
    maskreg=ones(samprate+1,1,'single');
  end

  if iscell(mask_in)
    for i=1:length(mask_in)
      mask{i} = false(size(mask_in{i}));
      % conv only operates on numbers, not logicals
      tmpmask = single(mask_in{i});
      for jj=1:size(mask_in{i},2)
        mask{i}(:,jj) = logical(conv(tmpmask(:,jj), maskreg, 'same'));
      end
    end
  else
    mask = false(size(mask_in));
    tmpmask = single(mask_in);
    for jj=1:size(mask_in,2)
      mask(:,jj) = logical(conv(tmpmask(:,jj), maskreg, 'same'));
    end
  end

  return

%%
function x=destep(xin,ich,idx,sg,eg,winsize)

  x=xin;

  % Construct a window on each side of glitch
  cpre=(idx<sg & sg-idx<=winsize);
  cpost=(idx>eg & idx-eg<=winsize);

  % Simple DC offset
  delta=nanmean(x(cpost,ich))-nanmean(x(cpre,ich));

  % Adjust following part of trace for delta
  if isfinite(delta)
    x(eg:end,ich)=x(eg:end,ich)-delta;
  end

  return

%%
function maskds=make_destep_mask(maskfj,kernlen,p,ch,maxsteps)

  % Identify channels that have stepped, but not too many times.
  dochans=sum(maskfj,1)>0 & sum(maskfj,1)<=maxsteps*kernlen;

  maskds=false(size(maskfj));
  for ich=1:length(ch)
    if ~dochans(ich)
      continue
    end

    % Determine the same col and rows for the mce.  Note that same-row
    % also requires same-readout-card, since we don't get crosstalk 
    % across RCs.
    ich_samerow=find(p.mce==p.mce(ch(ich)) & p.mce_row==p.mce_row(ch(ich)) ...
      & floor(p.mce_col/8)==floor(p.mce_col(ch(ich))/8));
    ich_samecol=find(p.rx==p.rx(ch(ich)) & p.mce_col==p.mce_col(ch(ich)));

    % only correct channels included in ch set passed to deglitch
    % and find indices of chinduse within ch
    ch_ind=find(ismember(ch,[ich_samecol(:)' ich_samerow(:)']));

    % And add to mask
    maskds(:,ch_ind)=maskds(:,ch_ind) | repmat(maskfj(:,ich),1,length(ch_ind));
  end

  return

%%
function mask=make_partner_mask(maskin,p,ch)

  mask=maskin;

  for ich=1:length(ch)
    ch_partner=find(p.rx==p.rx(ch(ich)) & p.tile==p.tile(ch(ich)) ...
                  & p.det_row==p.det_row(ch(ich)) & p.det_col==p.det_col(ch(ich)));
    ch_partner=ch_partner(ch_partner~=ch(ich));
    if length(ch_partner)==1
      ich_partner=find(ch==ch_partner);
      mask(:,ich)=mask(:,ich) | mask(:,ich_partner);
    end
  end

  return

