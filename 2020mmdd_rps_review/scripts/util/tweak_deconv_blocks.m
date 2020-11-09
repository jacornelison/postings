function b=tweak_deconv_blocks(d)
% b=tweak_deconv_blocks(d)
%
% tweak the deconv data blocks to retain as much of the
% the elnods as possible.  find data regions of constant
% samprate.
%
% begin/end regions are NaNd in deconv_scans, and the
% timing on the feature bits is not exact (hence the need
% for tweak_elnods).  this ends up cutting off portions
% of the elnods.  expand the deconv blocks out to a change
% in samprate to preserve as much data as possible.  also,
% try to combine en & fs blocks if they are contiguous.
% the feature bit drops to 0 between them, so the blocks
% are found separately.
%
% Note: expanding the deconv blocks may end up including
% large steps in the detector data before the elnods.  this is ok
% for deconv since the deconv kernels are FIR, and that data will
% not be selected in the en/fs blocks.  but, this is bad for deglitch
% which will may try to level those steps.  so, use separate data
% blocks for deconv than for deskip/delitch.
%
% this is valid since it's ok to deconv across a period of
% constant samprate.

% form the data block for deskip/deglitch/deconv
% most en are contiguous with fs but in some b2 first season data they are not
% identify en that are adjacent to fs and combine if possible

sampratio=length(d.t)/length(d.ts);

% begin by finding the regions where feature bit is set
% but exclude partial load curves (bit=2^14)
% this should create at most 3 distinct data blocks
bb=d.array.frame.features>0 & d.array.frame.features<2^14;

% the feature bit drops to zero between fs and en blocks for a few samples
% assume that drops to 0 for fewer than 10 samples are these regions
% find them and remove the drops to 0 from bb
trans=find(abs(diff(bb))==1);
indx=find(diff(trans)<=10);
for ii=1:length(indx)
  bb(trans(indx(ii))+1:trans(indx(ii)+1))=1;
end

% create the standard block
b=find_blk(bb,sampratio);

% expand these regions out to where the samprate changes
% don't expand out to where mce data isn't being acquired
% i.e. where all d.mce0.header.clock_counter=0 because
% of init_TES or combine off transition.
nsnap=d.array.frame.nsnap;
nf=reshape(repmat(nsnap',sampratio,1),[],1);
clock_counter=d.mce0.header.clock_counter;
has_clock=any(clock_counter>0,2);
ratechange=find(0~=diff((nf>1) & has_clock));

for ii=1:length(b.sf)
  % step ahead one frame to make sure we're not already past the rate change
  b.sf(ii)=b.sf(ii)+sampratio;

  % expand backward to the samprate change
  % if b.sf is smaller than all the ratechange just
  % expand to the first sample
  if nansum(b.sf(ii)<ratechange)==length(ratechange);
    b.sf(ii)=1;
  else
    % find the largest samprate change that it's
    % greater than and expand to down to it
    indx=find(b.sf(ii)>=ratechange);
    % and make sure you're past a patch where the mce data
    % is exactly 0 during init_TES
    datzero=find(nansum(d.mce0.data.fb,2)==0);
    datindx=find(b.sf(ii)>=datzero);
    if ~isempty(datindx)
      b.sf(ii)=max([ratechange(indx(end))+2 datzero(datindx(end))+1]);
    else
      b.sf(ii)=ratechange(indx(end))+2;
    end
  end
  % step back one frame for safety
  b.ef(ii)=b.ef(ii)-sampratio;

  % expand forward to the samprate change
  % find the smallest samprate change that it's
  % smaller than and expand up to it
  % If there is no ending rate change, the final elnod might have been cut off
  % and we'll leave b.e the same.
  indx=find(b.ef(ii)<=ratechange);
  if ~isempty(indx)
    % want end index to be right before ratechange so that the block defined
    % by b.sf and b.ef has constant rate.
    b.ef(ii)=ratechange(indx(1))-1;
  end

  % now adjust the slow indices
  b.s(ii)=floor((b.sf(ii)-1)/sampratio)+1;
  b.e(ii)=floor((b.ef(ii)-1)/sampratio)+1;
end
