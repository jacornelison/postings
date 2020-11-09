function S=identify_skips(d,fs)
% S=identify_skips(d,fs)
%
% Identifies skips in TOD data due to errors in the data acquisition process.
% Currently, this identifies 5 classes of errors:
%
%   1. Skipped samples ('skipped samples')
%   2. Filter failures ('filter failures')
%   3. Dropped frames in the mediator ('mediator dropped frames')
%   4. Dropped samples in an MCE ('drop in mce %d')
%   5. Bad synchronization in an MCE ('bad sync in mce %d')
%
% This was originally a subfunction of plot_skips() named find_skips().
%
% INPUTS
%   d     Standard TOD structure
%
%   fs    Standard block start/stop structure containing the sf and ef
%         fields. Length of sf and ef must be 1.
%
% OUTPUTS
%
%   S     A structure containing information that indicates where skips
%         occurred within the data stream. Fields are:
%
%           .name    Name of problem, matching cases above in parentheses
%           .sf      Start index of skip, from start of TOD
%           .ef      End index of skip, from start of TOD
%           .desc    Human-readable description of skip event
%
% EXAMPLE
%
%   load('data/real/201603/20160315F10_dk068_tod.mat');
%   skips = identify_skips(d, fsb);
%

  S = struct('name',{{}}, 'sf',{[]}, 'ef',{[]}, 'desc',{{}});

  % Calculate the frame length from the data
  framelen = size(d.mce0.syncBox.sampleNumber,1) / ...
      size(d.mce0.syncBox.status,1);
  error_on_frames = false;
  % Make sure the frame length is an integer. If not, warn about the situation
  % and be prepared to bail out later if we have to.
  if framelen-floor(framelen) > framelen*sqrt(eps())
    warning('Non-integer frame length encountered.')
    error_on_frames = true;
    framelen = round(framelen);
  end


  % Case 1: Skipped samples
  mask = d.antenna0.syncBox.sampleNumber(fs.sf:fs.ef)==0 ...
      & d.t(fs.sf:fs.ef)~=0;
  S = add_skip(S, mask, 'skipped samples');

  % Case 2: Filter failures
  mask = d.antenna0.syncBox.sampleNumber(fs.sf:fs.ef)==0 ...
      & all(d.mce0.syncBox.sampleNumber(fs.sf:fs.ef,:)==0, 2);
  S = add_skip(S, mask, 'filter failures');
  filtmask = mask; % needed for MCE droped samples later

  % Case 3: Dropped frames in the mediator
  fidx = double(d.array.filter.idx(fs.sf:fs.ef));
  mask = [0;diff(fidx)]>1 | ([0;diff(fidx)]<0 & fidx>0);
  S = add_skip(S, mask, 'mediator dropped frames');

  nmce = size(d.mce0.syncBox.sampleNumber,2);

  % Case 4: Dropped samples in an MCE
  for ii=1:nmce
    mask = d.mce0.syncBox.sampleNumber(fs.sf:fs.ef,ii)==0 & ~filtmask;
    S = add_skip(S, mask, sprintf('drop in mce %d', ii-1));
  end

  % Case 5: Bad synchronization in an MCE
  for ii=1:nmce
    stat = repmat(d.mce0.syncBox.sync_status(:,ii), 1, framelen);
    stat = reshape(stat', 1, [])';
    mask = stat(fs.sf:fs.ef) > 4;

    % Bail out and make the user deal with a situation where the calculated
    % frame length was initially non-integer and we've also identified a
    % problem using per-frame data. This is because the replication we did
    % in stat may be wrong.
    if error_on_frames && any(mask)
      error('A potential skip was identified, but the frame length is unreliable.')
    end

    S = add_skip(S, mask, sprintf('bad sync in mce %d', ii-1));
  end

  % All the indices added by add_skip are relative to fs.sf, so fix that up.
  S.sf = S.sf + fs.sf - 1;
  S.ef = S.ef + fs.sf - 1;

  % With the indices all fixed up, now add descriptions that include correct
  % offsets.
  S.desc = arrayfunc(@(sf,ef) fmtmsg(sf,ef,framelen), S.sf, S.ef);
end

function S=add_skip(S,mask,txt)
  [sf,ef] = cut2idx(mask);
  if length(sf) == 0
    return
  end

  ii = length(S.sf) + 1;
  jj = ii + length(sf) - 1;

  S.sf(ii:jj,1) = sf;
  S.ef(ii:jj,1) = ef;
  [S.name{ii:jj,1}] = deal(txt);
end

function [sf,ef]=cut2idx(mask)
% [sf,ef]=cut2idx(mask)
%
% Turns a mask of values into start and end index pairs.
%

  mask = cvec(mask);
  sf = find([mask(1); diff(mask)>0]);
  ef = find([diff(mask)<0; mask(end)]);
end

function desc=fmtmsg(sf,ef,framelen)
% desc=fmtmsg(sf,ef)
%
% Pretty-prints a message describing the number of samples/frames which
% are identified as bad.
%

  if mod(sf,framelen)==1 && mod(ef,framelen)==0
    if ef-sf+1 == framelen
      desc = sprintf('Single frame at %d', sf);
    else
      desc = sprintf('Run of %d frames at %d', (ef-sf+1)/framelen, sf);
    end
  else
    if ef==sf
      desc = sprintf('Single sample at %d', sf);
    else
      desc = sprintf('Run of %d samples at %d', (ef-sf+1), sf);
    end
  end
end

