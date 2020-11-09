function sqws = shiftchopref(sqw,shiftind)
%
% function sqws = shiftchopref(sqw,shiftind)
%
% TSG 2018-03-18
%
% Using circshift with a chop reference generally doesn't work since 
% it doesn't preserve the phase.  This function applies a shift
% while matching the phase at the beginning/end of the chop ref.
%
% INPUTS:
%   sqw = vector of chop ref, boolean
%   shiftind = # indices shift to apply, positive or negative integer 
%
% OUTPUT:
%   sqws = shifted chop ref, vector same length as sqw, boolean
%

% Have output shape match input shape.  Need to avoid CAT errors.
if size(sqw,1)==1
  sqw = sqw';
  flipflag = 1;
else
  flipflag = 0;
end

% Positive shift = move entries from the end to the beginning
% Negative shift = move entries from the beginning to the end
% Loop until we're at a point that matches the phase of the first 
% sample, copy that chunk to the beginning (or end if negative),
% then trim off appropriate number of points from the end.
if shiftind>=0
  for ii = 2:length(sqw)
    if sqw(ii)~=sqw(ii-1)
      break
    end
  end
  for jj = length(sqw)-1:-1:1
    if sqw(jj+1) - sqw(jj) == sqw(ii) - sqw(ii-1)
      sqw_wrap = sqw((jj-ii-shiftind+2):(jj-ii+1));
      break
    end
  end
  sqws = [sqw_wrap; sqw(1:end-shiftind)];
else 
  for ii = length(sqw)-1:-1:1
    if sqw(ii+1) ~= sqw(ii)
      break
    end
  end
  for jj = 2:length(sqw)
    if sqw(jj) - sqw(jj-1) == sqw(ii+1) - sqw(ii)
      sqw_wrap = sqw(jj+1:jj-shiftind);
      break
    end
  end
  sqws = [sqw(1-shiftind:end); sqw_wrap];
end

if flipflag
  sqws = sqws';
end

return