function idx = inrange(A,low,high,inc)
% Find the values or array A that are in the range of low and high
% inc -- inclusive flag.
%   true (default) -- inclusive low & high
%   false -- exclusive low & high
%   [bool bool] -- choose inc/exc for [low high] individually

if ~exist('inc','var')
    inc = true;
end

incdiff = diff(inc);

% if we designate low&hi, but they're both the same pretend we didn't.
if length(inc)==2 && incdiff~=0

    % If we're here and low isn't true, then high must be.
    if inc(1)
        idx = A>=low & A<high;
    else
        idx = A>low & A<=high;
    end
else

    % Use any() in case we're here because we designated both low & hi.
    if any(inc)
        idx = A>=low & A<=high;
    else
        idx = A>low & A<high;
    end
end

idx = logical(idx);
