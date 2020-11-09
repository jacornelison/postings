function chop = chop_deglitch(chop)
% chop = chop_deglitch(chop)
%
% Gets rid of single sample transitions in the chop reference
% (i.e. ..., high, high, low, high, high, ...) by extending them to two 
% samples. Note that this approach assumes that the transition is real, 
% but we are chopping a bit too fast for the sample rate. If the chop signal
% is legitimately noisy with false transitions, then you would want to do 
% something else.
%
% [Arguments]
%   chop  Input chop reference to be deglitched. This timestream should 
%         already be digitized, i.e. consists of only ones and zeros.
%
% [Returns]
%   chop  Output chop reference has the same length as input, just modified 
%         to eliminate single sample transitions.

% Last update: 2014-02-06 CAB

% Find all spots where we change by +-1 between samples i-1 and i, then 
% change by -+1 between samples i and i+1.
glitchfinder = diff(chop) - circshift(diff(chop), 1);

% Single 0 samples will show up as value of +2.
% This extra if statement is necessary for some reason.
if ~isempty(find(glitchfinder == 2))
  for i=find(glitchfinder == 2)
    % Flip a coin to decide whether to extend the 0 samples in the left or 
    % right direction.
    if (rand(1) > 0.5) || (i(1) == 1)
      chop(i+1) = 0;
    else
      chop(i-1) = 0;
    end
  end
end

% Single 1 samples will show up as value of -2.
% This extra if statement is necessary for some reason.
if ~isempty(find(glitchfinder == -2))
  for i=find(glitchfinder == -2)
    % Flip a coin to decide whether to extend the 1 samples in the left or 
    % right direction.
    if (rand(1) > 0.5) || (i(1) == 1)
      chop(i+1) = 1;
    else
      chop(i-1) = 1;
    end
  end
end
