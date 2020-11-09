% SQUARE_DEMOD  Demodulate a square-wave signal
%   from a chopped source
%
% [C, S] = SQUARE_DEMOD (X, B)
%
%    X = modulated detector signal
%    B = chop signal (1 for source on, <=0 for source off)
%    C = demodulated cosine part
%    S = demodulated sine part
%
% C and S are returned at one value per chop cycle,
% effectively downsampled to the chope rate.  (See
% below for an option to expand the outputs to the 
% same length as the input.)
%
% [C, S, IC, IS] = SQUARE_DEMOD(...) also returns IC
% and IS, the index at the middle of each chop cycle
% for cosine and sine demodulation.
%
% X may be a matrix of N rows by M columns, representing
% N samples in each of M channels.  In this case, the
% length of B must be N.  C and S will each have M
% columns.
%
% SQUARE_DEMOD (..., 'deglitch') removes glitches from
% the chop reference signal.  It is assumed that the
% true chop is a clean square wave with a fixed period.
%
% [C, S] = SQUARE_DEMOD (..., 'expand') expands C and S
% to have the same number of samples as the input X.  This
% is done by repeating the demodulated values over a chop
% cycle.  It's also possible to interpolate using a spe-
% cified method, e.g. 'linear' or 'spline'.
%
% [C, S] = SQUARE_DEMOD (..., 'keep_single') keeps parts
% of the chop reference where a single sample is high or
% low.  Otherwise, these are interpreted as bad reference
% values and result in a NaN in the demodulated output.
%
% UPDATE 20180601
% [C, S] = SQUARE_DEMOD (..., 'highres')
% Chop signal B is treated as high-resolution chop reference, 
% i.e. from cleanchopref.m.  B will be downsampled 
% to resolution of X but will use locations of transitions in 
% high-res chop ref when evaluating time-domain kernel.
% This replaces the use of rand, making SQUARE_DEMOD deterministic.
% ONLY USE THIS OPTION IF YOU ARE SUPPLYING HIGH RES CHOP REF.
% Example:
% >> b = d.antenna0.pmac.fast_aux_input(:,4);
% >> b = b > mean(b);
% >> ch = cleanchopref(b);
% >> [c,s] = square_demod(x,ch.refclean_hires,'highres');

function [c, s, ic, is] = square_demod (x, b, varargin)

% Parse extra options
DO_DEGLITCH = 0;
DO_INTERP = '';
NAN_SINGLE = true;
HIGH_RES = 0;
if (nargin > 2)
  for (ii = 1:length(varargin))
    switch (lower (varargin{ii}))
      case 'deglitch', DO_DEGLITCH=1;
      case 'expand', DO_INTERP='nearest';
      case {'linear','spline','nearest'}, DO_INTERP = lower(varargin{ii});
      case 'keep_single', NAN_SINGLE = false;
      case 'highres', HIGH_RES = 1;
      otherwise, error (['Unknown option ' varargin{ii}]);
    end;
  end;
end;

if max(b)>1 || min(b)<0 || any(~isfinite(b))
  error(['Chop reference b should be in range 0-1.']);
end

b = (b(:) >= 0.5);
if size(x,1)==1
  x = x(:);
end;

% If input chop ref is high-resolution, downsample it and find transitions
if (HIGH_RES)
  b_old = b;
  i_hires = linspace(1,size(x,1),length(b))';
  dif_hires = diff([b(1); b]);
  ind_du = i_hires(dif_hires>0);
  ind_ud = i_hires(dif_hires<0);
  b = interp1(i_hires,b,1:size(x,1));
  b = (b(:) > 0.5);
end

% Deglitch the chop reference, if requested
if (DO_DEGLITCH)
  b_old = b;
  b = (b - 0.5) *2;
  [p f] = time2psd (b, 1);
  p (1:4) = 0;
  p ((end-3):end) = 0;
  [pmax imax] = max (p);
  fchop = f (imax);
  chcos = cos (2*pi*fchop * (1:length(b))');
  chsin = sin (2*pi*fchop * (1:length(b))');
  [bflt aflt] = butter (2, 2*fchop/5, 'low');
  bcos = filtfilt (bflt, aflt, chcos .* b);
  bsin = filtfilt (bflt, aflt, chsin .* b);
  ph = atan2 (bsin, bcos);
  b = cos (2*pi*fchop * (1:length(b))' - ph) > 0;
  % b90 = sin (2*pi*fchop * (1:length(b))' - ph) > 0;
end;

% Note: chop cycle phases are as follows
%
% phase 1 low
% phase 2 high
% phase 3 high
% phase 4 low
% phase 1 low
% phase 2 high
% ...
%
% cosine component is -1 + 2 + 3 - 4
% sine component is   -2 + 3 + 4 - 1
%
% if first chop sample is high, it's 1->2
% if first chop sample is low, it's 4->1
% if last chop sample is high, last+1 is 2->3 or 3->4 (but useless either way)
% if last chop sample is low, last+1 is either 1->2 or 4->1 (finish with sine or cos?)

% Make vector of on/off indices
j{2} = find([b;0] & [1;~b]);    % 1->2, low->high
j{4} = find([~b;0] & [0;b]);    % 3->4, high->low

% handle case with last sample low and choose 1->2
if ~b(end)
  nlow = length(b) - j{4}(end) + 1;
  nave = mean(diff(j{2}));
  if nlow >= nave/4 && nlow >= 2
    j{2} = [j{2}; length(b)+1];
    % Make same adjustment for high-res option
    if (HIGH_RES); ind_du = [ind_du; length(b)+1]; end
  end
end

if length(j{2})<2 || length(j{4})<2
  error(['Need at least two chop cycles.']);
end

% Make vector of 90 degree offset indices
% Use knowledge of high-res transitions if given
if (HIGH_RES)
  % Make sure we stick with convention that it's 1->2 if first sample high
  if b(1); ind_du = [1; ind_du]; end
  tmp1 = [j{2}(j{2}<j{4}(end)),j{4}(j{4}>j{2}(1))];
  tmp2 = [ind_du(ind_du<ind_ud(end)),ind_ud(ind_ud>ind_du(1))];
  j{3} = mean(tmp1,2);             % 2->3, high->high
  j{3} = round(j{3} + 1e-5*(mean(tmp2,2)-j{3}));
  tmp1 = [j{4}(j{4}<j{2}(end)),j{2}(j{2}>j{4}(1))];
  tmp2 = [ind_ud(ind_ud<ind_du(end)),ind_du(ind_du>ind_ud(1))];
  j{1} = mean(tmp1,2);             % 4->1, low->low
  j{1} = round(j{1} + 1e-5*(mean(tmp2,2)-j{1}));
else
  tmp = [j{2}(j{2}<j{4}(end)),j{4}(j{4}>j{2}(1))];
  j{3} = mean(tmp,2);             % 2->3, high->high
  j{3} = round(j{3} + 1e-5*(rand(size(j{3}))-0.5));
  tmp = [j{4}(j{4}<j{2}(end)),j{2}(j{2}>j{4}(1))];
  j{1} = mean(tmp,2);             % 4->1, low->low
  j{1} = round(j{1} + 1e-5*(rand(size(j{1}))-0.5));
end
% handle case with first sample low
if ~b(1)
  j{1} = [1;j{1}];
end

% handle case with last sample low and choose 4->1
if ~b(end) && j{2}(end)~=length(b)+1
  j{1} = [j{1}; length(b)+1];
end

if length(j{1})<2 || length(j{3})<2
  error(['Need at least two chop cycles.']);
end

% Construct cosine phase demod signal
[c,ic] = demod1(j,[1 2 3 4],x,NAN_SINGLE);

% Construct sine phase demod signal
[s,is] = demod1(j,[2 3 4 1],x,NAN_SINGLE);

% Interpolate to expand output to full length, if requested
if ~isempty(DO_INTERP)
  ctmp = zeros (size(x));
  stmp = zeros (size(x));
  jj = 1:size(x,1);
  wst = warning('query','MATLAB:interp1:NaNinY');
  warning('off','MATLAB:interp1:NaNinY');
  for i=1:size(x,2)
    ctmp(:,i) = interp1 (ic, c(:,i), jj + 1e-5*(rand(size(jj))-0.5), DO_INTERP, 'extrap');
    stmp(:,i) = interp1 (is, s(:,i), jj + 1e-5*(rand(size(jj))-0.5), DO_INTERP, 'extrap');
  end
  warning(wst);
  c = ctmp;
  s = stmp;
end

return

% Demodulate a sine or cosine part
% j = j{1,2,3,4} is list of indices where chop phase changes
% ph = order of phases to demodulate
% x = undemodulated data
function [z,iz] = demod1(j,ph,x,nan_singletons)
  jj = zeros(length(j{ph(1)})-1,5);
  for i=1:length(ph)
    tmpj = j{ph(i)}(j{ph(i)}>=j{ph(1)}(1) & j{ph(i)}<=j{ph(1)}(end));
    jj(:,i) = tmpj(1:size(jj,1));
  end
  jj(:,5) = j{ph(1)}(j{ph(1)}>j{ph(1)}(1));
  z = zeros(size(jj,1),size(x,2));
  sgn = [-1,1,1,-1];
  mul = zeros(size(x,1),1);
  for i=1:4
    js = jj(:,i);
    je = jj(:,i+1);
    n = zeros(size(x,1)+1,1);
    n(js) = 1./max(1,je-js);
    n(je) = n(je)-n(js);
    mul = mul + sgn(i) * cumsum(n(1:(end-1)));
  end
  hasnans = false;
  for k=1:size(x,2)
    tmpx = x(:,k);
    % Handle NaNs -- want to end up with NaN only in affected cycle;
    % don't let cumsum push the NaNs through whole time stream.
    % But also don't check isfinite on entire input x; for some reason,
    % isfinite is super expensive / slow.
    if hasnans
      cnan = ~isfinite(tmpx);
      if any(cnan)
        tmpx(cnan) = 0;
      end
    end
    tmpx = cumsum(tmpx .* mul);
    z(:,k) = diff([0;tmpx(je-1)]);
    if ~hasnans && ~isfinite(sum(z(:,k)))
      hasnans = true;
      tmpx = x(:,k);
      cnan = ~isfinite(tmpx);
      tmpx(cnan) = 0;
      tmpx = cumsum(tmpx .* mul);
      z(:,k) = diff([0;tmpx(je-1)]);
    end
    % NaN out affected cycles
    if hasnans && any(cnan)
      cnan = cumsum(cnan);
      z(diff([0;cnan(je-1)])>0,k) = NaN;
    end
  end

  % Demodulation fails whenever there are singleton chop samples,
  % i.e. high or low for only one sample.  Can handle this in either of
  % two ways -- either fudge around by counting the sample twice, or
  % else NaN out this chop cycle.
  ising = find(any(diff(jj,[],2)==1,2));
  % if ~isempty(ising)
  %   disp('Found some singletons');
  % end
  if nan_singletons
    for k=1:length(ising)
      z(ising(k),:) = NaN;
    end
  elseif ~isempty(ising)
    % Slower method that lets you safely overlap parts of chop
    % phases
    z(ising,:) = 0;
    for i=1:4
      js = jj(ising,i);
      je = jj(ising,i+1);
      % Fudge for single sample high or low
      if i==1 || i==3
        js(je==js) = js(je==js)-1;
      else
        je(je==js) = je(je==js)+1;
      end
      if any(js<1) || any(je>size(x,1))
        error(['Calculation ran out of bounds in singleton fix.']);
      end
      for k=1:length(ising)
        z(ising(k),:) = z(ising(k),:) ...
          + sgn(i) * mean(x(js(k):(je(k)-1),:),1);
      end
    end
  end
  
  z = z/2;
  iz=jj(:,3);

  return


