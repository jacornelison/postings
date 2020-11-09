% d=find_steps(x,n)
%
% Step-finding parameter to distinguish steps from
% transients.  Input data is x(samples,chans),
% and output is d(samples,chans), in the same units.
%
% Specify the number of samples to consider on each
% side of the step: n=3 by default.
function d=find_steps(x,n)

  if(nargin<2 || isempty(n))
    n=3;
  end
  if(nargin<3 || isempty(ww))
    ww=1;
  end

  x1=x(1,:);
  xend=x(end,:);

  % Rising-edge step
  y1=Inf;
  y2=-Inf;
  for i=1:n
    y1=min(y1,[x((i+1):end,:);repmat(xend,i,1)]);
    y2=max(y2,[repmat(x1,i,1);x(1:(end-i),:)]);
  end
  d=y1-y2;

  % Falling-edge step
  y1=-Inf;
  y2=Inf;
  for i=1:n
    y1=max(y1,[x((i+1):end,:);repmat(xend,i,1)]);
    y2=min(y2,[repmat(x1,i,1);x(1:(end-i),:)]);
  end
  d=max(y2-y1,d);

  d(~isfinite(d))=NaN;

  return
