function [bincenter,n] = hfill(vals,nbin,low,high,weights,opt)
% [bincenter,n] = hfill(vals,nbin,low,high,weight,opt)
%
% Histogram data into equally sized bins specified by number
% and range.
%
% vals are the events to be binned
% nbin defaults to 10
% low edge of the bins defaults to min(vals)
% high edge of the bins defaults to max(vals)
% weight of the events default to 1
% opt to be passed through to hplot when plotting
% 
% If no output specified plots the resulting histogram.
%
% eg: hfill(randn(100),100,-3,3,rand(100));

if(~exist('nbin','var'))
  nbin=[];
end
if(~exist('low','var'))
  low=[];
end
if(~exist('high','var'))
  high=[];
end
if(~exist('weights','var'))
  weights=[];
end
if(~exist('opt','var'))
  opt=' ';
end

% Set weights if not provided
if(isempty(weights))
  weights=ones(size(vals));
end

if((~isreal(vals))||(~isreal(weights)))
  error('Data (and weights) must be real');
end

% Max and min require vector data
vals=vals(:)'; weights=weights(:)';

if(isempty(nbin))
  nbin=10;
end
if(isempty(low))
  low=min(vals);
end
if(isempty(high))
  high=max(vals);
  % Make max just fall in top bin
  high=high+(high-low)*1e-9;
end

% stop it from crashing if all values same
if(low==high)
  low=low*0.9; high=high*1.1;
end
if(low==0&high==0)
  low=-1; high=+1;
end

% If there is non-finite data or weights remove
if(any(~isfinite(vals))|any(~isfinite(weights)))
  ind=isfinite(vals)&isfinite(weights);
  vals=vals(ind);
  weights=weights(ind);
  %warning('Non-finite input data (removed)');
end

% C code for speed
[bincenter,n]=hfillc(double(vals),double(nbin),double(low),double(high),double(weights));

if(nargout==0)
  hplot(bincenter,n,opt);
end
