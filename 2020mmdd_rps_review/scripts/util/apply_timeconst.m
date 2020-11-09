function d=apply_timeconst(d,tau,n,fsamp)
% d=apply_timeconst(d,tau,n,fsamp)
%
% d     = TOD
% tau   = time constants in seconds, not ms
%         can be single value or nchan array
% n     = size of FIR filter in samples (default 100)
% fsamp = sample rate in samples/s (default is output of get_samprate_ratio)

if ~exist('fsamp') || isempty(fsamp)
  fsamp=get_sampratio_rate(d);
end
% Assume sample rate is in Hz.  Should use
% "samprate" output from get_sampratio_rate.
if ~isscalar(fsamp)
  error(['fsamp should be a scalar sample rate, as returned by get_sampratio_rate.']);
end

% Default to 100 samples
if ~exist('n','var') || isempty(n)
  n=100;
end

% Make sure we have the correct number and
% shape of time constants
if isscalar(tau)
  tau=tau*ones(1,size(d.mce0.data.fb,2));
end
if numel(tau)~=size(d.mce0.data.fb,2)
  error(['Number of time constants should match number of channels.']);
end
tau=reshape(tau,1,[]);

% Loop over unique values, generate FIR filter,
% apply for those that aren't zero.
tau(~isfinite(tau) | tau<=0) = 0;
[tau_u i j]=unique(tau);
for k=1:length(tau_u)
  if tau_u(k)==0
    continue
  end
  disp(['tau=' num2str(tau_u(k)) ', ' num2str(sum(j==k)) ' channels.']);
  ck=(j==k);
  b=get_timeconst_fir(fsamp,tau_u(k),n);
  % y(:,ck)=conv(y(:,ck),b,'same');
  tmpy=filter2(b,d.mce0.data.fb(:,ck),'full');
  d.mce0.data.fb(:,ck)=tmpy(1:size(d.mce0.data.fb,1),:);
end

return

%%%%%%%%%%%%%%%%%%%%%%%
function b=get_timeconst_fir(fsamp,tau,n)

% tt must be a column vector for kabbalistic reasons
tt=(0:(n-1))'/fsamp;
b=exp(-tt/tau);
b=b/sum(b);
% filter2 convolves in time reversed sense...  use this
% reversal for filter2, not for conv.
b=b(end:-1:1);

return
