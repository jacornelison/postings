function [psd,f]=mean_psd(d,fs)
% [psd,f]=mean_psd(d,fs)
%
% Take the mean psd's of the half scans.
% one can regard the output as one-sided PSD in units of ADU^2/Hz
% or two-sided PSD in units of ADU^2*seconds
%
% For many years it said above that the output units of this function
% are ADU^2/Hz - That is incorrect if one assumes the output PSD
% is double-sided - which is perhaps the default assumption.
%
% Note that in the Fourier Transform delta f = 1/(2*delta t))
% 
% The practical outcome is this:
% - If you want a two-sided PSD in units of ADU^2/Hz multiply the output
% of this function by 2.
% - If you want a two-sided PSD in units of ADU^2*seconds do nothing.
%
% (When working with pair sum/diff data note that there is potentially
% an extra layer - if you want the per single detector equivalent
% in ADU^2*seconds then multiply the output of this function by 2.)
%
% NB: input must have at least mean (or median) removed on half scan
% basis or else window causes mixing of DC into other freq
%
% e.g.
% load data/20110505C10_dk068_tod
% [p,ind]=get_array_info('20110501');
% d=sumdiff_pairs(d,p,fs,ind.a,ind.b);
% d=filter_scans(d,fs,'p0');
% [psd,f]=mean_psd(d,fs,samprate);
% loglog(f,psd(:,ind.rglb)); axis tight; ylim([1e-4,1e4]);

disp('mean_psd...');

[sampratio,samprate]=get_sampratio_rate(d);

% all blocks must be equal length - code will choke if not
n=fs.ef(1)-fs.sf(1)+1;

% use a mild window  
w=tukeywin(n,0.2);
  
% loop over half-scans
% (go backwards to cause alloc of large averaging array at full size
% on first iteration)
for i=length(fs.sf):-1:1
  s=fs.sf(i); e=fs.ef(i);
  v=d.mce0.data.fb(s:e,:);
  
  % apply window
  v=v.*repmat(w,[1,size(v,2)]);
  
  % pad to get more efficient fft
  n=2^nextpow2(n);
  
  % take fft
  ft=fft(v,n);
  
  % undo normalization change due to window application (and padding)
  ft=ft*sqrt(n/sum(w.^2));
  
  % throw away symmetric (neg freq) part
  nup=ceil((n+1)/2);
  ft=ft(1:nup,:);

  % scale so result not a function of input length
  ft=ft/sqrt(n);
  
  % generate auto spectra
  hs=ft.*conj(ft);
  
  % force result real instead of fold-over-and-add which does
  % the same thing
  hs=real(hs);

  % store for average
  hsa(:,:,i)=hs;
  
end

% take average over half-scans
psd=nanmean(hsa,3);

% convert from ADU^2/samprate to ADU^2/Hz
psd=psd/samprate;

% calc the half freq axis
switch rem(n,2)
  case 0
    f=[0:n/2]'*samprate/n;
  case 1
    f=[0:n/2-0.5]'*samprate/n;
end

return
