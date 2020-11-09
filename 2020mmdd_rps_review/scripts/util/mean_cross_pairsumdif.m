function [sSS,sSD,sDS,sDD,f]=mean_cross_pairsumdif(d,fs,p,ind)
% [sSS,sSD,sDS,sDD,f]=mean_cross_pairsumdif(d,fs,p,ind)
% sSS: spectra Sum x Sum, sSD: spectra Sum x Dif
% sDS: spectra Dif x Sum, sDD: spectra Dif x Dif
%
% Take the mean cross spectra of pair dif/sum half scans.
% Units are ADU^2/Hz
%
% NB: input must have at least mean (or median) removed on half scan
% basis or else window causes mixing of DC into other freq
%
% e.g.
% load data/20110505C10_dk068_tod
% [p,ind]=get_array_info('20110501');
% d=sumdiff_pairs(d,p,fs,ind.a,ind.b);
% d=filter_scans(d,fs,'p0');
% [sSS,sSD,sDS,sDD,f]=mean_cross_pairsumdif(d,fs,p,ind);

disp('mean_cross_pairsumdif...');

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
  % pair dif timestream
  v=d.mce0.data.fb(s:e,ind.rglb);
  % pair sum timestream
  u=d.mce0.data.fb(s:e,ind.rgla);
  
  % apply window
  v=v.*repmat(w,[1,size(v,2)]);
  u=u.*repmat(w,[1,size(u,2)]);
  
  % pad to get more efficient fft
  n=2^nextpow2(n);
  
  % take fft
  ftv=fft(v,n);
  ftu=fft(u,n);
  
  % undo normalization change due to window application (and padding)
  ftv=ftv*sqrt(n/sum(w.^2));
  ftu=ftu*sqrt(n/sum(w.^2));
  
  % throw away symmetric (neg freq) part
  nup=ceil((n+1)/2);
  ftv=ftv(1:nup,:);
  ftu=ftu(1:nup,:);

  % scale so result not a function of input length
  ftv=ftv/sqrt(n);
  ftu=ftu/sqrt(n);
  
  % generate auto and cross spectra
  hdd=ftv.*conj(ftv);
  hds=ftv.*conj(ftu);
  hsd=ftu.*conj(ftv);
  hss=ftu.*conj(ftu);
  
  % force result real instead of fold-over-and-add which does
  % the same thing
  hdd=real(hdd);
  hds=real(hds);
  hsd=real(hsd);
  hss=real(hss);

  % store for average
  hdda(:,:,i)=hdd;
  hdsa(:,:,i)=hds;
  hsda(:,:,i)=hsd;
  hssa(:,:,i)=hss;
  
end

% take average over half-scans
sDD=nanmean(hdda,3);
sDS=nanmean(hdsa,3);
sSD=nanmean(hsda,3);
sSS=nanmean(hssa,3);

% convert from ADU^2/samprate to ADU^2/Hz
sDD=sDD/samprate;
sDS=sDS/samprate;
sSD=sSD/samprate;
sSS=sSS/samprate;

% calc the half freq axis
switch rem(n,2)
  case 0
    f=[0:n/2]'*samprate/n;
  case 1
    f=[0:n/2-0.5]'*samprate/n;
end

return
