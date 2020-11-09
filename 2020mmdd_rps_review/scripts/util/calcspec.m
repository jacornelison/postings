function [l,Cs_l,n,be]=calcspec(ad,ft1,ft2,w,bintype)
% [l,Cs_l,n,be]=calcspec(ad,ft1,ft2,w,bintype))
%
% calculate (cross) power spectrum 
%
% ft1 and ft2 are the fourier planes to be crossed
% w is weight map in fourier plane
% bintype specfies annular binning

% To make the assumption that l(l+1)C_l/2pi is piecewise flat we need
% to multiply the F plane by the sqrt of this factor before we take
% the mean square within annuli - not after as I was doing in the
% past. The practical effect of this is very small for narrow bands.
l=2*pi*ad.u_r;
sf=sqrt(l.*(l+1)/(2*pi)*prod(ad.del_u));
ft1=ft1.*sf;
ft2=ft2.*sf;

[be,n]=get_bins(bintype);

%if bintype is phi, then need a phi as ad.u_r rather than ell
if(strcmp(bintype, 'phi'))
  %setup so we remove factor of l from radius
  r=ad.u_r*2*pi;
  ft1=ft1./sqrt(r);
  ft2=ft2./sqrt(r);
  %and replace ad.u_r with phi
  [xx,yy]=meshgrid(ad.u_val{1}.*2*pi, ad.u_val{2}.*2*pi);
  ad.u_r=atan2(yy,xx)./(2*pi);
end

% l=2pi*u
% take real as imag component will sum to zero due to each val
% cancelled by its complex conjugate
[l,Cs_l,n]=calcspec_sub(ad.u_r*2*pi,real(ft1.*conj(ft2)),w,be);

% conv out to column vectors
l=cvec(l);
Cs_l=cvec(Cs_l);

return

function [bincenter,mu,n]=calcspec_sub(x,y,w,binedge)
% [bincenter,mu,n]=calcspec_sub(x,y,w,binedge)
%
% Calc the mean of y within each bin in x defined by arbitrary bin
% edge vector binedge
%
% If weights provided make weighted mean
%
% reuses mapping from x into bins from previous call if
% global variable calcspec_ind exists

% if weight vector empty equal weighting
if(isempty(w))
  w=ones(size(y));
end

% Vectorize data
x=x(:)'; y=y(:)'; w=w(:)';

nbin=length(binedge)-1;

% loop over bins
for i=1:nbin
  ind=find(x>=binedge(i)&x<=binedge(i+1));
  % nansum returns 0 if all input elements are nan, 
  % have sxw be nan instead
  % an all nan map then results in a nan spectrum
  if all(isnan(y(ind)))
    sxw(i)=nan;
  else
    % allow for nan weights to indicate modes not to use
    sxw(i)=nansum(y(ind).*w(ind));
  end
  sw(i)=nansum(w(ind));
  n(i)=sum(~isnan(w(ind)));
  bincenter(i)=mean(binedge(i:(i+1)));
end

% divide by sum of weights
mu=sxw./sw;

return
