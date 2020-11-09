function cons1d_stats(x,y)
% cons1d_stats(x,y)
%
% Given2 a 1d constraint curve x,y print stats to screen
%
% This script used for BKP paper - see hpd_interval also
%
% e.g.:
% x=0:0.1:10; y=gauss([1,3,1],x); cons1d_stats(x,y)

if(any(x<0))
  error('x must start from zero and increase');
end

if(x(1)~=0)
  warning('x doesnt start from zero - attempting to pad');
  d=x(2)-x(1);
  xx=0:d:x(1);
  xx=xx(1:end-1);
  x=[xx,rvec(x)]; y=[zeros(size(xx)),rvec(y)];
end

% look up the zero to peak ratio
r=y(1)/max(y);

% this is Wilk's theorem - see:
% http://en.wikipedia.org/wiki/Likelihood-ratio_test#Distribution:_Wilks.27s_theorem
chi2=-2*log(r);

% For a likelihood ratio of zero this gives chi2=0 which has PTE 1
% - this is correct when the observed value can fall either
% side of the hypothesis value - the chance of the likelihood curve
% peaking exactly at the true value is zero
% However if the hypothesis is that the true value is zero then
% even if the ratio is epsilon less than one then the PTE is
% already less than 50% - half of the likelihood curves in the null
% case peak at zero ("try" to peak below it).
% So it seems that this is correct
pte=0.5*(1-chi2cdf(chi2,1));

% To convert a PTE to sigmas we do the following
sig=norminv(1-pte);
% note that this is different to likeratio2sigmas as used in B2 PRL
% but it seems correct - when converting a PTE to sigmas it ought
% to be possible to get a negative result if PTE>0.5 - of course in
% this case it never can be for a curve which peaks above zero.

disp(sprintf('L0/Lm=%.3e (Lm/L0=%.1f) pte=%.3e sig=%.1f',r,1/r,pte,sig));

% make the cdf
yc=cumsum(y);
yc=yc./yc(end);
% prevent interp1 from failing
[yc,i]=unique(yc); xp=x(i);
% find the 95% point
lim95=interp1(yc,xp,0.95);

disp(sprintf('95 percent upper limit=%.3e',lim95));

% to find the level which encompases 68% of the likelihood it is
% necessary to have a high-res curve

% increase the resolution
xp=linspace(0,x(end),10000);
yp=interp1(x,y,xp,'spline');

% find the peak position
[mv,i]=max(yp); p=xp(i);

% norm to unit integral
yp=yp./sum(yp);
% working down...
l=linspace(max(yp),0,10000);
% find the level which encompasses 68% of the integral
for i=1:length(l)
  if(sum(yp(yp>l(i)))>0.6827) % diff(normcdf([-1,1]))
    break
  end
end
% look up crossing levels
lim=find(yp>l(i));
ll=xp(lim(1)); ul=xp(lim(end));

if(0)
  clf
  plot(xp,yp);
  line(xlim,[l(i),l(i)]);
  yl=ylim;
  line([ll,ll],yl);
  line([ul,ul],yl);
  line([lim95,lim95],yl,'color','r');
end

disp(sprintf('peak posn=%f+%f-%f (%f sym)',...
             p,ul-p,p-ll,(ul-ll)/2));

return
