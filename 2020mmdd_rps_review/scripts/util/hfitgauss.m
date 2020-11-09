function [p,pe,prob,stat]=hfitgauss(x,n,lims,opt,freepar)
% [p,pe,prob,stat]=hfitgauss(x,n,lims,opt,freepar)
%
% A driver for matmin to make it as easy as possible
% to fit a Gaussian.
% Very similar to hfit except pre-calcs the starting
% values for Gaussian fit
%
% x,n is histogram data
% lims is optional x range to fit over
% opt is plotting options:
%   'L' = log scale
%   '0' = don't draw
% optional freepar specifies which parameters are free
%
% p,pe are parmeters and errors
% stat is the return status of the fit.
% The fit is done log likelihood style considering all bins.
% To allow a rough goodness of fit assessment chi-squared is
% calculated using non-empty bins only, and converted to the
% corresponding point prob on the chi2cdf.
%
% eg: [x,n]=hfill(randn(1,1000),100,-3,3);
%     [p,pe]=hfitgauss(x,n)
%     [p,pe]=hfitgauss(x,n,[-2,2],'L',[1,1,0])

if(~exist('lims','var'))
  lims=[];
end
if(~exist('opt','var'))
  opt=' ';
end
if(~exist('freepar','var'))
  freepar=ones(1,3);
end

% vectorize data (if not already)
x=x(:)';
n=n(:)';

if(isempty(lims))
  lims(1)=x(1);
  lims(2)=x(end);
end

ind=find(x>=lims(1)&x<=lims(2));
xl=x(ind);
nl=n(ind);

% Calc starting values for fit
inpar(1)=max(nl);
N=sum(nl);
if(freepar(2)==0)
  inpar(2)=0;
else
  inpar(2)=(1/N)*sum(xl.*nl);
end
inpar(3)=sqrt((1/N)*sum((xl-inpar(2)).^2.*nl));

% Do the fit
[p,pe,gof,stat]=matmin('logl',inpar,freepar,'gauss',nl,xl);
e=sqrt(nl);
chi2 = chisq(p,'gauss',nl,e,xl);
npt=length(e(e~=0));
redchi=chi2/npt;
prob=chi2cdf(chi2,npt);
disp(sprintf('\nChiSq/npnt = %.2f / %.0f = %.2f (%.1f%s)\n\n',...
    chi2,npt,redchi,prob*100,'%'));

if(~any(opt=='0'))
  hplot(x,n);
  
  span=x(end)-x(1); bw=x(2)-x(1);
  
  hold on;
  
  % Plot the fit function over the range used
  x2=xl(1)-bw/2:span/100:xl(end)+bw/2;
  plot(x2,gauss(p,x2),'r');
  
  % Plot the fit function over the full range dashed
  x2=x(1)-bw/2:span/100:xl(1)-bw/2;
  hold on; plot(x2,gauss(p,x2),'g');
  x2=xl(end)+bw/2:span/100:x(end)+bw/2;
  hold on; plot(x2,gauss(p,x2),'g');
  
  hold off;
  
  if(any(opt=='L'))
    ylim([3e-1,max(n)*2]);
    set(gca,'YScale','log');
  else
    ylim([0,max(n)*1.1]);
  end
end
