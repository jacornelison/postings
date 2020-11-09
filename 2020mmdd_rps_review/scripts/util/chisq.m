function chi2 = chisq(pars,func,y,e,varargin)
% chi2 = chisq(pars,func,y,e,x...)
%
% Calculate the ChiSq between a set of points predicted by function func
% with parameters pars evaluated at x, and the observed values y
% with uncertainties e. 
%  
% e can either be a vector of uncertainties or a covariance matrix.
%
% Note that for fitting count histograms should use hfit to get
% correct Poisson count errors.
%
% eg: x=-5:0.1:5; y=gauss([10,0,1],x); y=y+randn(size(y)); e=ones(size(y));
%     plot(x,y,'.');
%     p=fminsearch('chisq',[10,0,1],[],'gauss',y,e,x)
%     [p,pe]=matmin('chisq',[10,0,1],[],'gauss',y,e,x)
%     hold on; plot(x,gauss(p,x),'r'); hold off
%
% See also logl

% If no data provided
if(isempty(y))
  chi2=0;
  return
end

if(isempty(e))
  e=ones(size(y));
end

% a large uncertainty to make zero uncertainties have
% no effect in the fit:
large_err=1e99;

% the calculation is here split in two cases:
% first, the caculation with sigma uncertainties only
% second, the caculation with the covariance matrix
% we keep this split since the using the covariance matrix
% introduces an overhead when going to huge sets of points
% this is for a vector of sigmas:
if (min(size(e))==1)
  % Zero errors can result when assume e=sqrt(n) but they make no sense
  % so following paw I set the error huge so the point it excluded.
  e(e==0)=large_err;

  % Evaluate the function
  f=feval(func,pars,varargin{:});

  % This is the correct behavior when fitting to complex visibilities...
  if(~isreal(y))
    y=[real(y),imag(y)];  f=[real(f),imag(f)];
    % Allow for real or complex error
    if(~isreal(e))
      e=[real(e),imag(e)];
    else
      e=[e,e];
    end
  end

  % Calculate chi2
  chi2=nansum(((y(:)-f(:))./e(:)).^2);
else
% this is for a covariance matrix:
  
  % Zero errors on the diagonal can result when assume e=sqrt(n) but they 
  % make no sense so following paw I set the error huge so 
  % the point it excluded.
  % this is the diagonal of the covariance matrix:
  de = e(1:size(e,1)+1:end);
  % set it to an equivalent value as in the upper calculation:
  de(de==0)=large_err^2;
  % and fill back onto the diagonal:
  e(1:size(e,1)+1:end) = de;

  % Evaluate the function
  f=feval(func,pars,varargin{:});

  % This is the correct behavior when fitting to complex visibilities...
  if(~isreal(y))
    y=[real(y),imag(y)];  f=[real(f),imag(f)];
    % Allow for real or complex covariances, split up the covariance matrix:
    % not sure if this is meaningful:
    if(~isreal(e))
      e = blkdiag(real(e),imag(e));
    else
      e = blkdiag(e,e);
    end
  end

  % Calculate chi2
  chi2=(y-f)*inv(e)*(y-f)';
end
