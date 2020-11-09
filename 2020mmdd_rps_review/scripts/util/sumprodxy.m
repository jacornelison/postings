function y=sumprodxy(x,sigma,rho,ndof)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% y=sumprodxy(x,sigma,rho,ndof)
%
% GPT - 2014 Feb 20
%
% transformed variance-gamma distribution
% related to the sum of the product of normally distributed variables X, Y
%
% inputs:
%   x = points where the distribution is to be evaluated
%   sigma = geometric mean of std(X) and std(Y)
%   rho = correlation coefficient of X and Y
%   ndof = number of degrees of freedom
%
% If rho=1 (i.e. X=Y), this is the chi^2 distribution.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[errorcode,x,sigma,rho,ndof] = distchck(4,x,sigma,rho,ndof);
if errorcode > 0
  error('Requires non-scalar arguments to match in size');
end

x(x==0)=eps;
sigma(sigma<1e3*eps)=1e3*eps;

y=zeros(size(x));

% round to chi^2 if correlation is very high
v=find(rho>0.999);
y(v)=chi2pdf(x(v)./sigma(v).^2,ndof(v))./sigma(v).^2;
v=find(rho<-0.999);
y(v)=chi2pdf(-x(v)./sigma(v).^2,ndof(v))./sigma(v).^2;

v=find(y==0);
if isempty(v)
  % everything is chi^2 and we are done
  return
end

% normalization factor
normfac=sqrt(pi*2.^(ndof-1).*(1-rho.^2)).*gamma(ndof/2).*sigma.^(ndof+1);
% argument of the exponential term
u=rho.*x./(sigma.^2.*(1-rho.^2));
% argument of the besselk term
z=abs(x)./(sigma.^2.*(1-rho.^2));

% try Matlab internal functions first for good behavior near the origin
y(v)=(1./normfac(v)).*abs(x(v)).^((ndof(v)-1)/2).*exp(u(v)).*besselk((ndof(v)-1)/2,z(v));

v=find(y<=0|~isreal(y)|~isfinite(y));

% asymptotic series expansion of logarithm of besselk
N=ceil(nanmax(ndof(:)))+10; % plenty of terms
c=get_besselk_coeff((ndof-1)/2,N);
besselk_series=ones(size(y));
for jj=1:N
  % converge the asymptotic series as well as we can
  this_term=real(exp(c{jj}(v)-jj*log(z(v))));
  this_term(~isfinite(this_term))=0;
  besselk_series(v)=besselk_series(v)+this_term;
end

% logarithm of pdf
logy=zeros(size(y));
logy(v)=-(ndof(v)/2)*log(2)-ndof(v).*log(sigma(v))-gammaln(ndof(v)/2)...
    +((ndof(v)-2)/2).*log(abs(x(v)))+u(v)-z(v)+log(besselk_series(v));
% re-exponentiate
y(v)=exp(logy(v));

y(y<0|~isreal(y)|~isfinite(y))=0;

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function coeff_table=get_besselk_coeff(nu,N)

% see Abramowitz & Stegun 9.7.2
% actually calculate the logarithm

coeff_table=repmat({zeros(size(nu))},[N,1]);

mu=4*nu.^2;
coeff_table{1}=log((mu-1)/8);

for jj=2:N
  coeff_table{jj}=coeff_table{jj-1}+log(mu-(2*jj-1)^2)-log(8*jj);
end

return