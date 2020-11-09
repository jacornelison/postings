function r=sumprodxyrnd(sigma,rho,ndof,sizeOut)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% r=sumprodxyrnd(sigma,rho,ndof,sizeOut)
%
% GPT - 2014 Feb 24
%
% random samples from a transformed variance-gamma distribution
% related to the sum of the product of normally distributed variables X, Y
%
% inputs:
%   sigma = geometric mean of std(X) and std(Y)
%   rho = correlation coefficient of X and Y
%   ndof = number of degrees of freedom
%   sizeOut = [m,n,...] size of output matrix (default=[1,1])
%
% If rho=1 (i.e. X=Y), this is the chi^2 distribution.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% generate random chi^2
c=chi2rnd(ndof,sizeOut);

% normal variance-mean mixture
r=sigma^2*(rho*c+sqrt(1-rho^2)*sqrt(c).*randn(sizeOut));