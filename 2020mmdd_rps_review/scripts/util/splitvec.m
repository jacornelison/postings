function x2=splitvec(x,n)
% x=splitvec(v,n)
%
% For vector v, split up into chunks of size n and return as cell array of vectors.
% i.e. x=splitvec(1:100,7)

nx=numel(x);
clear x2

k=1;
s=1;
e=n;

nx2=ceil(nx/n);
for k=1:nx2
  x2{k}=x(s:min(e,nx));
  s=s+n;e=e+n;
end

return
