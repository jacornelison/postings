function cd=get_chol(x)

% Not sure if this is needed anymore
if(all(rvec(x)==0))
  cd=zeros(size(x));
  return
end

% take eigen values/vectors
[V,D] = eig(x);

% set small/negative eigenvalues to a small positive value
D=diag(D);
D=max(D,(max(D)/1e12)*ones(size(D)));
    
% regenerate the matrix which should now be positive definite
z=V*diag(D)*V';

% get rid of small imag component on diag that appears to be caused by
% numerical precision limits
z=(z+z')/2;

% take the chol decomposition
cd=chol(z);
    
% look at effect of "pos-def fixup"
if(0)
  subplot(1,2,1)
  imagesc(x); colorbar
  title('input fourier cov matrix')
  subplot(1,2,2)
  imagesc(x-z); colorbar
  title('delta of input and output') 
end

return
