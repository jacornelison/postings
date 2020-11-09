function r=diag_bandpow(r,n)
% r=diag_bandpow(r)
%
% diagonalize bandpowers according to their cov mat
% see diagonalization.m for notes on diagonalization of bandpowers
%
% n is number of off diags to include - set to -1 for full cov mat

if(~exist('n','var'))
  n=2;
end

% diagonalize each spectrum individually
for i=1:length(r)
  for j=1:6
    
    % get cov mat for this spec
    C=r(i).cov{j};
    
    % strip to close diag if wanted
    if(n>=0)
      C=triu(tril(C,n),-n);
    end
    
    % take inverse of cov mat
    F=inv(C);
    % xform matrix can be sqrtm or chol etc - see diagonalization.m
    Z=sqrtm(F);
    % make xform matrix "scale preserving" - i.e. set column sums to unity 
    Z=Z./repmat(sum(Z),[size(Z,1),1]);
    
    % apply xform to bandpowers    
    q=squeeze(r(i).sim(:,j,:))';
    q=q*Z;
    r(i).sim(:,j,:)=q';
    
    q=squeeze(r(i).real(:,j,:))';
    q=q*Z;
    r(i).real(:,j,:)=q';
    
    % apply xform to cov mat
    C=Z'*C*Z;
    r(i).cov{j}=C;
    
    % recalc diag errors
    r(i).derr(:,j)=sqrt(diag(C));
  end
end

return
