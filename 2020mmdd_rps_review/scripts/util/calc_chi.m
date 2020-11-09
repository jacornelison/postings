function r=calc_chi(r,chibins,n,meansim)
% r=calc_chi(r,chibins,n,meansim)
%
% Calculate chisq of real and sim spectra
%
% It's always better to recalc cov mat excluding each sim when calc
% chi2 of each sim - this makes the sims directly comparable to
% data. It does mean that both move outside the theory range when
% including the full cov mat and nsim reasonable. When including only
% main + 2 offdiags the difference is small for nsim>100.
%
% If x factor exist they will be used for TT/EE/BB
%
% If abscal and beam uncertainties exist they will be used to make
% appropriate additions to the bandpower covariance matrix
%
% e.g:
% load final/0751/real_a_filtp3_weight3_gs_dp1102_jack0_matrix_directbpwf_rbc
% chibins=get_chibins();
% r=calc_chi(r,chibins,[],1)
% r=add_sys_uncer(r,[],1)
% r=calc_chi(r,chibins)

if(~exist('chibins','var'))
  chibins=[];
end
if(~exist('n','var'))
  n=[];
end
if(~exist('meansim','var'))
  meansim=[];
end

if(isempty(chibins))
  for i=1:9
    chibins{i}=1:size(r(1).real,1);
  end
end
if(isempty(n))
  n=1; % changed 2>1 by CLP 140205 - this will change all numbers a
       % bit!
end
if(isempty(meansim))
  meansim=false;
end

if(meansim)
  % set expectation value to mean of s+n sims
  for j=1:length(r)
    r(j).expv=nanmean(r(j).sim,3);
  end
end

% should probably use add_sys_uncer to pre-calc the diag stripped,
% syserr augmented cov mat before entering this loop - should be
% identical result but would avoid code repetition

for j=1:length(r)
  for i=1:size(r(j).real,2)
  
    % take chisq for real spectrum
    d=r(j).real(chibins{i},i);
    
    % model expectation values
    m=r(j).expv(chibins{i},i);

    % get the cov mat of s+n realizations
    s=squeeze(r(j).sim(chibins{i},i,:));
    c=nancov(s');

    % if wanted strip to diag region
    if(n>=0)
      c=triu(tril(c,n),-n);
    end
    
    % add abscal and beam uncertainties to cov mat
    % (and remember gu and su for use with sims below)
    % use the model to make the addition    
    if(isfield(r(j),'abscal_uncer'))
      gu=r(j).abscal_uncer(chibins{i},i);
      c=c+(cvec(gu)*rvec(gu)).*(cvec(m)*rvec(m));
    end
    if(isfield(r(j),'beam_uncer'))
      su=r(j).beam_uncer(chibins{i},i);
      c=c+(cvec(su)*rvec(su)).*(cvec(m)*rvec(m));
    end
    %[s,co]=cov2corr(c); co
    
    % take inverse of cov mat
    c=inv(c);
      
    % apply offset lognormal transformation for TT/EE/BB
    if(ismember(i,[1,3,4]) & isfield(r(j),'xfac'))
      [d,m,c]=xform_cov(d,m,c,r(j).xfac(chibins{i},i));  
    end
    
    % finally take the chi2
    r(j).rchisq(i)=rvec((d-m))*c*cvec((d-m)); 
   
    % calc likelihood
    %r(j).rlike(i)=(1/sqrt(det(inv(c))))*exp(-r(j).rchisq(i)/2);
    r(j).rlike(i)=exp(-r(j).rchisq(i)/2);
    
    % take chisq for each sig+noi sim
    if(isfield(r,'sim'))
      for k=1:size(r(j).sim,3)
    
        % take chisq for sim spectrum
        d=r(j).sim(chibins{i},i,k);
        
        % model expectation values
        m=r(j).expv(chibins{i},i);
        
        % recalc cov mat from sig+noi sims excluding the realization we
        % are about to calc chi2 for - this is then properly equivalent
        % to what we do with real data
        % get the sims excluding this one
        x=true(1,size(r(j).sim,3));
        x(k)=false;
        s=squeeze(r(j).sim(chibins{i},i,x));
        % take covariance
        c=nancov(s');
        
        % if wanted strip to diag region
        if(n>=0)
          c=triu(tril(c,n),-n);
        end
        
        % add abscal and beam uncertainties to cov mat
        % NB: does this make sense? shouldn't we only "increase the
        % error bars" for the real data? If we do it for sims as
        % well won't the pte versus sims stay basically the same?
        if(exist('gu','var'))
          c=c+(cvec(gu)*rvec(gu)).*(cvec(m)*rvec(m));
        end
        if(exist('su','var'))
          c=c+(cvec(su)*rvec(su)).*(cvec(m)*rvec(m));
        end
        
        % take inverse of cov mat
        c=inv(c);
        
        % apply offset lognormal transformation for TT/EE/BB
        if(ismember(i,[1,3,4]) & isfield(r(j),'xfac'))
          [d,m,c]=xform_cov(d,m,c,r(j).xfac(chibins{i},i));
        end
        
        % finally take the chi2
        r(j).schisq(i,k)=rvec((d-m))*c*cvec((d-m));
        
        % calc likelihood
        %r(j).slike(i,k)=(1/sqrt(det(inv(c))))*exp(-r(j).schisq(i,k)/2);
        r(j).slike(i,k)=exp(-r(j).schisq(i,k)/2);

      end
    end
    
    % calc probability to exceed versus theoretical chisq distribution
    
    % the "proper" way
    %r(j).ptet(i)=1-chi2cdf(r(j).rchisq(i),length(chibins{i}));
    
    % This should be equivalent for sensible chi2 and continues to
    % give non-zero value for outrageously large values
    r(j).ptet(i)=sum(chi2pdf(r(j).rchisq(i):0.01:1000,length(chibins{i})))/100;
    
    % and again for simulation set
    if(isfield(r,'sim'))
      r(j).ptes(i)=sum(r(j).rchisq(i)<r(j).schisq(i,:))/size(r(j).sim,3);
    end
    
    % convert chi2 pte to "sigmas"
    % - this is two sided defn as apparently used in Kovac et al
    % paper - PDG may prefer one sided... need to check
    % (To reproduce page 11 -norminv(8.46e-7/2,0,1)=4.92)
    r(j).sigmat(i)=-norminv(r(j).ptet(i)/2,0,1);
    
    % copy chibins
    r(j).chibins{i}=chibins{i};
  end
end

% following multi spectrum code currently disabled - needs fixup
% before use again so matches the algorithm above
if(0)
% do it again for multiple spectra at once
for j=1:length(r)
    
  % get the cov mat of s+n realizations
  % WARNING: This appears to assume chibins is the same for all spectra!
  sz=size(r(j).sim(chibins{i},1:4,:));
  s=reshape(r(j).sim(chibins{i},1:4,:),[prod(sz(1:2)),sz(3)]);
  c=nancov(s');
  
  % if wanted strip to diag region
  if(n>=0)
    c=zero_offdiags(n,4,c);
  end
  
  % take inverse
  c=inv(c);
  
  % get data and model values
  d=r(j).real(chibins{i},1:4);
  m=r(j).expv(chibins{i},1:4);

  % take chisq for real spectrum
  r(j).rchisqa=rvec((d-m))*c*cvec((d-m));
  % record the ndf
  r(j).chisqa_ndf=length(c);
  
  % take chisq for each sig+noi sim
  for k=1:size(r(j).sim,3)
    
    % recalc cov mat from sig+noi sims excluding the realization we
    % are about to calc chi2 for - this is then properly equivalent
    % to what we do with real data
    
    % get the sims excluding this one
    x=true(1,size(r(j).sim,3));
    x(k)=false;
    s=reshape(r(j).sim(chibins{i},1:4,x),[prod(sz(1:2)),sz(3)-1]);
    
    c=nancov(s');
    
    % if wanted strip to diag region
    if(n>=0)
      c=zero_offdiags(n,4,c);
    end
    
    c=inv(c);
      
    d=r(j).sim(chibins{i},1:4,k);    
    r(j).schisqa(k)=rvec((d-m))*c*cvec((d-m));
  end
    
  % calc probability to exceed versus theoretical chisq distribution
  r(j).pteat=1-chi2cdf(r(j).rchisqa,r(j).chisqa_ndf);
  
  % and again versus simulation set
  r(j).pteas=sum(r(j).rchisqa<r(j).schisqa)/size(r(j).sim,3);
end
end

return

% do offset-lognormal transformation
function [d,m,c]=xform_cov(d,m,c,x);

c=c.*(cvec(d+x)*rvec(d+x));
m=real(log(m+x));
d=real(log(d+x));

return

function bpcmp=zero_offdiags(nod,nspec,bpcm)
% bpcmp=zero_offdiags(nod,nspec,bpcm)
%
% Force offdiags to zero in multispectrum band power covariance matrix
%
% nod = "n off diags" to keep
% nspec = "number of spectra" in bpcm
%
% Note input is TT/TE/EE/BB

n=length(bpcm)/nspec;

bpcmp=zeros(size(bpcm));

for i=0:(nspec-1)
  for j=0:(nspec-1)
    % pull out the submatrix
    x=bpcm(i*n+1:(i+1)*n,j*n+1:(j+1)*n);
    
    y=zeros(size(x));

    % only select secondary off diags
      
    % keep intra spectra main + n
    if(j==i)
      y=triu(tril(x,nod),-nod);
    end
    
    % keep TT/TE and TE/EE main + 1 only
    if((i==1 & j==0) | (i==2 & j==1) | (j==1 & i==0) | (j==2 & i==1))
      y=triu(tril(x,1),-1);
    end
    
    % re-insert sub matrix
    bpcmp(i*n+1:(i+1)*n,j*n+1:(j+1)*n)=y;
  end
end
