function r=add_sys_uncer(r,noffdiag,bw,bwu,meansim)
% r=add_sys_uncer(r,bw,bwu)
%
% Adds abscal and beam uncertainty fields to r structure
%
% use one code to calculate beam_uncer, allow input bw & bwu for now
% optional input:  bw (in deg) bwu (% uncert)
%
% 2014-03-11 UPDATED to include abscal and beam width uncertainties
% appropriate for BICEP2 (still hard-coded). See 2014-03-11
% posting. (CAB)

if(~exist('noffdiag','var'))
  noffdiag=1;
end
if(~exist('bw','var'))
  bw=[];
end
if(~exist('bwu','var'))
  bwu=[];
end
if(~exist('meansim','var'))
  meansim=0;
end

if ~isfield(r,'abscal_uncer')
  r.abscal_uncer=zeros(size(r.real));
  for i=1:size(r.real,2)
    switch i
      case 1
        % Chiang et al Sec 9.2 gives different values for TT, TE/TB
        % and EE/BB/EB
        % Values now updated for BICEP2.
        r.abscal_uncer(:,i)=0.026; % in power units
      case {2,5}
        r.abscal_uncer(:,i)=0.027;
      otherwise
        r.abscal_uncer(:,i)=0.027;
    end
  end
end

if ~isfield(r,'beam_uncer')
  % Randol's thesis / instrument paper say 31.22 arcmin FWHM
  if isempty(bw)
    % Convert from FWHM to sigma, and arcmin to degrees.
    bw = 31.22 / 60 / sqrt(8 * log(2));
  end
  if isempty(bwu)
    % Beam width fractional uncertainty.
    % Value updated for BICEP2.
    bwu=0.011;
  end
  bu=exp((bw*pi/180)^2*(bwu^2+2*bwu)*r.l.*(r.l+1))-1;

  % Chiang et al says subtract the mean in range 56<ell<265 due to way
  % in which abscal done
  % Paper doesn't mention taking absolute value although this is done
  % data release file
  % In fact it seems wrong to do so - but stays in for now as we are
  % trying to reproduce numbers in paper

  % We decide Remove the absolute value is correct -- IDB 2013-03-25
  % Copy this correction from B1 to B2 -- IDB 2014-02-25
  % Updated abscal bins to match BICEP2 abscal procedure.
  abscal_bins = [4:9]; % Ell bins that were used to calculate abscal.
  bu=(bu-mean(bu(abscal_bins)));

  r.beam_uncer=zeros(size(r.real));
  for i=1:size(r.real,2)
    r.beam_uncer(:,i)=bu;
  end

  % Also calculate "beam uncertainty weight", which are used if we want to 
  % analytically correct the variance of the rho statistic. See BICEP2 
  % posting from 2014-03-04. This term is equivalent to gamma*epsilon from 
  % the posting.
  r.beam_uncer_weight = -2 * (bw * pi / 180)^2 * bwu * ...
      (r.l .* (r.l + 1) - mean(r.l(abscal_bins) .* (r.l(abscal_bins) + 1)));

end

% add it to the covariance matrix. Always recalculate the cm
% first to avoid adding it several times:
r=get_bpcov(r);
% trim the bpcov and the add the systematic uncertainty which
% is not limited by mc statistic:
if exist('noffdiag','var')
  r=trim_bpcov(r,noffdiag);
end

if(meansim)
  % set expectation value to mean of s+n sims
  for j=1:length(r)
    r(j).expv=mean(r(j).sim,3);
  end
end

% add the uncertainty:
for i=1:size(r.real,2)
  % add sys uncer based on the model
  m=r.expv(:,i);
  
  gu=r.abscal_uncer(:,i);
  r.cov{i}=r.cov{i}+(cvec(gu)*rvec(gu)).*(cvec(m)*rvec(m));
  
  su=r.beam_uncer(:,i);
  r.cov{i}=r.cov{i}+(cvec(su)*rvec(su)).*(cvec(m)*rvec(m));

  % retake diag err bars
  r.derr(:,i)=sqrt(diag(r.cov{i}));
end

return
