function [simopt p]=get_dithers(simopt,ind)
% [simopt p]=get_dithers(simopt,ind)
%
% Generate dither parameters that will be added to output of get_array_info. The fields
% of output p are:
%
% rndeps - Nchan x Nrlz
% rndchi - Nchan x Nrlz
% rndcen - Nchan x Nrlz x 2(ra/dec) of B
% rndwid - Nchan x Nrlz x 2(fwhm_maj/fwhm_min) of B
% rndwid - Nchan x Nrlz of B
% abgain - Nchan x Nrlz
%
% If these fields in simopt, have size 2x2 then random numbers are
% generated and the full Nchan x Nrlz x (1 or 2) are created.
% If these fields already have length Nchan, nothing is done.

% record state or random number generator
if ~isfield(simopt,'ditherstate')
  simopt.ditherstate=randn('state');
end

% epsilon
simopt=get_dither_simopt(simopt,ind,'rndeps');
if isfield(simopt,'rndeps')
  p.epsilon=simopt.rndeps;
end

% chi (pol angle)
simopt=get_dither_simopt(simopt,ind,'rndchi');
if isfield(simopt,'rndchi')
  p.chi=simopt.rndchi;
end

% rndcen of B
simopt=get_dither_simopt(simopt,ind,'rndcen');
if isfield(simopt,'rndcen')
  p.ra_off_dos=simopt.rndcen(:,:,1);
  p.dec_off=simopt.rndcen(:,:,2);
end

% fwhm_maj and fwhm_min of B % prepare for removal (STF)
%  simopt=get_dither_simopt(simopt,ind,'rndwid');
%  if isfield(simopt,'rndwid')
%    p.fwhm_maj=simopt.rndwid(:,:,1);
%    p.fwhm_min=simopt.rndwid(:,:,2);
%  end

% alpha (ellipse angle) of B % prepare for removal (STF)
%  simopt=get_dither_simopt(simopt,ind,'rndalpha');
%  if isfield(simopt,'rndalpha')
%    p.alpha=simopt.rndalpha;
%  end

% sigma of sigma,c,p
simopt=get_dither_simopt(simopt,ind,'rndsig');
if isfield(simopt,'rndsig')
  p.sigma=simopt.rndsig;
end

% p of sigma,c,p
simopt=get_dither_simopt(simopt,ind,'rndelp');
if isfield(simopt,'rndelp')
  p.p=simopt.rndelp;
end

% c of sigma,c,p
simopt=get_dither_simopt(simopt,ind,'rndelc');
if isfield(simopt,'rndelc')
  p.c=simopt.rndelc;
end

% A/B gain mismatch
simopt=get_dither_simopt(simopt,ind,'abgain');
if isfield(simopt,'abgain')
  p.abgain=simopt.abgain;
end

if ~exist('p','var')
  p=[];
end

return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function simopt=get_dither_simopt(simopt,ind,fld)

if isfield(simopt,fld)
  x=getfield(simopt,fld);
else
  disp(sprintf('field %s not found in simopt - no dithers created',fld));
  return
end

% keep track of the observing freqs for a given experiment
bands=get_bands_from_ind(ind);

% only alter simopt if the dither params aren't already specified for each detector
% separately 
if size(x,1) ~= numel(ind.e)
    
  if any(rvec(x))
    % if doing dithering, expand to multiple realizations
    nrlz=length(simopt.rlz);
  else
    % if not, don't
    nrlz=1;
  end
  
  % special case for abgain, for which [1,0] indicates no dithering
  if strcmp(fld,'abgain')
    if ~any(rvec(bsxfun(@minus, x, [1,0])))
      nrlz=1;
    end
  end
  
  switch fld
   case 'rndeps'
    simopt.rndeps=zeros(numel(ind.e),nrlz);
    for ff=1:length(bands)
      simopt.rndeps(ind.(['l' bands{ff}]),:)=repmat(x(ff,1),numel(ind.(['l' bands{ff}])),nrlz)+x(ff,2)*randn(numel(ind.(['l' bands{ff}])),nrlz);
    end
   case 'rndchi'
    simopt.rndchi=zeros(numel(ind.e),nrlz);
    for ff=1:length(bands)
      simopt.rndchi(ind.(['l' bands{ff}]),:)=repmat(x(ff,1),numel(ind.(['l' bands{ff}])),nrlz)+x(ff,2)*randn(numel(ind.(['l' bands{ff}])),nrlz);
    end
    
   case 'rndcen'
    simopt.rndcen=zeros(numel(ind.e),nrlz,2);
    simopt.rndcen(ind.lb,:,1)=repmat(x(1,1),numel(ind.lb),nrlz)+x(1,2)*randn(numel(ind.lb),nrlz);
    simopt.rndcen(ind.lb,:,2)=repmat(x(2,1),numel(ind.lb),nrlz)+x(2,2)*randn(numel(ind.lb),nrlz);
    % Copy over into A pairs to dither pair centroid, leaving differential pointing unchanged
    simopt.rndcen(ind.la,:,1)=simopt.rndcen(ind.lb,:,1);
    simopt.rndcen(ind.la,:,2)=simopt.rndcen(ind.lb,:,2);
    
%     case 'rndwid' %prepare for removal (STF)
%      simopt.rndwid=zeros(numel(ind.e),nrlz,2);
%      simopt.rndwid(ind.lb,:,1)=repmat(x(1,1),numel(ind.lb),nrlz)+x(1,2)*randn(numel(ind.lb),nrlz);
%      simopt.rndwid(ind.lb,:,2)=repmat(x(2,1),numel(ind.lb),nrlz)+x(2,2)*randn(numel(ind.lb),nrlz);    

%     case 'rndalpha' %prepare for removal (STF)
%      simopt.rndalpha=zeros(numel(ind.e),nrlz);
%      simopt.rndalpha(ind.lb,:,1)=repmat(x(1,1),numel(ind.lb),nrlz)+x(1,2)*randn(numel(ind.lb),nrlz);

   case 'rndsig'
    simopt.rndsig=zeros(numel(ind.e),nrlz);
    simopt.rndsig(ind.la,:,1)=repmat(x(1,1),numel(ind.la),nrlz)+x(1,2)*randn(numel(ind.la),nrlz);
    % keep A = B beam when only the non-differential beamwidth is given, e.g. [1,0.01;0,0]
    simopt.rndsig(ind.lb,:,1)=simopt.rndsig(ind.la,:,1);
    % convert the \delta\sigma = (sA^2-sB^2)/(sA^2+sB^2) for the differential beamwidth into
    % \delta\sigma' = (sA-sB)/(sA+sB) which can easily be applied
    % catch the case where the calculation of the meandiff needs a limit discussion and set it to the proper answer:
    if(x(2,1)==0)
      meandiff = 0;
      sigdiff  = x(2,2)/2;
    else % use the regular caclulation
      meandiff = (1-sqrt(1-x(2,1)^2))/x(2,1);
      sigdiff  = (1-sqrt(1-x(2,1)^2))/(x(2,1)^2*sqrt(1-x(2,1)^2)) * x(2,2);
    end
    randiff  = repmat(meandiff,numel(ind.la),nrlz)+sigdiff*randn(numel(ind.la),nrlz);
    simopt.rndsig(ind.lb,:,2) = 1-randiff;
    simopt.rndsig(ind.la,:,2) = 1+randiff;    

   case 'rndelp'
    simopt.rndelp=zeros(numel(ind.e),nrlz);
    simopt.rndelp(ind.la,:)=repmat(x(1,1),numel(ind.la),nrlz)+x(1,2)*randn(numel(ind.la),nrlz);
    simopt.rndelp(ind.lb,:)=simopt.rndelp(ind.la,:,1);
    randiff = repmat(x(2,1),numel(ind.la),nrlz)+x(2,2)*randn(numel(ind.la),nrlz);
    simopt.rndelp(ind.la,:)= simopt.rndelp(ind.la,:,1)+randiff;
    simopt.rndelp(ind.lb,:)= simopt.rndelp(ind.lb,:,1)-randiff;

   case 'rndelc'
    simopt.rndelc=zeros(numel(ind.e),nrlz);
    simopt.rndelc(ind.la,:)=repmat(x(1,1),numel(ind.la),nrlz)+x(1,2)*randn(numel(ind.la),nrlz);
    simopt.rndelc(ind.lb,:)=simopt.rndelc(ind.la,:);
    randiff = repmat(x(2,1),numel(ind.la),nrlz)+x(2,2)*randn(numel(ind.la),nrlz);
    simopt.rndelc(ind.la,:)= simopt.rndelc(ind.la,:)+randiff;
    simopt.rndelc(ind.lb,:)= simopt.rndelc(ind.lb,:)-randiff;
    
   case 'abgain'
    % We want to set ab gains such that a/b is distributed according to a
    % Gaussian with given mean and sigma, and we want the overall gain
    % within a pair unchanged:
    % a/b=r and (a+b)/2=1 so a=2/(1+1/r) and b=2/(1+r)
    simopt.abgain=ones(numel(ind.e),nrlz);
    for ff=1:length(bands)
      r=repmat(x(ff,1),numel(ind.(['l' bands{ff} 'a'])),nrlz)+x(ff,2)*randn(numel(ind.(['l' bands{ff} 'a'])),nrlz);
      simopt.abgain(ind.(['l' bands{ff} 'a']),:)=2./(1+1./r);
      simopt.abgain(ind.(['l' bands{ff} 'b']),:)=2./(1+r);
    end
  end
  
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function bands=get_bands_from_ind(ind)
% Keep track of the observing freqs for a given experiment
% Get the frequency bands from ind

indfieldnames=fieldnames(ind);
bands={};
for ii=1:length(indfieldnames)
  % look for which frequencies were used by looking for fields
  % such as .l100 .l150 etc.
  % Take the frequency values and store them
  band_used=regexp(indfieldnames{ii},'^l(\d+)$','tokens');
  if ~isempty(band_used)
    bands{end+1}=band_used{1}{1};
  end
end

return

