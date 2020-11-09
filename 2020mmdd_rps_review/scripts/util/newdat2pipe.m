function [r,bpwf]=newdat2pipe(ndfile,beamerrfile,abscalerr)
%
% take a newdat file and convert to pipeline 'r' structure 
%
% input: ndfile = newdat file
%        beamerrfile = optional file with beam err vs l
%        abscalerr = optional vec specifying gain err
%                  abscalerr=[TT TE EE BB TB EB]
%
% output: pipeline real data structure r
%         optional bpwf output in usual format
%
% r=newdat2pipe('newdat/BICEP_20090618.newdat');
% r=newdat2pipe('newdat/BICEP_20090618.newdat','newdat/BICEP_20090618_beam_err.txt')
% r=newdat2pipe('newdat/BICEP_20090618.newdat','newdat/BICEP_20090618_beam_err.txt',[.04 .045 .056 .056 .045 .056])

% read in newdat file and put in bicep spectra order
dat=read_newdat(ndfile,1);

% set flag in r struct so we know it's from a newdat file
r.newdat=1;

nspec=length(dat.ss);
nbins=dat.nbins(1);

% add l bins
r.l=dat.lm{1}';

% add real and derr
% add bpcov & corr matrices
% add xfactors
for ii=1:nspec
  r.real(:,ii)=dat.Cl{ii}';
  % in principle the upper/lower could be different, take mean
  r.derr(:,ii)=mean([dat.Cleu{ii}' dat.Clel{ii}'],2);

  r.cov{ii}=dat.bpcm{ii};
  r.corr{ii}=dat.corr{ii};

  r.xfac(:,ii)=dat.xfac{ii}';
end

% add beam/gain errors, use single value from newdat
% if not specified in input arguments
% pass beam width & uncert to add_sys_uncer to calculate
if ~exist('beamerrfile','var')
  bw=dat.beamsize/60;  % convert from arcmin to deg
  bwu=dat.beamuncer/dat.beamsize;   % convert to %err
  dum=add_sys_uncer(r,bw,bwu);
  be=dum.beam_uncer;
else
  beam=textread(beamerrfile,'','commentstyle','shell');
  if sum(r.l==mean([beam(:,1) beam(:,2)],2))~=length(r.l)
    disp('beam file contains different bands than newdat file!')
    return
  end    
  be=repmat(beam(:,3),1,nspec);
end 
if ~exist('abscalerr','var')
  abse=repmat(dat.abscaluncer,nbins,nspec);
else
  abse=repmat(abscalerr,nbins,1);
end
r.beam_uncer=be;
r.abscal_uncer=abse;

if nargout>1
  % ret bpwf in format n_ell x bins x (TT/TE/EE/BB/(EE->BB)/(BB->EE)
  % input bpwf from read_newdat in format dat.bpwf struct with 
  % (nbins*nspec) fields, (nbinsTT)(nbinsTE)(nbinsEE)(nbinsBB)
  % each field has 5 columns of data, [ell TT TE EE BB]
  bpwf=[];
  % first get the TT/TE/EE/BB bandpowers
  for ii=1:4
    for jj=1:nbins
      bpwf.l(:,jj,ii)=dat.bpwf{jj+nbins*(ii-1)}.l;
      bpwf.Cs_l(:,jj,ii)=dat.bpwf{jj+nbins*(ii-1)}.Cs_l(:,ii);
    end
  end
  % if EE->BB & BB->EE exist, they should be in this config (I think)
  for jj=1:nbins
    % EE->BB
    bpwf.Cs_l(:,jj,5)=dat.bpwf{jj+nbins*(4-1)}.Cs_l(:,3);
    % BB->EE
    bpwf.Cs_l(:,jj,6)=dat.bpwf{jj+nbins*(3-1)}.Cs_l(:,4);
  end
  bpwf.l=dat.bpwf{1}.l;
end

return

