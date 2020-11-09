function dat=read_newdat(nd_file,reorderspec)
% dat=read_newdat(newdat,reorderspec)
%
% Read back newdat file
% optional input: 
%                 reorderspec = return dat file in order
%                   of bicep spectra TT TE EE BB TB EB
%
% e.g.
% dat=read_newdat('~dbarkats/public_html/bicep2009/data/BICEP_May09.newdat')

if ~exist('reorderspec','var');
  reorderspec=0;
end

% standard order of bicep spectra, for reorderspec
bspeclist={'TT','TE','EE','BB','TB','EB'};

% open the file
fh=fopen(nd_file,'r');

% get the window file base name
x=fgetl(fh);
wfbn=sscanf(x,'%s');
dat.wfbn=wfbn;

% get the number of bins
x=fgetl(fh);
nbins=sscanf(x,'%f',6);
dat.nbins=nbins;

% read "BAND_SELECTION" line
x=fgetl(fh);
x=sscanf(x,'%s',1);
if(~strcmp(x,'BAND_SELECTION'))
  error('next line should be BAND_SELECTION')
end

% read the band selections
for i=1:6
  x=fgetl(fh);
  x=sscanf(x,'%f',2); binl(i)=x(1); binh(i)=x(2);
end

% read the abs cal line
x=fgetl(fh);
x=sscanf(x,'%f',3);
dat.recalswitch=x(1); dat.recalfac=x(2); dat.abscaluncer=x(3);

% read the beam size/uncer line
x=fgetl(fh);
x=sscanf(x,'%f',3);
dat.beamswitch=x(1); dat.beamsize=x(2); dat.beamuncer=x(3);

% read liketype
x=fgetl(fh);
dat.liketype=sscanf(x,'%f',1);

% find bins in included spectra
binsuse=nbins(nbins~=0);

for i=1:sum(nbins>0)
  
  % check says TT, EE etc
  x=fgetl(fh);
  spectag=sscanf(x,'%s',1);
  %if(~strcmp(spectag,dat.ss{i}))
  %  error('read_newdat failed to read expected spec tag')
  %end
  dat.ss{i}=spectag;
 
  % read lines giving bandpower data
  for j=1:binsuse(i)
    x=fgetl(fh);
    x=sscanf(x,'%f');
    dat.bn{i}(j)=x(1);
    dat.Cl{i}(j)=x(2);
    dat.Clel{i}(j)=x(3);
    dat.Cleu{i}(j)=x(4);
    dat.xfac{i}(j)=x(5);
    dat.ll{i}(j)=x(6);
    dat.lu{i}(j)=x(7);
    dat.lm{i}(j)=mean(x(6:7));
    if(length(x)>7)
      dat.ltyp{i}(j)=x(8);
    else
      dat.ltyp{i}(j)=0;
    end
  end
  
  % now scan the block diagonal normalized cov matrices
  dat.corr{i}=fscanf(fh,'%f',[binsuse(i),binsuse(i)]);
  
  % seems necessary to skip to next line
  x=fgetl(fh);
end

% read the all-with-all bpcm
a=sum(nbins);
bpcm=fscanf(fh,'%f',[a,a]);
dat.bpcmall=bpcm;

% now just return the "diagonal" spectra bpcm
for i=1:sum(nbins>0)
  range=(1+binsuse(i)*(i-1)):(binsuse(i)+binsuse(i)*(i-1));
  dat.bpcm{i}=bpcm(range,range);
end

fclose(fh);

% now read the bpwf

% find the dir name
x=strfind(nd_file,'/');
bpwf_file=nd_file(1:x(end));

for i=1:sum(nbins)
  if exist(sprintf('%s/windows/%s%02d',bpwf_file,wfbn,i))
    %x=textread(sprintf('%s/windows/%s%02d',bpwf_file,wfbn,i));
    x=load(sprintf('%s/windows/%s%02d',bpwf_file,wfbn,i));
  else
    %x=textread(sprintf('%s/windows/%s%d',bpwf_file,wfbn,i));
    x=load(sprintf('%s/windows/%s%d',bpwf_file,wfbn,i));
  end
  dat.bpwf{i}.l=x(:,1); 
  % cosmomc newdat files put bpwf in as w_l^b/l
  % as per Chiang et al eqn 19 normalized as
  % sum((l+.5)/(l+1))*(w_l^b/l)=1
  % to get to Cs_l as per eqn 18
  % Cs_l=(l+.5)/(l+1)*(w_l^b/l)
  dat.bpwf{i}.wlb_div_l=x(:,2:end);
  dat.bpwf{i}.Cs_l=x(:,2:end).*repmat((x(:,1)+.5)./(x(:,1)+1),1,size(x(:,2:end),2));
end

if reorderspec
  nspec=sum(nbins>0);
  for ii=1:length(bspeclist)
    tmp=strmatch(bspeclist{ii},dat.ss);
    if ~isempty(tmp)
      reordnum(ii)=tmp;
    end
  end
  fields=fieldnames(dat);
  for ii=1:length(fields)
    tmpfld=getfield(dat,fields{ii});
    tmpfldnew={};
    if size(tmpfld,2)==nspec
      for jj=1:nspec
        tmpfldnew{jj}=tmpfld{reordnum(jj)};
      end
      dat=setfield(dat,fields{ii},tmpfldnew);
    end
  end
  % reorder the bpwf
  tmp=reshape(dat.bpwf,nbins(1),nspec);
  for ii=1:nspec
    tmpnew(1:nbins(1),ii)=tmp(1:nbins(1),reordnum(ii));
  end
  tmpnew=reshape(tmpnew,1,nbins(1)*nspec);
  dat.bpwf=tmpnew;
end

return
