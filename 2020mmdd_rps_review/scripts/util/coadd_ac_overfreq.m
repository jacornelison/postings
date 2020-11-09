function acc=coadd_ac_overfreq(ac,coaddopt)
% function ac=coadd_ac_overfreq(ac,coaddopt)
%
% written to coadd saved acs into frequencies
% no channel/rx selection is done
% 
% example:
% load('1315/real_2012_filtp3_weight3_gs_dp1100_jack01.mat')
% ac=coadd_ac_overfreq(ac)

if ~exist('ac','var')
  error('Need to input ac!');
end
if ~exist('coaddopt','var')
  error('Need to input coaddopt!');
end

% coaddopt from maps combined over different observations
if iscell(coaddopt)
  coaddopt=coaddopt{1};
else
  coaddopt=coaddopt;
end

% Check for coaddtypes...
% if coaddtype is 0, no need to do anything
if coaddopt.coaddtype==0
  acc=ac;
  return
end

%combcomap just concatenates ac, so have to loop through those
if iscell(ac)
  for kk=1:length(ac)
    switch coaddopt.coaddtype
    case 5
      acc{kk}=addovertile(ac{kk},coaddopt.p,coaddopt.ind);
    case {2,3}
      acc{kk}=addoverpixel(ac{kk},coaddopt.p,coaddopt.ind);
    case 1
      acc{kk}=addoverrx(ac{kk},coaddopt.p,coaddopt.ind);
    case 4
      acc{kk}=addoverdet(ac{kk},coaddopt.p,coaddopt.ind);
    end
  end
else
  switch coaddopt.coaddtype
  case 5
    acc=addovertile(ac,coaddopt.p,coaddopt.ind);
  case {2,3}
    acc=addoverpixel(ac,coaddopt.p,coaddopt.ind);
  case 1
    acc=addoverrx(ac,coaddopt.p,coaddopt.ind);
  case 4
    acc=addoverdet(ac,coaddopt.p,coaddopt.ind);
  end
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
function acc=addovertile(ac,p,ind)

%find number of tiles
ntile=length(unique(p.tile(ind.a)));
%combine to a single index
tiles=p.rx(ind.la)*ntile+p.tile(ind.la); % light because darks have NaN tile  
ntile=length(unique(p.tile(ind.a)));
freqs=unique(p.band(ind.la)); %light a because the darks have 0 band

for nf=1:length(freqs)

  %find which tiles are at the frequency
  tile=unique(tiles(p.band(ind.la)==freqs(nf)));

  %initialize ac
  acc(nf,:)=ac(tile(1),:);

  %add in each of the following rxs for each jackhalf
  for jj=1:size(ac,2)
    for ii=tile(2:end)
      acc(rx,jj)=addac(acc(rx,jj),ac(ii,jj));
    end
  end
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function acc=addoverpixel(ac,p,ind)

%find the rxs
freqs=unique(p.band(ind.la)); %light a because the darks have 0 band

for nf=1:length(freqs)

  %find which pixels are at the freqeuncy
  pix=find(p.band(ind.a)==freqs(nf));
  %initialize ac
  acc(nf,:)=ac(pix(1),:);

  %add in each of the following pixels for each jackhalf
  for jj=1:size(ac,2)
    for ii=pix(2:end)'
      acc(nf,jj)=addac_nan(acc(nf,jj),ac(ii,jj));
    end
  end
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function acc=addoverrx(ac,p,ind)

%find the rxs
freqs=unique(p.band(ind.la)); %light a because the darks have 0 band

for nf=1:length(freqs)

  %find which rxs are at the frequency
  rx=unique(p.rx(p.band==freqs(nf)));
  %rx is indexed from 0, fix it for matlab to start from 1:
  rx = rx+1;
  %initialize ac
  acc(nf,:)=ac(rx(1),:);

  %add in each of the following pixels for each jackhalf
  for jj=1:size(ac,2)
    if numel(rx) > 1
      for ii=rx(2:end)'
        acc(nf,jj)=addac(acc(nf,jj),ac(ii,jj));
      end
    end
  end
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function acc=addoverdet(ac,p,ind)

%find the rxs
freqs=unique(p.band(ind.la)); %light a because the darks have 0 band
acord=[ind.a ind.b]; %kind of strange ordering the ac is saved in

for nf=1:length(freqs)

  %find which dets are at the frequency
  deta=find(p.band(ind.a)==freqs(nf));
  [lia detloc]=ismember([ind.a(deta) ind.b(deta)],acord);
  
  %initalize ac
  acc(nf,:)=ac(detloc(1),:);

  %add in each of the following pixels for each jackhalf
  for jj=1:size(ac,2)
    for ii=detloc(2:end)
      acc(nf,jj)=addac_nan(acc(nf,jj),ac(ii,jj));
    end
  end
end

return






