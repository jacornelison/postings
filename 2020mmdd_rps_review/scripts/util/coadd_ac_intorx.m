function acc=coadd_ac_intorx(ac,coaddopt)
% function ac=coadd_ac_intorx(ac,coaddopt)
%
% written to coadd saved acs over receiver
% no channel/rx selection is done
% assumes type5 keck
% 
% example:
% load('1315/real_2012_filtp3_weight3_gs_dp1100_jack01.mat')
% ac=coadd_ac_intorx(ac)

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
% if coaddtype is 0 or 1, no need to do anything
if coaddopt.coaddtype==0 || coaddopt.coaddtype==1
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
    end
  end
else
  switch coaddopt.coaddtype
  case 5
    acc=addovertile(ac,coaddopt.p,coaddopt.ind);
  case {2,3}
    acc=addoverpixel(ac,coaddopt.p,coaddopt.ind);
  end
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
function acc=addovertile(ac,p,ind)

%find number of tiles/rxs
ntile=length(unique(p.tile(ind.a)));
nrx=length(unique(p.rx(ind.a)));

for rx=1:nrx
  %initialize ac
  rxb=(rx-1)*ntile; %first tile.  assume rx and tile are in order...
  acc(rx,:)=ac(rxb+1,:);

  %add in each of the following rxs for each jackhalf
  for jj=1:size(ac,2)
    for ii=rxb+2:rxb+4
      acc(rx,jj)=addac(acc(rx,jj),ac(ii,jj));
    end
  end
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function acc=addoverpixel(ac,p,ind)

%find the rxs
rxs=unique(p.rx(ind.a))';

for rx=rxs
  %find which pixels are in teh rx
  pix=find(p.rx(ind.a)==rx);
  acc(rx+1,:)=ac(pix(1),:);

  %add in each of the following pixels for each jackhalf
  for jj=1:size(ac,2)
    for ii=pix(2:end)'
      acc(rx+1,jj)=addac_nan(acc(rx+1,jj),ac(ii,jj));
    end
  end
end

return








