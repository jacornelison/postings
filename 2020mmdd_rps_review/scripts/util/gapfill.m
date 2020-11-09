function d=gapfill(d,fs,ind,scancheck,pl)
% d=gapfill(d,fs,ind,scancheck,pl)
%
% Gap-fill TODs using Cynthia's algorithm
% pl=1 means plot

disp(sprintf('gapfilling...'));

% plot?
if(~exist('pl','var'))
  pl=0;
end

if(~exist('scancheck','var'))
  scancheck=ones(size(d.antenna0.bolo.mag,2),size(fs.sf,1));
end


%For each halfscan
for i=1:length(fs.sf)
  s=fs.sf(i); e=fs.ef(i);
  
  % pull out channels that are both bad (as per scancheck) and called
  % for (as per ind). 
  badscans=find(scancheck(:,i)>0);
  doscans=intersect(badscans,ind);

  if(~isempty(doscans))    
    
    dat=d.antenna0.bolo.mag(s:e,doscans);

    % determine which data are glitches
    [badind,datsm,rms]=get_glitchind(dat,pl);

    % mirror the TOD to avoid edge effects
    dat2=[datsm(end:-1:1,:);datsm;datsm(end:-1:1,:)];
    badind2=[badind(end:-1:1,:);badind;badind(end:-1:1,:)];
    
    % cubic spline interpolate over bad points
    x=cvec(1:size(dat,1));
    x2=cvec(1:size(dat2,1));
    for j=1:size(dat,2)
      xx=x2(~badind2(:,j));
      yy=dat2(~badind2(:,j),j);
      xi=x(badind(:,j));
      xi2=xi+size(dat,1);
      yi=interp1(xx,yy,xi2,'linear');
      temp=dat(:,j);
      n=randn(size(yi))*rms(j);
      dat(xi,j)=yi+n;
      if(j==1 & pl==1)
        figure(1);subplot(3,1,3)
        setwinsize(gcf,900,700)
        plot(temp,'r');hold on;plot(dat(:,j),'blue');title(num2str(doscans(j)));hold off
        
      end
    end   
    
    % replace timestream with gap-filled data
    d.antenna0.bolo.mag(s:e,doscans)=dat;    
  end  
end

return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ind,datsm,rms]=get_glitchind(dat,pl)

% smoothing kernel
ker=ones(1,7);
ker=ker/sum(ker(:));

% pad array with mirrored data the size of the kernel (so filter doesn't
% do the padding with zeros!)
dat=[dat(7:-1:1,:);dat;dat(end:-1:(end-6),:)];
datsm=conv2(dat,cvec(ker),'same');
datsm=datsm(8:1:(end-7),:);
dat=dat(8:1:(end-7),:);

% median referenced rms
dif=dat-datsm;
rms=sqrt(median(dif.^2)-median(dif).^2);

% plotting
if(pl)
  figure(1)
  subplot(3,1,1)
  plot(dat(:,1)); hold on;
  plot(datsm(:,1), 'r');hold off;
  subplot(3,1,2)
  plot(abs(dif(:,1)));hold on; plot([0,size(dat,1)],[5*rms(1),5*rms(1)],'r');hold off;
end

% find where abs(dif)>5*rms and exclude a window of 25 centered on each
% glitch 
ind=(abs(dif)-repmat(5*rms,[size(dat,1),1]))>0;
badind=find(ind==1);
[y,x]=ind2sub(size(dif),badind);
nsample=size(dat,1);
for i=1:length(badind)
  maskind=max(y(i)-12,1):min(y(i)+12,nsample);
  ind(maskind,x(i))=1;
end


return


