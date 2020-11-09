function d=ground_subtract(d,fs,ind,psm)
% d=ground_subtract(d,fs,ind,psm)
%
% function to subtract ground (az) fixed signal from each scanset.  
% use ind.l
%
% ex.
% d=ground_subtract(d,fs,ind.l);
% 
% jul 28 2008 (dB): modified mean to nanmean on line 43,44.

if(~exist('psm','var'))
  psm=[];
end
if(isempty(psm))
  psm=false(size(d.mce0.data.fb));
end

disp('ground_subtract...');

% construct Denis' damn fs.src
if(~isfield(fs,'src'))
  for i=1:length(fs.s)
    fs.src{i,1}=d.antenna0.tracker.source(fs.s(i),:);
  end
end

nscanset=size(unique(fs.set),1);
nch=length(ind);

% for each scanset
for i=1:nscanset

  %find halfscans in that scanset
  q=find(fs.set == i);
  nhalfscans=length(q);
  npoints=fs.ef(q(1))-fs.sf(q(1))+1;

  %find out if the source is az-fixed. 
  if isempty(strfind(fs.src{q(1)},'CMB4')) && ...
     isempty(strfind(fs.src{q(1)},'CMB5')) && ...
     isempty(strfind(fs.src{q(1)},'AZCIRC'))
    disp('not scanning a CMB az fixed region, no GS possible')
    continue % go to next scanset
  end
  
  %cut fs structure to only include this scanset.
  fs_tmp=structcut(fs,q);
  mapind=make_mapind(d,fs_tmp); 
  
  %take light bolo only in these halfscans
  if isfield(d.antenna0,'bolo')
    v=double(d.antenna0.bolo.mag(mapind,ind));
  else
    v=double(d.mce0.data.fb(mapind,ind));
  end
  v0=v;

  %mask out timestream
  v(psm(mapind,ind))=NaN;
  numNaN=zeros(size(v));
  numNaN(isnan(v))=1;

  %reshape
  v=reshape(v,npoints,nhalfscans,nch);
  v0=reshape(v0,npoints,nhalfscans,nch);  
  numNaN=reshape(numNaN,npoints,nhalfscans,nch);

  %split left and right going
  l=find(fs_tmp.inc == 1);
  r=find(fs_tmp.inc == 0);
  vl=v(:,l,:);
  vr=v(:,r,:);
  vl0=v0(:,l,:);
  vr0=v0(:,r,:);
  numNaNl=numNaN(:,l,:);
  numNaNr=numNaN(:,r,:);

  %take mean
  vl_m=nanmean(vl,2);
  vr_m=nanmean(vr,2); 
  numNaNl=sum(numNaNl,2);
  numNaNr=sum(numNaNr,2);  

  %mask out ground fixed signal where > .25 of the halfscans at each azimuth were NaN
  vl_m(numNaNl>nhalfscans*.25)=NaN;
  vr_m(numNaNr>nhalfscans*.25)=NaN;
  
  %interpolate over missing regions in the ground fixed signal  
  az=reshape(d.azoff(mapind),npoints,nhalfscans);
  azl=mean(az(:,l),2);
  azr=mean(az(:,r),2);
  
  for j=1:size(vl_m,3)
      vl_tmp=vl_m(:,1,j);      
      vr_tmp=vr_m(:,1,j);    
      badindl=isnan(vl_tmp);
      badindr=isnan(vr_tmp);
      numbadl=length(find(badindl==1));
      numbadr=length(find(badindr==1));
            
    if(numbadl>0 & numbadl<npoints)
      x=azl(~badindl);
      y=vl_tmp(~badindl);
      xi=azl(badindl);
      
      % mirror sky
      x=[x-360;x;x+360];
      y=[y;y;y];

      vil=interp1(x,y,xi,'cubic');
      vl_tmp(badindl)=vil;
      vl_m(:,1,j)=vl_tmp;
    end    
    if(numbadr>0 & numbadr<npoints)  
      x=azr(~badindr);
      y=vr_tmp(~badindr);
      xi=azr(badindr);
      
      % mirror sky
      x=[x-360;x;x+360];
      y=[y;y;y];

      vir=interp1(x,y,xi,'cubic');
      vr_tmp(badindr)=vir;
      vr_m(:,1,j)=vr_tmp;
      
      %plot(x,y,'.');hold on;
      %plot(xi,vir,'.r');hold off
      %pause     
    end
    
  end
  
  %subtract mean from each halfscan
  vl0=vl0 - repmat(vl_m,1, size(l,1));
  vr0=vr0 - repmat(vr_m,1, size(r,1));

  %put it back into d
  v0(:,l,:)=vl0;
  v0(:,r,:)=vr0;
  v0=reshape(v0,nhalfscans*npoints,nch);

  if isfield(d.antenna0,'bolo')
    d.antenna0.bolo.mag(mapind,ind)=v0;
  else
    d.mce0.data.fb(mapind,ind)=v0;
  end
end


return
    
 
