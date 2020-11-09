function [d,coeff]=filter_scans(d,fs,filttype,ind,psm,polytype)
% d=filter_scans(d,fs,filttype,ind,psm)
%
% Filter half scans in the specified manner
%
% Optional arg psm is point source mask
% For polynomial sub fit is made to masked data and then subtracted
% from unmasked. This gets rid of the stripes either side of a bright
% point source.
%
% optional argument type determines whether to do a normal filtering
% or a filtering with the edge point pegged ( most useful with psm
% masking)
% polytype='pegged' only works in polynomial filtering for now
%
% e.g. to filter sum/diff with different poly orders after
% sumdiff_pairs do:
% d=filter_scans(d,fs,'p3',ind.gla)
% d=filter_scans(d,fs,'p0',ind.glb)

if(~exist('ind','var'))
  ind=1:size(d.mce0.data.fb,2);
end

if(~exist('psm','var'))
  psm=false(size(d.mce0.data.fb));
end

if(~exist('polytype','var'))
  polytype=[];
end
if(isempty(polytype))
  polytype='normal';
end

disp(sprintf('filter_scans...with type %s',polytype));

% suppress some warnings from an mldivide in regress_no_arg_check
warning('off','MATLAB:nearlySingularMatrix')

coeff = [];
switch filttype(1)
  case 'n'
    % do nothing
    
  case {'p', 'd', 's'}	% 'd' for BICEP dif/sum
    % subtract poly of required order
    [d,coeff]=polysub_scans(d,fs,str2num(filttype(2:end)),ind,psm,polytype);
    
  case 'q'
    % subtract poly using point into scan as regressor rather than
    % azoff
    d=polysub_scans_pnt(d,fs,str2num(filttype(2:end)),ind,psm);
    
  case 'b'
    % apply Butterworth high pass filter
    d=butfilt_scans(d,fs,str2num(filttype(2:end)),ind,psm);
    
  case 'm'
    % median subtract
    d=mediansub_scans(d,fs,ind);

  case 't'
    % percentile filtering
    d=percentilesub_scans(d,fs,str2num(filttype(2:end)),ind);
    
  case 'c'
    % remove chebyshev polynomial modes
    [d,coeff]=chebpolysub_scans_pnt(d,fs,str2num(filttype(2:end)),ind,psm);

end

warning('on','MATLAB:nearlySingularMatrix')

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function d=mediansub_scans(d,fs,ind)

for i=1:length(fs.sf)
  s=fs.sf(i); e=fs.ef(i);
  
  v=double(d.mce0.data.fb(s:e,ind));

  v=bsxfun(@minus,v,nanmedian(v,1));
  
  d.mce0.data.fb(s:e,ind)=v;
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function d=butfilt_scans(d,fs,cutfreq,ind,psm)

samprate=10;
[b,a]=butter(3,cutfreq/(samprate/2),'low');

scanset=unique(fs.set);

% for each scanset
for i=1:length(scanset)
  
  q=find(fs.set == scanset(i));
  fscut=structcut(fs,q);
  mapind=make_mapind(d,fscut);
  v=double(d.mce0.data.fb(mapind,ind));
  azz=d.azoff(mapind);
  x=linspace(1,size(v,1),size(v,1));
  m=psm(mapind,ind);
 
  % copy the  signal
  vp=v;
  
  for j=1:size(v,2)
  %remove the galactic signal by interpolating through the gal points
  % we could also remove a signal from WMAP resampled
  q=find(m(:,j) == 0);
  qp=find(m(:,j) == 1);
  
  vp(:,j)=interp1(x(q),vp(q,j),x,'linear',mean(vp(q,j)));
  
  if(j == 0)
  plot(v(:,j));hold on
  plot(x(q),vp(q,j),'r') ;pause
  end
  end
   
  % apply the highpass filter
  vp=filtfilt(b,a,vp);
  v=v-vp;
 
  d.mce0.data.fb(mapind,ind)=v;
end

return
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [d,pcoeff]=polysub_scans(d,fs,porder,ind,psm,polytype)

pcoeff=zeros(length(fs.sf), length(ind), porder+1);

% for each half scan
for i=1:length(fs.sf)
  
  s=fs.sf(i); e=fs.ef(i);
  azz=d.azoff(s:e);
  
  % 060130: need double typecast or precision not sufficient when
  % scans have DC offset added in
  v=double(d.mce0.data.fb(s:e,ind));
   % take copy of the data and mask out point sources in the copy
  vp=v;
  vp(psm(s:e,ind))=NaN;
 
   
  % When cross obs disagree we can't make pointing correction and
  % azoff becomes NaN - this makes regress below fail.
  % Ignoring such scans here is OK as they cannot contribute to maps
  % since their pointing info is NaN.
  if(~any(isnan(azz)))
  
    % contruct the regressor
    clear X
    for j=0:porder
      X(:,j+1)=azz.^j;
    end
  
    % loop over all channels
    for j=1:size(v,2)
      if(any(~isnan(v(:,j))))
        nscan=length(vp(:,j));
        nnan=length(find(isnan(vp(:,j))));
        
      	% regress against the point source masked data
        switch polytype
         case 'normal'
          b=regress_no_arg_check(vp(:,j),X);
          baseline=X*b;

         % Use softened extrapolation in hanging PSM region
         case 'psmsoft'
          b=regress_no_arg_check(vp(:,j),X);
          baseline=X*b;
          baseline_orig=baseline;
          if isnan(vp(1,j))
            k0=1; k1=find(~isnan(vp(:,j)),1,'first');
            c=false(size(vp(:,j)));
            c(k0:k1)=true;
            az0=azz(k0); az1=azz(k1);
            if az0~=az1
              if porder<1
                bsoft(3)=0;
              else
                bsoft(3)=1/2/(az1-az0) * polyval((porder:-1:1)'.*(b(end:-1:2)),az1);
              end
              bsoft(2)=-2*bsoft(3)*az0;
              bsoft(1)=polyval(b(end:-1:1),az1)-bsoft(2)*az1-bsoft(3)*az1.^2;
              baseline_soft=bsoft(3)*azz.^2+bsoft(2)*azz+bsoft(1);
              baseline(c)=baseline_soft(c);
            end
          end
          if isnan(vp(end,j))
            k0=length(vp(:,j)); k1=find(~isnan(vp(:,j)),1,'last');
            c=false(size(vp(:,j)));
            c(k1:k0)=true;
            az0=azz(k0); az1=azz(k1);
            if az0~=az1
              if porder<1
                bsoft(3)=0;
              else
                bsoft(3)=1/2/(az1-az0) * polyval((porder:-1:1)'.*(b(end:-1:2)),az1);
              end
              bsoft(2)=-2*bsoft(3)*az0;
              bsoft(1)=polyval(b(end:-1:1),az1)-bsoft(2)*az1-bsoft(3)*az1.^2;
              baseline_soft=bsoft(3)*azz.^2+bsoft(2)*azz+bsoft(1);
              baseline(c)=baseline_soft(c);
            end
          end
          
         case 'pegged'
          n=~isnan(vp(:,j));
          end1=mean(v(1:10,j)) ;end2=mean(v(end-10:end,j));
          c=mmpolyfit(azz(n),vp(n,j),porder,'Point',[[azz(1),end1];[azz(end),end2]]); 
          baseline=c(4)+c(3)*azz+c(2)*azz.^2+c(1)*azz.^3;
        end
        
       if(0)
        if (nnan >= nscan-porder);
          %sprintf('%d  %d',nscan, nnan)
          %subplot_stack(10,11,mod(j-1,110)+1);       
          plot(azz,v(:,j),'c');title(num2str(j));
          hold on;
          plot(azz,vp(:,j),'b');
          vp=v;
          ncuts1(j)=1;
          plot(azz,baseline,'r'); 
          plot(azz(1),end1,'*r'); plot(azz(end),end2,'*r');
          set(gca,'XTickLabel',{''});
          set(gca,'YTickLabel',{''});
          hold off;
        end
       end
        % subtract from unmasked data
        v(:,j)=v(:,j)-baseline;
   
	% store the poly coefficients
	pcoeff(i,j,:)=b;
      end
    end  % end loop over channels
     
  end

  % lscov can fit multiple columns in single shot
  % but annoyingly does not treat NaNs as missing values like regress
  % checked that one chan at a time regress above gives numerically
  % identical results
  %b=lscov(X,v);
  %v=v-X*b;

  % put poly sub data back in mce0.data.fb
  d.mce0.data.fb(s:e,ind)=v;
  
end % end loop over half scans

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function d=percentilesub_scans(d,fs,pct,ind)

for i=1:length(fs.sf)
  s=fs.sf(i); e=fs.ef(i);
  
  v=double(d.mce0.data.fb(s:e,ind));
  
  v=bsxfun(@minus,v,prctile(v,pct,1));
    
  d.mce0.data.fb(s:e,ind)=v;
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function d=polysub_scans_pnt(d,fs,porder,ind,psm)
%
% Use the point into scan as regressor rather than azoff to see if
% this gets rid of weird ghost of hitmap effects in scan dir jackknife

% take copy of the data
adcData2=d.mce0.data.fb;
% mask out point sources in the copy
adcData2(psm)=NaN;

% for each half scan
for i=1:length(fs.sf)
  
  s=fs.sf(i); e=fs.ef(i);
  
  % 060130: need double typecast or precision not sufficient when
  % scans have DC offset added in
  v=double(d.mce0.data.fb(s:e,ind));
  vp=double(adcData2(s:e,ind));

  % When cross obs disagree we can't make pointing correction and
  % azoff becomes NaN - this makes regress below fail.
  % Ignoring such scans here is OK as they cannot contribute to maps
  % since their pointing info is NaN.
  if(~any(isnan(d.azoff(s:e))))
  
    % contruct the regressor
    clear X
    for j=0:porder
      X(:,j+1)=[1:length(vp)].^j;
    end
    
    % loop over all channels
    for j=1:size(v,2)
      if(any(~isnan(v(:,j))))
	% regress against the point source masked data
	b=regress_no_arg_check(vp(:,j),X);
	% subtract from unmasked data
	v(:,j)=v(:,j)-X*b;
      end
    end
  end

  d.mce0.data.fb(s:e,ind)=v;
  
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [d,pcoeff]=chebpolysub_scans_pnt(d,fs,porder,ind,psm)
%
% Use orthogonal Chebyshev polynomials

% take copy of the data
adcData2=d.mce0.data.fb;
% mask out point sources in the copy
adcData2(psm)=NaN;

X=[];
pcoeff=zeros([length(fs.sf),length(ind),porder+1]);

% for each half scan
for i=1:length(fs.sf)
  
  s=fs.sf(i); e=fs.ef(i);
  n=e-s+1;
  
  % 060130: need double typecast or precision not sufficient when
  % scans have DC offset added in
  v=double(d.mce0.data.fb(s:e,ind));
  vp=double(adcData2(s:e,ind));

  % When cross obs disagree we can't make pointing correction and
  % azoff becomes NaN - this makes regress below fail.
  % Ignoring such scans here is OK as they cannot contribute to maps
  % since their pointing info is NaN.
  if(~any(isnan(d.azoff(s:e))))
  
    % contruct the regressor if necessary
    if(size(X,1)~=n)
      clear X
      x=linspace(-1,1,n);
      for j=0:porder
	X(:,j+1)=chebpoly(j,x);
      end
    end
    
    % loop over all channels
    for j=1:size(v,2)
      if(any(~isnan(v(:,j))))
	% regress against the point source masked data
	b=regress_no_arg_check(vp(:,j),X);
	% subtract from unmasked data
	v(:,j)=v(:,j)-X*b;
	
	% store the poly coefficients
	pcoeff(i,j,:)=b;
      end
    end
    
  end

  d.mce0.data.fb(s:e,ind)=v;
  
end

return

