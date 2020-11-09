function d=clean_blips(d,p,a,ind)
% d=clean_blips(d,p,a,ind)
%
% Seems there should be a general solution for how best to clean
% source crossing scans to remove ramps and atmosphere.
% At the moment in reduc_initial for crossobs and rowcal we do
% simplest possible fit with no cleaning.
% In reduc_findtc which is basically just a special rowcal we do
% something a little more sophisticated.
% Up to now radpoint_reduc just subtracted another channel without any
% scaling (no regress).
% Notice near some of pointing sources there are other bright blobs
% which mess up trying to make a template.
% At the moment this code not used.

% for each scan
for i=1:length(a.s)
%for i=1
  
  % what is the deck angle for this scan set
  dk=round(d.tracker.actual(a.s(i),3));

  % rotate array to this deck angle
  q=rotarray(p,dk);
    
  s=a.sf(i); e=a.ef(i);
  
  clear y yp
  
  % this is the copy we will subtract from
  y=double(d.lockin.adcData(s:e,:));

  % this is the flagged copy we will punch to NaN when close to source
  % and then use to fit to and make templates from
  yf=y;
  
  subplot(1,4,1)
  for k=ind.gl
    yp(:,k)=y(:,k)+k*0.1;
  end
  plot(d.tf(s:e),yp(:,ind.gl),'.','MarkerSize',1); axis tight; 
  
  % NaN when close to source
  for j=ind.gl
    
    % find distance to source
    r=pyth(d.azoffsky(s:e)+q.ra_off_dos(j),d.eloff(s:e)-q.dec_off(j));
    
    % if source close to feed NaN the timestream
    yf(r<0.2,j)=NaN;
  end

  subplot(1,4,2)
  for k=ind.gl
    yp(:,k)=yf(:,k)+k*0.1;
  end
  plot(d.tf(s:e),yp(:,ind.gl),'.','MarkerSize',1); axis tight; 
  
  % remove linear in time from each channel
  % (helps to clean up ramping channels)
  clear X
  X(:,2)=[1:size(y,1)]'; X(:,1)=1;
  for j=ind.gl
    if(any(~isnan(yf(:,j))))
      b=regress(yf(:,j),X); % regress against flagged version
      yf(:,j)=yf(:,j)-X*b;  % subtract from flagged
      y(:,j)=y(:,j)-X*b;    % subtract from un-flagged
    end
  end
  
  subplot(1,4,3)
  for k=ind.gl
    yp(:,k)=y(:,k)+k*0.1;
  end
  plot(d.tf(s:e),yp(:,ind.gl),'.','MarkerSize',1); axis tight; 
  
  % make templates using flagged version
  clear t100 t150
  t100(:,2)=nanmean(yf(:,ind.gl100),2); t100(:,1)=1;
  t150(:,2)=nanmean(yf(:,ind.gl150),2); t150(:,1)=1;
  
  for k=1:size(yf,1)
    z=yf(k,ind.gl100);
    z=z(~isnan(z));
    t100(k,2)=mean(z);
  end
  
  % regress template from channels at that freq
  for j=ind.gl100
    if(any(~isnan(yf(:,j))))
      b=regress(yf(:,j),t100); % regress against flagged version
      y(:,j)=y(:,j)-t100*b;    % subtract from un-flagged
    end
  end
  for j=ind.gl150
    if(any(~isnan(yf(:,j))))
      b=regress(yf(:,j),t150);
      y(:,j)=y(:,j)-t150*b;
    end
  end

  subplot(1,4,4)
  for k=ind.gl
    yp(:,k)=y(:,k)+k*0.1;
  end
  plot(d.tf(s:e),yp(:,ind.gl),'.','MarkerSize',1); axis tight; 
  for k=ind.gl
    text(d.tf(s),k*0.1,q.channel_name(k),'HorizontalAlignment','right','color','b');
  end

  gtitle(sprintf('src=%s dk=%d sc=%d',d.tracker.source(a.s(i),:),dk,i));
  
  pause
  
  % put un-flagged version back in
  d.lockin.adcData(s:e,:)=y;
end

return

