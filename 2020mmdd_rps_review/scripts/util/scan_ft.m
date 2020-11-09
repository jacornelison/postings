function spec=scan_ft(d,fs,nset,ind,scancheck)

disp(sprintf('scan_ft...'));

% if scancheck not explicity passed, use all half scans
if(~exist('scancheck','var'))
  scancheck=false(size(d.antenna0.bolo.mag,2),size(fs.sf,1));
end

samprate=round(1/(d.t(2)-d.t(1))/24/3600);	
nscanset=length(fs.sf)/round(nset);
nfb=12; % db 080317 defines number of freq bins here.

% for each scanset (50mn between el step)
for i=1:nscanset
  clear dat;
  z=1; 

  % for each half-scan in the set
  for j=1:nset
    
    % get the half-scan
    k=(i-1)*nset+j;		
    s=fs.sf(k); e=fs.ef(k);		
    v=double(d.antenna0.bolo.mag(s:e, ind.all));
    % only use this half-scan for pow spec estimation if no channels
    % were lost due to glitch removal
    % (a glitch in a single channel can cause resulting covariance
    % matrix to be non-posdef)
    if(~any(scancheck(:,k)))   
   
     % do the ft
      % see http://www.mathworks.com/support/tech-notes/1700/1702.html
      if(size(v,1)<=200) 
	n=200; 
      else 
	n=2^nextpow2(size(v,1)); 
      end
    
      ft=fft(v,n);
    
      
      % scale so result not a function of unpadded input length
      % Spectra will be in "resolution" independent form
      % suitable for binning up in freq and resampling at any freq resolution
      ft=ft/sqrt(size(v,1));
      
      % generate auto/cross spectra
      if(z==1)
	%allocate Nans for speed up for loop.(dB 080317)
        dat=NaN(n/2+1,size(v,2),size(v,2),nset);
      end

      for l=1:size(v,2)
        for m=l:size(v,2)	
          % make spectrum	 
          s=ft(:,l).*conj(ft(:,m));
          dat(:,l,m,z)=2*s(1:n/2+1);
        end
      end
      
      z=z+1;
    end
  end

  disp(strcat(num2str(i),'=i, half-scans used=',num2str(z-1)));

  if(z==1)
    disp('No half-scans without glitch in at least one channel for this half scan. Skipping noisemodel for this phase.')

spec.ps=zeros(nscanset,size(v,2),size(v,2),nfb);
spec.nhit=zeros(nscanset,nfb);
spec.bs=zeros(nscanset, nfb);
spec.be=zeros(nscanset,nfb+1);
    return
  end
 
  % calc the freq axis	
  f=linspace(0,samprate/2,size(dat,1))'; % 0 to Nyquist for sample rate
  
  % throw out zero freq comp
  dat=dat(2:end,:,:,:);
  f=f(2:end);

  % log freq
  lf=log10(f);
      
  % expand freq to size required for hprof
  x=repmat(lf,[1,size(dat,4)]);	
  
  if(i==1)
    %allocate mu and nu
    mu=zeros(nscanset,size(v,2),size(v,2),nfb);
    nu=mu;
    sig=mu;
    nhit=zeros(nscanset,nfb);
  end
  
  % Bin up spectra over the half-scan set
  % for each channel pair bin up using hprof (profile histogram)
  for l=1:size(v,2)
    for m=l:size(v,2)
      [bincen(i,:),mu(i,l,m,:),sig(i,l,m,:),dummy,dummy,nhit(i,:)] = ...
	  hprof(x,real(dat(:,l,m,:)),nfb,lf(1)-1e-9,lf(end)+1e-9);
      [bincen(i,:),nu(i,l,m,:),sig(i,l,m,:),dummy,dummy,nhit(i,:)] = ...
	  hprof(x,imag(dat(:,l,m,:)),nfb,lf(1)-1e-9,lf(end)+1e-9);
      binedge(i,:)=linspace(lf(1)-1e-9,lf(end)+1e-9,nfb+1);
    end
  end
  
  clear dummy sig;
  
  % There is slight wrinkle when binning one sided psd like this -
  % zero and Nyquist freq actually have "half the normalization" in
  % one sided psd - so in principle cannot be binned with other freq.
  % Here we throw out zero freq comp, and Nyquist is binned with many
  % other freq into bin which is already supressed by low pass filter
  % so effect is minimal. Could get rid of the problem by dividing psd
  % 2:n/2 by 2, which is actually what we want when come to apply the
  % spectra, but for the moment keep psd normalization conventional
  % and ignore the slight error...

  if(0)
  % plot 1,1 spectra to illustrate
  clf
  set(gcf,'DefaultAxesColorOrder',jet(nset));
  loglog(10.^x,squeeze(dat(:,34,56,:)),'.');
  be=binedge(i,:);
  for j=1:length(be)-1
    line(10.^be(j:j+1),[mu(i,34,56,j),mu(i,34,56,j)],'LineWidth',5);
  end
  title(sprintf('mean channel 1-1 spectra as derived from scan set %d',i));
  pause
  end
end

% keep in power units
%spec.ps=mu;				% mu(i,l,m, bin average)
spec.ps=complex(mu,nu)/f(1);		% V^2/Hz

% go to Hz
spec.bs=10.^bincen;
spec.be=10.^binedge;
spec.nhit=nhit;

return
