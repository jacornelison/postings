function d=deconv_scans(d,p,mask,el)
% d=deconv_scans(d,p,mask,el)
%
% Deconvolve bolometer/filter transfer functions
%
% this version does deconv directly in the time domain by convolving
% the data with the deconvolution kernel generated by get_xferfunc.  
%
% the deconvolution kernel includes a low pass filter
%
% A few issues:
%
% arg p is from get_array_info, specifies which channels are from which
% receiver
%
% optional arg mask specifies that start/end regions which are invalid
% should be masked out as NaN
%
% optional arg el specifies additional elements of d structure which
% are to be low pass filtered, apply the same filter convolved with
% deconv kernel, given in the tf.lpf structure

if(~exist('mask','var'))
  mask=[];
end

if(~exist('el','var'))
  el=[];
end

if(isempty(mask))
  mask=0;
end

tf=d.tf;

disp('deconv_scans...')

% for each data period divide by transfer function
% transfer function could be different for each period
% the different periods are stored in tf.dc
% since field scans and load curves taken with different settings
for j=1:size(tf,2)

  % loop over number of MCEs as defined in get_xferfunc.m
  for imce=1:size(tf,1)

    % get the deconv blocks from tf.dc
    s=tf(imce,j).dc.sf; e=tf(imce,j).dc.ef;

    % fetch the data
    dat=d.mce0.data.fb(s:e,1+p.mce==imce);

    % get the time domain deconv kernel, tf.deconvkern  
    % if it has 1 column, then no empirical xfer functions were used
    dct=tf(imce,j).deconvkern;

    % convolve the data with the deconv kernel
    % do this in the time domain to keep the invalid regions
    % truly zero outside the kernel length
    datdc=convn(dat,dct,'same');

    if(mask==1) && ~isempty(dct)
      % set the invalid regions to NaN
      % since deconv in the time domain with FIR kernel
      % cut invalid regions based on central sample
      [a,centsamp]=max(dct);

      % Mask width is half-width of deconvolution kernel + half-width of GCP filter post-downsample
      nmask1=centsamp+ceil(length(tf(imce,j).gcp_taps)/tf(imce,j).gcpds/2);
      datdc(1:nmask1,:)=NaN;
      nmask2=(length(dct)-centsamp)+ceil(length(tf(imce,j).gcp_taps)/tf(imce,j).gcpds/2);
      datdc(end-nmask2+1:end,:)=NaN;

      % Store masked sample information
      dcmask=[];
      dcmask.sf=s-1+[1;size(datdc,1)-nmask2+1];
      dcmask.ef=s-1+[nmask1;size(datdc,1)];
      d.tf(imce,j).dcmask=dcmask;
    end

    % store result
    d.mce0.data.fb(s:e,1+p.mce==imce)=datdc;
    
    % plot
    if(0)
      clf
      setwinsize(gcf,1000,800);
      subplot(2,1,1)
      x=[0:length(dat)-1]*1/samprate;
      plot(x,dat);
      hold on; plot(x,datdc,'r'); hold off
      axis tight
      xl=xlim; xlim([xl(1)-1,xl(2)+1]);
      %yl=ylim; yl(yl>10)=10; yl(yl<-10)=-10; ylim(yl);
      grid; xlabel('time (sec)');
      title(sprintf('deconv channel gcp%d',p.gcp(i)));
      %legend({'input','output'},'Location','SouthWest');
      legend({'input','output'},'Location','Best');
      
      % generate t domain deconv kernel which applies the deconv/lowpass
      dc=fftshift(ifft(ff./impf,'symmetric'));
      tt=fftshift(t); 
      if(0)
        % do t domain deconv - this produces exactly the same result
        datdct=convn(dat,dc(abs(tt)<15),'same');
        hold on; plot(x,datdct,'g'); hold off
      end
      
      subplot(2,2,3); plot(tt,dc); xlim([-0.5,1.5]); grid
      xlabel('time (sec)');
      title(sprintf('deconv kernel'))
      
      subplot(2,2,4); plot(tt,dc);
      plot(x,dat);
      hold on; plot(x,datdc,'r'); hold off
      axis tight
      zoom xon; zoom(10);
      grid; xlabel('time (sec)');
      
      pause
      %drawnow
      
      if(0)
        mkgif(sprintf('deconv/%02d.gif',i));
        fh=fopen(sprintf('deconv/%02d.html',i),'w');
        fprintf(fh,'<a href="../index.html">up</a> ------ ');
        if(i==ind(1))
          z=ind(end);
        else
          z=ind(k-1);
        end
        fprintf(fh,'<a href="%02d.html">last</a> ',z);
        if(i==ind(end))
          z=ind(1);
        else
          z=ind(k+1);
        end
        fprintf(fh,'<a href="%02d.html">next</a> ',z);
        fprintf(fh,'<p><img src="%02d.gif"> ',i);
        fclose(fh);
      end
    end
  end
  
  % also low pass filter additional structure fields basically pointing info
  % apply the same low pass filter convolved with deconv kernel in get_xferfunc

  % make sure we have a good low-pass
  lpft=[];
  % loop over number of MCEs as defined in get_xferfunc.m
  for imce=1:size(tf,1)
    if isfield(tf(imce,j),'lpf')
      lpft=tf(imce,j).lpf;
    end
    if ~isempty(lpft)
      break;
    end
  end

  % If lpft is empty, that means get_xferfunc didn't define a low pass filter,
  % so don't do anything here
  if ~isempty(lpft)
    for ii=1:length(el)
      eval(sprintf('dat=d.%s(s:e,:);',el{ii}));   
      datdc=convn(dat,lpft,'same');
      eval(sprintf('d.%s(s:e,:)=datdc;',el{ii}));
    end
  end
end

return
