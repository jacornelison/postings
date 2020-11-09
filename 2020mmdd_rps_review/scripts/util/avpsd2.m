function [freq,spec,nrep]=avpsd2(input, fsample, fs, med, deg)
%function [freq,spec,nrep]=avpsd(input, fsample, fs, med, deg)
%make an average psd by averaging the psd from small segments. The resolution is 
% set here by 1 halfscan which is the longest continuous data.

if(~exist('med'))
  med=[];
end
if isempty(med)
  med=0;
end
                    
if(~exist('deg'))
  deg=[];
end
if isempty(deg)
  deg=-1;
end
 
n = length(input);
nout = size(fs.sf(20):fs.ef(20),2);
nrep = size(fs.set,1);
x=cvec(linspace(1,nout,nout));

psdarr = zeros(nrep, nout / 2 +1);

if (deg ~=-1)
    for i = 1:nrep 
      y=input(fs.sf(i):fs.ef(i));
      p=polyfit(x,y, deg);
      baseline=polyval(p,x);
      y=y-baseline;
      [freq,ps]=psd(y, fsample);
      psdarr(i, :) = ps;
    end
else
    for i = 1:nrep
      y=input(fs.sf(i):fs.ef(i));
      [freq,ps]=psd(y, fsample);
      psdarr(i, :) = ps;
    end
end

if(med)
  spec = median(psdarr,1);
else
  spec=sqrt(mean(psdarr.^2,1));
end


return


end

