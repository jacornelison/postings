function [freq,spec,nrep]=avpsd(input, fsample, fres, med, deg)
%function [freq,spec,nrep]=avpsd(input, fsample, fres, med, deg)
%make an average psd by averaging the psd from small segments of resolution fres
%pb if the data is not continuous on long timescales
% if dealing  with data cut up in halfscans, use avpsd2.m
  
if(~exist('med','var'))
  med=[];
end
if isempty(med)
  med=0;
end
                    
if(~exist('deg','var'))
  deg=[];
end
if isempty(deg)
  deg=-1;
end
 
n = length(input);
nout = floor(fsample / fres);
nrep = floor(n / nout);
x=cvec(linspace(1,nout,nout));

psdarr = zeros(nrep, nout / 2 -1);

if (deg ~=-1)
    for i = 1:nrep 
      y=input(1 + (i-1)*nout : i * nout );
      p=polyfit(x,y, deg);
      baseline=polyval(p,x);
      y=y-baseline;
      [freq,ps]=psd(y, fsample);
      psdarr(i, :) = ps;
    end
else
    for i = 1:nrep
       y=input(1 + (i-1)*nout : i * nout );
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

