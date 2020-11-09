function [v,n,t]=cut_incomp_chanblk(v,fs,ind,frac)
% [v,n,t]=cut_incomp_chanblk(v,fs,ind,frac)
%
% Set any channel/period which has greater than frac fraction of NaN's
% in it to all NaN
%
% This is called from reduc_applycal to get rid of scan sets which are
% heavily NaN'd due to ADC rail filtering - these end up with
% under-estimated noise in reduc_makesim due to filling in NaN's with
% constant number.
% Also called from reduc_makecomap we decline to add incomplete
% half-scans into the map as there's no need (they are tiny fraction
% of total data) and might cause problems.

disp('cut_incomp_chanblk');

% default to zero NaN tolerance
if(~exist('frac','var'))
  frac=0;
end

n=0;
t=0;
for i=1:length(fs.sf)
  s=fs.sf(i); e=fs.ef(i);
  l=e-s+1;
  
  for j=1:length(ind)
    t=t+1;
    
    f=sum(isnan(v(s:e,ind(j))));
    
    if(f/l>frac)
      %[i j]
      %plot(v(s:e,ind(j))); pause
      
      v(s:e,ind(j))=NaN;
      n=n+1;
    end
  end
end

return
