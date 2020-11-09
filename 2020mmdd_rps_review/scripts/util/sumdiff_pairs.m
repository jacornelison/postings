function [d,p]=sumdiff_pairs(d,p,fs,inda,indb)
% [d,p]=sumdiff_pairs(d,p,fs,inda,indb)
%
% Go from A/B to sum/diff
%
% e.g:
% [d,p]=sumdiff_pairs(d,p,fs,ind.rgla,ind.rglb)

disp('sumdiff_pairs...');

% for each period
for i=1:length(fs.sf)
  s=fs.sf(i); e=fs.ef(i);
  v=d.mce0.data.fb(s:e,:);
  
  vp=v;
  
  vp(:,inda)=(v(:,inda)+v(:,indb))/2; % A+B in A channels
  vp(:,indb)=(v(:,inda)-v(:,indb))/2; % A-B in B channels

  d.mce0.data.fb(s:e,:)=vp;
end

% massage p array to reflect sum/diff
for i=1:length(inda)
  a=inda(i); b=indb(i);
  
  [r,theta]=average_pairs(p.r([a,b]),p.theta([a,b]));
  p.r([a,b])=r; p.theta([a,b])=theta;
  
  p.pol{a}='sum'; p.pol{b}='diff';
end

return