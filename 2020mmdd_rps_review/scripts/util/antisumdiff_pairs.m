function d=antisumdiff_pairs(d,fs,ind)
% d=antisumdiff_pairs(d,fs,ind)
%
% Go from sum/diff to A/B 
disp('antisumdiff_pairs...');
% for each half scan
for i=1:length(fs.sf)
  s=fs.sf(i); e=fs.ef(i);
  v=d(s:e,:);
  
  vp=v;
  
  vp(:,ind)  =(v(:,ind)+v(:,ind+1)); % A+B
  vp(:,ind+1)=(v(:,ind)-v(:,ind+1)); % A-B

  d(s:e,:)=vp;
end

% do not touch p array because the sum/diff was only in the noise tod
% generation and there was nothing in p that reflected that, so leave p alone.


return