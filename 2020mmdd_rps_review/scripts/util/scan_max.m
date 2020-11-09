function mv=scan_max(d,fs)
% mv=scan_max(d,fs)
%
% Find maximum value occuring in each scan

% for each scan
mv=[];

for i=1:length(fs.sf)
  s=fs.sf(i); e=fs.ef(i);
  mv(i,:)=nanmax(d.mce0.data.fb(s:e,:));
end

return
