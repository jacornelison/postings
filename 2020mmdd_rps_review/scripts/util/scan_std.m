function sd=scan_std(d,fs)
% sd=scan_std(d,fs)
%
% Find standard deviation of scans

% for each scan
sd=[];

for i=1:length(fs.sf)
  s=fs.sf(i); e=fs.ef(i);
  sd(i,:)=nanstd(d.mce0.data.fb(s:e,:));
end

return