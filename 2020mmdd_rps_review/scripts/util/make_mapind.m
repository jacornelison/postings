function mapind=make_mapind(d,fs)
% mapind=make_mapind(d,fs)
%
% Build index array pointing from scan start/stop pointers

mapind=false(length(d.t),1);
for i=1:length(fs.sf)
  s=fs.sf(i); e=fs.ef(i);
  mapind(s:e)=true;
end

return
