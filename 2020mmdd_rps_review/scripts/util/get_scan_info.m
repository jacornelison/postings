function fs=get_scan_info(scanname,fs)
% fs=get_scan_info(scanname,fs)
%
% scan name is stored in archive but constant velocity scan throw and
% rate are not - therefore must get from an external file.
%
% Note that when use a new scan profile must add it to the external
% file - otherwise this function returns NaN for scan throw and rate.

% read in scan profile info
fid=fopen('aux_data/scansource.txt');
scans=textscan(fid,'%s %n %n','commentStyle','#');
fclose(fid);

fields={'name', 'throw', 'rate'};
scans=cell2struct(scans, fields, 2);
scans.throw=scans.throw*2;

for i=1:size(scanname,1)
  j=strmatch(deblank(scanname(i,:)),scans.name,'exact');
  if(~isempty(j))
    fs.throw(i,1)=scans.throw(j);
    fs.rate(i,1)=scans.rate(j);
  else
    fs.throw(i,1)=NaN;
    fs.rate(i,1)=NaN;
  end
end

return
