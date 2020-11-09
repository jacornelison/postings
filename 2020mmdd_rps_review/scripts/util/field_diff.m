function [d,fs,psm]=field_diff(d,fs,psm)
% [d,fs,psm]=field_diff(d,fs,psm)
%
% difference fields, store the result in the lead field
%
% result is normalized (x-y)/sqrt(2) to leave CMB spectrum unchanged
%
% Optional arg psm is logical point source mask - if present this
% or'd across lead/trail

% find number of half scans between lead/trail
% This is somewhat kludgy but better than before where this number was
% hard coded in get_scan_info and passed to field_diff
src=unique(d.tracker.source(fs.s,:),'rows');
src_hit1=strmatch(src(1,:),d.tracker.source(fs.s,:),'exact');
src_hit2=strmatch(src(2,:),d.tracker.source(fs.s,:),'exact');
fs1=structcut(fs,src_hit1);
fs2=structcut(fs,src_hit2);
x=fs1.s<fs2.s(1);
x=find(x==1);
n=x(end);

ind=false(size(fs.s));
for i=0:n*2:length(fs.s)-1
  for j=1:n
    s1=fs.sf(i+j);   e1=fs.ef(i+j);
    s2=fs.sf(i+j+n); e2=fs.ef(i+j+n);
    d.lockin.adcData(s1:e1,:)=(d.lockin.adcData(s1:e1,:)-d.lockin.adcData(s2:e2,:))/sqrt(2);
    if(exist('psm','var'))
      psm(s1:e1,:)=psm(s1:e1,:)|psm(s2:e2,:);
    end
    ind(i+j)=true;
  end
end

fs=structcut(fs,ind);

return
