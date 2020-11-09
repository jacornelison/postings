function v=exp_scan_val(d,fs,val,samenan)
% v=exp_scan_val(d,fs,v,samenan)
%
% expand per scan value array val so that explicitly
% lists values across each scan

if ~exist('samenan','var') || isempty(samenan)
  samenan = true;
end

% make output same type as input
v=zeros(length(d.t),size(val,2),class(val));
    
% expand val to cover each scan
for i=1:length(fs.sf)
  s=fs.sf(i); e=fs.ef(i);
  v(s:e,:)=repmat(val(i,:),[e-s+1,1]);
end

if samenan
  % where d.mce0.data.fb is NaN, v should be NaN also
  v(isnan(d.mce0.data.fb))=NaN;
end

return
