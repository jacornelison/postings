function c=compress_perhs(c,opt)
% c=compress_perhs(c,opt)
%
% Compress structure of cut parameters or cut masks over the
% half-scans (first) dimension

if ~exist('opt','var')
  opt='mean';
end

fn=fieldnames(c);
for i=1:length(fn)
  if ismember(fn{i},{'tag','nhs','nch','nrx','nmce'})
    continue
  end
  tmp=c.(fn{i});
  if size(tmp,1)<c.nhs
    continue
  end
  switch(lower(opt))
    case 'max', tmp=max(tmp,[],1);
    case 'min', tmp=min(tmp,[],1);
    case 'mean', tmp=mean(tmp,1);
  end
  c.(fn{i})=tmp;
end

return
