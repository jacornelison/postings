function c=expand_cuts(c,p)
% c=expand_cuts(c,p)
%
% Expand all the cuts in structure c to be
% half-scans x channels cutmask cm
% p is used to find which channels belong to which rx

fields=fieldnames(c)';
fields=setdiff(fields,{'nhs','nch','nrx','nmce'});

% expand all fields to n-half-scan x n-channels
for f=fields
  f=f{1};

  % fetch the cut array
  x=c.(f);
  s=size(x);
  
  % if needed expand first dim to n half-scans
  switch s(1)
    case 1
      x=repmat(x,[c.nhs,1]);
  end
  
  % if needed expand second dim to n channels
  switch s(2)
    case 1
      x=repmat(x,[1,c.nch]);
    case c.nrx
      for r=unique(p.rx)'
        rp=find((p.rx)==r);
        y(:,rp)=repmat(x(:,r+1),[1,length(rp)]);
      end
      x=y;
    case c.nmce
      for r=unique(p.mce)'
        rp=find((p.mce)==r);
        y(:,rp)=repmat(x(:,r+1),[1,length(rp)]);
      end
      x=y;
  end

  % put the cut array back
  c.(f)=x;
end

return
