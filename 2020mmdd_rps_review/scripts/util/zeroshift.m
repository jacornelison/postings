function x=zeroshift(x,n,val)
% ZEROSHIFT Shift array, pad missing values
% B = ZEROSHIFT(A,SHIFTSIZE,PADVALUE)
%
% Behaves like circshift, but instead of wrapping the array, pad values with PADVAL
% (default zero).

% First, apply circshift
x=circshift(x,n);

% Now replace wrapped values with pad value
if ~exist('val','var') || isempty(val)
  val=0;
end

ndim=numel(size(x));

for k=1:ndim

  sz=size(x,k);
  
  if n(k)>0
    s=1;
    e=min(n(k),sz);
  elseif n(k)<0
    s=max(sz+n(k)+1,1);
    e=sz;
  else
    continue
  end

  % This is super kludgy, but I can't think of a nice way to do this sort of array
  % indexing without a priori knowing the dimensionality of the array, so here it is.
  y='';
  for j=1:ndim
    if j==k
      y=[y sprintf('%d:%d',s,e)];
    else
      y=[y ':'];
    end
    
    if j<ndim
      y=[y ','];
    end
  end
  cmd=sprintf('x(%s)=val;',y);
  eval(cmd);
  
end

return


