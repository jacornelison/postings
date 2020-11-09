function [C,I]=minmd(A)
% [minval,minind]=minmd(x)
%
% Find the smallest element of a multi dimensional array
% and it's indexes, working over all dimensions, not just
% the first singleton as does normal min(A) function.
%
% Warning - will not currently work if more than one
% identical min values

[C,I]=min(A(:));

str='[';
for i=1:ndims(A)
  str=strcat(str,sprintf('I(%d),',i));
end
str=str(1:end-1);
eval([str,']=ind2sub(size(A),I);']);

return

% Old code from before I noticed ind2sub function

[C,ind]=min(A(:));
dims=size(A);
for i=length(dims):-1:2
  p=prod(dims(1:i-1));
  [q,r]=div(ind,p);
  if(r>0)
    I(i)=q+1;
    ind=r;
  else
    I(i)=q;
    ind=p;
  end
end
I(1)=ind;
