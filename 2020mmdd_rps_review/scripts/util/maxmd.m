function [C,I]=maxmd(A)
% [maxval,maxind]=maxmd(x)
%
% Find the largest element of a multi dimensional array
% and it's indexes, working over all dimensions, not just
% the first singleton as does normal min(A) function.
%
% Warning - will not currently work if more than one
% identical max values

[C,I]=max(A(:));

str='[';
for i=1:ndims(A)
  str=strcat(str,sprintf('I(%d),',i));
end
str=str(1:end-1);
eval([str,']=ind2sub(size(A),I);']);

return
