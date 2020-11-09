function s=structcut(s,x)
% structure=structcut(structure,i)
%
% I find I often want to shrink a structure
% of the form
%
%   s.a(n)
%   s.b(n)
%   s.c(n)
% 
% according to some matrix of indices i.
% For example it's a pain to do:
%
%   x=s.a>1;
%   s.a=s.a(x);
%   s.b=s.b(x);
%   s.c=s.c(x);
%
% and if I add another field I have to change
% the code in multiple places.
%
% This function automatically cuts all fields
% so the code becomes:
%
%   x=s.a>1;
%   s=structcut(s,x);
%
% and continues to work if fields are added.
%
% Now extended so that fields can be 2D - will try to figure
% out which dimension to cut over by looking for the dim which is
% equal length for all fields. If all are same cut will occur over
% rows.
%
% Extended again so fields can be up to 3D - cut dim must be 1st or 2nd

names=fieldnames(s);

% Get sizes of fields
for i=1:length(names)
  d(i,1)=size(getfield(s,char(names(i))),1);
  d(i,2)=size(getfield(s,char(names(i))),2);
end

if(all(d(:,1)==d(1,1)) & d(1,1)~=1) ...
  || (all(d(:,1)==1) & all(d(:,2)==1))
  for i=1:length(names)
    d=getfield(s,char(names(i)));
    switch ndims(d)
      case 2
	d=d(x,:);
      case 3
	d=d(x,:,:);
    end
    s=setfield(s,char(names(i)),d);
  end
elseif(all(d(:,2)==d(1,2)) & d(1,2)~=1)
  for i=1:length(names)
    d=getfield(s,char(names(i)));
    switch ndims(d)
      case 2
        d=d(:,x);
      case 3
        d=d(:,x,:);
    end      
    s=setfield(s,char(names(i)),d);
  end
else
  % do nothing
end
