function sout=structcat(dim,sin)
% s=structcat(dim,sin)
%
% Concatenate each field of a set of arrays
% Do sout.field=cat(dim,sin.field) for each field
%
% Modified so will cope with structures of structures
%
% eg: s=structcat(1,[s1,s2]);

% Keep consistent with earlier version which lacked dim arg
if(isstruct(dim))
  sin=[dim,sin];
  dim=2;
end

names=fieldnames(sin);

sout=[];

for i=1:length(names)
  if(isstruct(sin(1).(names{i})))
    x=cat(dim,sin.(names{i}));
    x=structcat(dim,x);
  else
    x=cat(dim,sin.(names{i}));
  end
  sout.(names{i}) = x;
end

return
