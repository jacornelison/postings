function d=framecut(d,ind)
% d=framecut(d,ind)
%
% Cutdown data structure by frames
%
% e.g. to remove all frames except when one or more flags set:
%    d=framecut(d,d.frame.features>0);

ind=logical(ind);

[sampratio,samprate]=get_sampratio_rate(d);

% Figure out if cutting over fast or slow reg
if(size(d.ts,1)==length(ind))
  inds=ind;
  indf=cvec(repmat(ind',[sampratio,1]));
else
  indf=ind;
  x=reshape(indf',sampratio,[]);
  inds=cvec(x(1,:));
  indf=cvec(repmat(inds',[sampratio,1]));
end

d=framecuts(d,inds,indf);

return

function d=framecuts(d,inds,indf)

names=fieldnames(d);
for i=1:length(names)
  if(eval(sprintf('isstruct(d.%s)',names{i})))
    eval(sprintf('d.%s=framecuts(d.%s,inds,indf);',names{i},names{i}));
  else
    if(size(eval(['d.',names{i}]),1)==size(inds,1))
      if(~isempty(getfield(d,names{i})))
	eval(sprintf('d.%s=d.%s(inds,:);',names{i},names{i}));
      end
    else
      if(~isempty(getfield(d,names{i})))
	eval(sprintf('d.%s=d.%s(indf,:);',names{i},names{i}));
      end
    end
  end
end

return
