function acc=multac(ac,m)
% function acc=multac(ac,m)
%
% multiply an ac structure - as if we had m times more data

for i=1:numel(ac)

  fnames=fieldnames(ac(i));

  for f=1:length(fnames)
    if(strfind(fnames{f},'max'))
      acc(i).(fnames{f})=ac(i).(fnames{f});
    else
      if ~iscell(ac(i).(fnames{f}))
        acc(i).(fnames{f})=ac(i).(fnames{f})*m(i);
      else
        for k=1:numel(ac(i).(fnames{f}))
          acc(i).(fnames{f}){k}=ac(i).(fnames{f}){k}*m(i);
        end
      end
    end
  end
end
  
acc=reshape(acc,size(ac));

return
