function s=summd(s,d)
% s=summd(s,d)
%
% take the sum of s over all *except* the specified dimensions

for i=1:length(size(s))
  if(~any(i==d))
    s=sum(s,i);
  end
end

s=squeeze(s);

return
