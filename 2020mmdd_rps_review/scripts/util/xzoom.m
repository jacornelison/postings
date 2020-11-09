function xzoom(xl)
% xzoom(n,m)
%
% Interactively zoom x on a multi panel plot

if(~exist('xl','var'))
  xl='auto';
end

% get vector of handles to all axes in current figure
c=get(gcf,'Children');

% loop forever waiting for user mouse clicks
while 1
  [x,y,b] = ginput(2);
  if(b(2)==3)
    % for each axes
    for i=1:length(c)
      xlim(c(i),xl);
    end
  else
    for i=1:length(c)
      xlim(c(i),x);
    end
  end
end

return
