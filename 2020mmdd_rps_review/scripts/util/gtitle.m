function h=gtitle(str,y,interpreter)
% function gtitle(str,y_placement,interpreter)
%
% Put a global title on the current figure
%
% str = string
% y_placement = normal coordinate vertical placement of text (default 0.96)
% interpreter = 'tex' (default LaTeX characters)
%             = 'none' (literal characters)

if(~exist('y','var'))
  y=[];
end
if(isempty(y))
  y=0.96;
end

if(~exist('interpreter','var'))
  interpreter='tex';
end

%ca=gca;
axes('Position',[0,0,1,1],'Visible','off')
h=text(0.5,y,str,'Units','Normalized','FontSize',12,...
    'HorizontalAlignment','center','Interpreter',interpreter);
%axes(ca);

return
