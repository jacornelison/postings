function ptitle(t,x,y,varargin)
% ptitle(t,x,y,varargin)
%
% simple function to put a label inside the axes
% varargin passed to ML's text()

if(~exist('x','var') || isempty(x))
  x=0.03;
end
if(~exist('y','var') || isempty(y))
  y=0.9;
end

xl=xlim; yl=ylim;

xp=xl(1)+x*diff(xl); yp=yl(1)+y*diff(yl);

text(xp,yp,t,varargin{:});

return
