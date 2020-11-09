function errorbox(x,y,e,c)
% errorbox(x,y,e,c)
%
% Make an errorbar plot with shaded rectangles instead of bars
%
% e.g.
% x = 1:10;
% y = sin(x);
% e = ones(size(x));
% errorbox(x,y,e)

if(~exist('c','var'))
  c='b';
end

x=rvec(x);
y=rvec(y);
e=rvec(e);

dx=(x(2)-x(1))/2;
xg=[x-dx;x-dx;x+dx;x+dx];
yg=[y-e;y+e;y+e;y-e];

patch(xg,yg,c);

box on;

return
