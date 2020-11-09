function y = egauss2(p,x1,x2)
% y = egauss2(p,x1,x2)
%
% Compute values on an elliptical Gaussian shaped surface
% with parameters p.
% 
% p 1 = peak height
% p 2 = center in the direction of x1
% p 3 = center in the direction of x2
% p 4 = sigma major
% p 5 = sigma minor
% p 6 = rot angle of major axis ccw from x1 axis (in radians)
%
% p 7 = optional zero offset

if(length(p)<7)
  p(7)=0;
end

% center on origin
x1=x1-p(2);
x2=x2-p(3);

% rotate by specified angle
s=sin(p(6)); c=cos(p(6));
x1r=+x1*c+x2*s;
x2r=-x1*s+x2*c;

y=p(1)*exp(-(x1r.^2/(2*p(4)^2)+x2r.^2/(2*p(5)^2)))+p(7);

return