function y = gauss(p,x)
% y = gauss(p,x)
%
% Compute values at x on a Gaussian shaped curve
% with parameters:
%
% p(1) = peak height
% p(2) = mean
% p(3) = sigma
%
% p(4) = optional zero offset 
%
% eg: x=[-3:0.1:3]; plot(x,gauss([10,0.5,1.3],x));
%
% NB - the parameter order is the reverse of
% that expected by Matlab funcs like quad - but swapping it doesn't
% help as quad refuses to pass forward anything but
% scalar parameters.

if(length(p)<4)
  p(4)=0;
end
  
y=p(1)*exp(-(x-p(2)).^2/(2*p(3)^2))+p(4);
