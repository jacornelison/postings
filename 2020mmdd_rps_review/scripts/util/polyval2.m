function y=polyval2(p,x)
% y=polyval2(p,x)
%
% I had poly function to return polynomial values
% Matlab uses poly for something else leading to conflict
% Matlab has polyval function which is same as mine except
% order of p vector is reversed - highest to lowest order.
% To avoid too much fussing add this for now...

y=polyval(flipud(p(:)),x);

return
