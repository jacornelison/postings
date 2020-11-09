function [x_tic,y_tic,n]=hfill2(x,y,nx,lx,hx,ny,ly,hy,w,opt)
% [x_tic,y_tic,n]=hfill2(x,y,nx,lx,hx,ny,ly,hy,w,opt)
%
% Histogram 2 dimensional data into bins specified by
% number and range.
%
% x,y are coords of events to be binned
% nx and ny are number of x and y bins (default to 10)
% lx,hx are lower and upper edges of x range (default to min/max of x)
% ly,hy are lower and upper edges of y range (default to min/max of t)
% w are event weights (default to 1)
%
% x_tic,y_tic are the bin centers
% n is the bin contents
%
% optional argument opt specifies to find max or min of w values in
% each cell rather than their sum
%
% If no output specified plots the resulting histogram.
%
% eg: hfill2(rand(1,1000),rand(1,1000),10,0,1,10,0,1);

if(nargin<2)
  error('Must provide x,y data to histogram');
end

if(~exist('nx','var'))
  nx=[];
end
if(~exist('lx','var'))
  lx=[];
end
if(~exist('hx','var'))
  hx=[];
end

if(~exist('ny','var'))
  ny=[];
end
if(~exist('ly','var'))
  ly=[];
end
if(~exist('hy','var'))
  hy=[];
end

if(~exist('w','var'))
  w=[];
end

% assume that NaN weights are flagged points
ind=isnan(w);
if(any(ind))
  x=x(~ind);
  y=y(~ind);
  w=w(~ind);
end

% min and max require vector data
x=x(:)';
y=y(:)';

if(isempty(nx))
  nx=10;
end
if(isempty(ny))
  ny=10;
end

if(isempty(lx))
  lx=min(x);
end
if(isempty(hx))
  hx=max(x);
end

if(isempty(ly))
  ly=min(y);
end
if(isempty(hy))
  hy=max(y);
end

if(isempty(w))
  w=ones(size(x));
end

if(isscalar(w))
  w=w*ones(size(x));
end

if(length(x)~=length(y) | length(x)~=length(w))
  error('number of elements in x, y (and w) must be equal');
end

if(exist('opt','var'))
  switch opt
    case 'max'
      opt=1;
    case 'min'
      opt=2;
    otherwise
      error('unknown option')
  end
else
  opt=0;
end

[x_tic,y_tic,n]=hfill2c(double(x),double(y),nx,lx,hx,ny,ly,hy,double(w),opt);

if(nargout==0)
  imagesc(x_tic,y_tic,n); axis xy; colorbar;
end

return
