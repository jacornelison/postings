function [x,y]=circle(x0,y0,r,npts,uflag)
% [x,y]=circle(x0,y0,r,npts,flag)
%
% Generate points on circles of radius r
% centered at x0,y0
% npts arg optional - default 50
% if flag==1 the first point of each circle is
% duplicated as the last (useful for plotting).

if(~exist('npts','var'))
  npts=[];
end

if(~exist('uflag','var'))
  uflag=[];
end

if(isempty(npts))
  npts=50;
end

if(isempty(uflag))
  uflag=0;
end

if(all(size(r)==1))
  r=r*ones(size(x0));
end
if(all(size(x0)==1))
  x0=x0*ones(size(r));
end
if(all(size(y0)==1))
  y0=y0*ones(size(r));
end
  

% Make sure x0,y0,r are row vectors
x0=x0(:)'; y0=y0(:)'; r=r(:)';

% Make column vector of angles
s=2*pi/npts; t=[0:s:2*pi-s]';

% Make arrays of both
x0=repmat(x0,size(t,1),1);
y0=repmat(y0,size(t,1),1);
r =repmat(r, size(t,1),1);
t =repmat(t, 1,size(x0,2));

[x,y]=pol2cart(t,r);
x=x+x0;
y=y+y0;

% Close the circles if wanted
if(uflag==1)
  x=[x;x(1,:)];
  y=[y;y(1,:)];
end
