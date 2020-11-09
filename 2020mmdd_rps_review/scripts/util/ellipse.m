function [x,y]=ellipse(x0,y0,a,b,phi,npts,uflag)
% ellipse(x0,y0,a,b,phi,npts,uflag)
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

if(all(size(a)==1))
  a=a*ones(size(x0));
end
if(all(size(b)==1))
  b=b*ones(size(x0));
end
if(all(size(phi)==1))
  phi=phi*ones(size(x0));
end
  
% Make sure x0,y0,a,b,phi are row vectors
x0=x0(:)'; y0=y0(:)'; a=a(:)'; b=b(:)'; phi=phi(:)'; 

% Make column vector of angles
s=2*pi/npts; t=[0:s:2*pi-s]';

% Make arrays of both
x0=repmat(x0,size(t,1),1);
y0=repmat(y0,size(t,1),1);
a =repmat(a, size(t,1),1);
b =repmat(b, size(t,1),1);
phi=repmat(phi, size(t,1),1);
t =repmat(t, 1,size(x0,2));

r=sqrt(a.^2.*b.^2./(b.^2.*cos(t-phi).^2+a.^2.*sin(t-phi).^2));
[x,y]=pol2cart(t,r);

x=x+x0;
y=y+y0;

% Close the circles if wanted
if(uflag==1)
  x=[x;x(1,:)];
  y=[y;y(1,:)];
end
