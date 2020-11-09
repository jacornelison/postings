function [ra,dec]=coord2point(ra0,dec0,r,theta_in)
%[ra,dec]=coord2point(ra0,dec0,r,theta)
%
%opposite of reckon; go from ra0 dec0 of bolometer to the 
%ra dec pointing of the array

% go from SPUD fp coordinates to BICEP theta
theta=theta_in+90;

r=r*pi/180;
theta=theta*pi/180;
decrad=cvec(dec0*pi/180);
rarad=cvec(ra0*pi/180);

ncoord=length(decrad);

maxdec = pi/2 - asin(sin(theta).*sin(r));


% Calculated quantities are a function of declination only so compute
% them for unique values of dec
[decrad,i,j]=unique(decrad);

%This is the basic operation, but diagonalizing a large vector takes 
%too much memory for a large array of coordinates, so...
%arg=diag(1./sin(pi/2-decrad))*repmat(sin(r).*sin(theta),ncoord,1);
aa=1./sin(pi/2-decrad);
bb=sin(r)*sin(theta);
arg=aa*bb;

B=asin(arg);
%maxdec=repmat(maxdec,ncoord,1);
%B(decgrid>maxdec)=NaN;
B(abs(arg)>1)=NaN;

%Napier's analogies for spherical trig
term(:,1)=tan((B-theta)/2);
term(:,2)=sin((r+(pi/2-decrad))/2);
term(:,3)=sin((r-(pi/2-decrad))/2);

%numerical precision stuff
term(abs(term)<1e-10)=0;
arg=term(:,1).*term(:,2)./term(:,3);

phi=pi-2*atan(arg);

%law of halversines
arg = sin((r-(pi/2-decrad))/2).^2 + sin(r).*sin(pi/2-decrad).*...
      sin(phi/2).^2;
a = 2*asin(sqrt(arg));

ra=(rarad-B(j))*180/pi;
dec=(pi/2-a(j))*180/pi;

sz=size(ra0);
ra=reshape(ra,sz);
dec=reshape(dec,sz);

return























