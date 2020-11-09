function alpha=chi2alpha(ra,dec,r,theta,chi,thetaref)
% alpha=chi2alpha(ra,dec,r,theta,chi,thetaref)
% 
% turn polarization angle chi into an angle w.r.t north (alpha) 
% given boresight pointing (r,theta), bolometer focal plane coordinates
% (r,theta). Since chi is defined w.r.t. a fiducial theta vector which
% can in principle be different than theta, a reference theta (thetaref) 
% must also be supplied.

% Polarization angle on sky for a given deck angle is a function of
% declination only, so only compute for unique declinations

szdec=size(dec);

r=cvec(r);
theta=cvec(theta);
chi=cvec(chi);
thetaref=cvec(thetaref);
dec=cvec(dec);
ra=cvec(ra);

sz=size(dec);

% Expand everything to be the same size
ndec=length(dec);
nra=length(ra);
nr=length(r);
ntheta=length(theta);
nchi=length(chi);
nthetaref=length(thetaref);

nmx=max([ndec,nra,nr,ntheta,nthetaref]);

if(nra==1)
  ra=repmat(ra,nmx,1);
end
if(ndec==1)
  dec=repmat(dec,nmx,1);
end
if(nr==1)
  r=repmat(r,nmx,1);
end
if(ntheta==1)
  theta=repmat(theta,nmx,1);
end
if(nchi==1)
  chi=repmat(chi,nmx,1);
end
if(nthetaref==1)
  thetaref=repmat(thetaref,nmx,1);
end

% Calulate ra/dec position of detector on the sky. Use best measured 
% bolometer coordinate, theta, not thetaref here. 
[dec_bol,ra_bol]=reckon(dec,ra,r,theta-90);

% Calculate the angle w.r.t. north of the theta vector *at the detector*
az=azimuth(dec_bol,ra_bol,dec,ra);
alpha=az-180+chi;

% Deal specifically with center pixel if r=0
ind=find(r==0);
alpha(ind)=-90+theta(ind)+chi(ind);

% Make alpha positive
alpha(alpha<0)=alpha(alpha<0)+360;

% Apply reference theta correction. This is probably a bit of an approximation.
% The reason for not using thetaref directly in the call to reckon above is
% that the angle of the theta vector w.r.t. north is a strong function of theta
% at high declinations. We want to make our best guess about this angle first before
% adding chi. Using thetaref instead of theta could give a very divergent result 
% near the poles even if the difference between it and theta is small. Adding a correction
% after the fact avoids this.
alpha=alpha+thetaref-theta;

return

