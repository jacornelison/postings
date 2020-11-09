function [az_ap el_ap]=sidelobe_parallax(az,el,dk,az_com,el_com,dk_com,p,ind,bs)
% bs is optional argument to center the maps on the boresight rather than
% the pixel center
if ~exist('bs','var')
  bs = [];
end
if isempty(bs)
  bs = 0;
end

% az and el are the boresight az and el as derived from the full pointing
% model

% az_com and el_com are the command coordinates

ntot=length(p.r);
nch=length(p.r);
nsamp=length(az);

% az and el are the topocentric horizontal coordinates of the boresight

% our origin is the intersection of the elevation and azimuthal axes
% in the following I assume the dk axis also intersects the origin.

% find the physical location of the pixel w.r.t. the origin when the 
% telescope is vertical and at dk=0:
plate_scale = 100; % degrees per m
[p0(:,1) p0(:,2) p0(:,3)]=pol2cart((90-p.theta)*pi/180,p.r/plate_scale,0);

% rotate that pixel's physical location to our current azimuth and elevation before reflection:
% Here we use the command coordinates:
alpha=-az_com;
beta=270+el_com;
gamma=dk_com;

% Construct the rotation matrix for the focal plane
R=make_rotations(alpha,beta,gamma);

% Calculate p0, the vector from the origin to the pixel centroid at this
% az,el,dk:
p0=repmat(p0,[1,1,nsamp]);
for ii=1:nch
  p0(ii,:,:)=multiprod(p0(ii,:,:),R);
end

% Calculate k hat, the unit vector of the beam centroid as it exits the
% telescope.  Here we use the full pointing model az and el:

if bs
  % center on boresight:
  p.r(:) = 0;
  p.theta(:) = 0;
end

k = calc_k_hat(az,el,p,dk);

% scale k by the source distance:
k = k * 9.1; % source distance in m 

k = p0 + k;

% Calculate the apparent source azimuth and elevation:

clear p0
[theta phi]=cart2sph(squeeze(k(:,1,:)),squeeze(k(:,2,:)),squeeze(k(:,3,:)));
% az_ap=zeros(nsamp,ntot); el_ap=zeros(nsamp,ntot);
az_ap=-theta'*180/pi;
el_ap=phi'*180/pi;

% prevent wrapping:
az_ap = mod(az_ap,360);
 
function k_hat = calc_k_hat(az,el,p,dk)
nch=length(p.r);
ntot=length(p.r);
nsamp=length(az);
% Calculate k hat, the unit vector of the beam centroid as it exits the
% telescope:
lat=zeros(nch,nsamp);
lon=zeros(nch,nsamp);

for ii=1:nch
  [lat(ii,:) lon(ii,:)]=reckon(el,az,p.r(ii),90-p.theta(ii)+dk);
end
k_hat=zeros([nch,3,nsamp]);
% phi_px=-lon;
% theta_px=90-lat;

k_hat(:,1,:)=cosd(-lon).*sind(90-lat);
k_hat(:,2,:)=sind(-lon).*sind(90-lat);
k_hat(:,3,:)=cosd(90-lat);

