function [az_ap el_ap pa]=projection_3d(az,el,dk,az_com,el_com,dk_com,p,ind,moon,do_pol,mirror)
% projection_3d(az,el,dk,az_com,el_com,dk_com,p,ind,moon,do_pol,experiment)

% projection_3d projects the position of each pixel onto a sphere a finite
% distance away from the origin.  To be used with beammapping data as part
% of a fully general, 3D pointing model.

% ver 1 - RWA

% Arguments:
% az, el, dk: boresight az and el as derived from the full pointing
% model (get_pointing_model) 
% az_com, el_com, dk_com: command coordinates
% p, ind: array info from get_array_info
% moon: set to 1 to indicate moon data
% do_pol: set to 1 to include polarization angle calculation
% experiment: Specify experiment.  Default 'bicep2'

if ~exist('experiment','var')
  experiment='bicep2';
end
if ~exist('do_pol','var')
  do_pol=0;
end


% az_com and el_com are the command coordinates

% Nov 2010 maps:
tilt=41.2;
roll=0.2;
roll=-0.2;

% Nov 2011 maps:
% measured the tilt to be 40.9 on nov 28
% tilt=40.9;

% roll = 0.155: tilt = 41.372: A(2) = .0059, A(3) = -.0341 @ dk = 90

% roll = .165 tilt = 41.372: A(2) = 0.1588, A(3) = -0.0407 dk = 0 
% roll = 0    tilt = 41.372: A(2) = -0.110, A(3) = -0.0403 @ dk = 0 
% roll = .083 tilt = 41.372: A(2) = 0.0251, A(3) = -0.0408 @ dk = 0 
% roll = .067 tilt = 41.372: A(2) = 0.0021, A(3) = -0.0408 @ dk = 0 
% roll = .067 tilt = 41.35: A(2) = -0.0014, A(3) = -0.1030 @ dk = 0 
% roll = .067 tilt = 41.387: A(2) = -0.0662, A(3) = -0.0401 @ dk = 0 px 28
% roll = .104 tilt = 41.402: A(2) = 0.0506, A(3) = .1386 @ dk = 0 px 307
% roll = .06 tilt = 41.316: A(2) = -0.0298, A(3) = -.0491 @ dk = 0 px 307
% roll = .06 tilt = 41.33: best guess at minimization for 307
tilt=41.33;
roll=0.06;
% compare dk 0 and dk 180
% Moon data, nov 28 2011
% roll = .06 tilt = 41.33 delta_x = -1.8 delta_y =  -10.3
% roll = .06 tilt = 41.30 delta_x = -3.4 delta_y =  -17.8
% roll = .06 tilt = 41.36 delta_x = -2.4 delta_y =  -3.6
% roll = .06 tilt = 41.375 delta_x = -2.4 delta_y =  0.31
% roll = .065 tilt = 41.375 delta_x = -2.4 delta_y =  0.31


tilt=41.375;
roll=0.07;

if ~exist('mirror','var')
  mirror=[];
end
if ~isempty(mirror)
  tilt = mirror.tilt;
  roll = mirror.roll;
end

arm=2.98; % distance from origin to where the boresight intersects the mirror, in m

ntot=length(p.r);
nch=length(ind.a);
nsamp=length(az);

% az and el are the topocentric horizontal coordinates of the boresight

% our origin is the intersection of the elevation and azimuthal axes
% in the following I assume the dk axis also intersects the origin.

% find the physical location of the pixel w.r.t. the origin when the 
% telescope is vertical and at dk=0:
plate_scale = 100; % degrees per m
[p0(:,1) p0(:,2) p0(:,3)]=pol2cart((90-p.theta(ind.a))*pi/180,p.r(ind.a)/plate_scale,0);
% rotate that pixel's physical location to our current azimuth and elevation before reflection:
% Here we use the command coordinates:
alpha=-az_com;
beta=270+el_com;
gamma=dk_com;

% Construct the rotation matrix for the focal plane
R=make_rotations(alpha,beta,gamma);
Rp = R;

% Calculate p0, the vector from the origin to the pixel centroid at this
% az,el,dk:
p0=repmat(p0,[1,1,nsamp]);
porg = p0;
for ii=1:nch
  p0(ii,:,:)=multiprod(p0(ii,:,:),R);
end

% Calculate b, the vector connecting the origin to the location at which
% the boresight intersects the mirror. Here we also use the command coordinates:
% found via solid model comparison
theta=-az_com;
phi=el_com;
b=arm*[cosd(theta).*cosd(phi) sind(theta).*cosd(phi) sind(phi)]';

% Calculate n_hat, the unit normal vector to the mirror:
alpha=180-az_com;
beta=180+tilt-el_com;
gamma=roll;
N=[1 0 0];

% Construct the rotation matrix for the mirror
R=make_rotations(alpha,beta,gamma);
% Perform the rotation:
n_hat=squeeze(multiprod(repmat(N,[1,1,length(az)]),R));

% Calculate k hat, the unit vector of the beam centroid as it exits the
% telescope.  Here we use the full pointing model az and el:

[k alpha_px beta_px gamma_px] = calc_k_hat(az,el,p,dk,ind);

% Construct the rotation matrix for the mirror
if do_pol
  % betax:
  betax = atand(k(:,1,:)./k(:,3,:));
  betay = atand(k(:,2,:)./k(:,3,:));

  P=zeros(ntot,3,nsamp);
  kt = zeros(nch,3,nsamp);
  z = zeros(size(gamma_px(ii,:)));
  kt(:,3,:) = 1;
  P(ind.a,1,:)=1;
  P(ind.b,2,:)=1;
  for ii=1:nch
%     Rp = make_rotations(alpha_px(ii),beta_px(ii),gamma_px(ii));
    P(ind.a(ii),:,:)=multiprod(P(ind.a(ii),:,:),Rp);
    P(ind.b(ii),:,:)=multiprod(P(ind.b(ii),:,:),Rp);
    dum = cross(P(ind.a(ii),:,:),P(ind.b(ii),:,:),2);
    cp = cross(dum,k(ii,:,:),2);
    norm = repmat((dot(cp,cp,2)).^(1/2),[1,3,1]);
    cp = cp./norm;
    thet = -acos(dot(dum,k(ii,:,:),2));
    R = arbrot(cp,thet);
    P(ind.a(ii),:,:) = multiprod(P(ind.a(ii),:,:),R);
    P(ind.b(ii),:,:) = multiprod(P(ind.b(ii),:,:),R);
%     Rdk=make_rotations(z,z,gamma_px(ii,:));
%     P(ind.a(ii),:,:)=multiprod(P(ind.a(ii),:,:),Rdk);
%     P(ind.b(ii),:,:)=multiprod(P(ind.b(ii),:,:),Rdk);
%     Rt = make_tilt(betax(ii,:),betay(ii,:));
% %     Rt = make_tilt(0,0);
%     P(ind.a(ii),:,:)=multiprod(P(ind.a(ii),:,:),Rt);
%     P(ind.b(ii),:,:)=multiprod(P(ind.b(ii),:,:),Rt);  
%     kt(ii,:,:)=multiprod(kt(ii,:,:),Rt); 
  end
end

% Calculate k vector, the line that connects the pixel (p0) to the point at
% which it intersects the mirror, along the direction k_hat:
for ii=1:nch
  if length(az)==1
    k(ii,:,:) = (dot((b - p0(ii,:)'),n_hat)/dot(k(ii,:)',n_hat)).*k(ii,:);
  else
    k(ii,:,:) = (dot((b - squeeze(p0(ii,:,:))),n_hat)/dot(squeeze(k(ii,:,:)),n_hat)).*k(ii,:,:);
  end
end

% source distance:
if moon
  range=4E8;
else
  range=200; % in meters
end

km=k;

% Reflect k about the mirror surface defined by n_hat and then scale it
% by the source distance
for ii=1:nch
  if length(az)==1
    ks=nnormn(k(ii,:)');
    km(ii,:)=range*nnormn(ks-2*dot(ks,n_hat',1)'*n_hat',1,2);
  else
    ks=squeeze(k(ii,:,:));
    ks=nnormn(ks,1,2);
    km(ii,:,:)=range*nnormn(ks-repmat(2*dot(ks,n_hat,1),3,1).*n_hat,1,2);
  end
end

% The polarization vectors as they leave the mirror:
if do_pol
  P = calc_pol(P,ind,k,n_hat,km,az);
end

% simple vector addition for full ray tracing:
km=km+k+p0;

if do_pol
  % Calculate the polarization projection:
  % Latitude and longitude of the pixel on the sphere is az_ap and el_ap
  P(ind.a,:,:) = P(ind.a,:,:) + km;
  P(ind.b,:,:) = P(ind.b,:,:) + km;
  [theta phi]=cart2sph(squeeze(P(:,1,:)),squeeze(P(:,2,:)),squeeze(P(:,3,:)));
  az_P=-theta'*180/pi;
  el_P=phi'*180/pi;
end

clear p0 k
[theta phi]=cart2sph(squeeze(km(:,1,:)),squeeze(km(:,2,:)),squeeze(km(:,3,:)));
az_ap=zeros(nsamp,ntot); el_ap=zeros(nsamp,ntot);
az_ap(:,ind.a)=-theta'*180/pi;
az_ap(:,ind.b)=az_ap(:,ind.a);
el_ap(:,ind.a)=phi'*180/pi;
el_ap(:,ind.b)=el_ap(:,ind.a);

% prevent wrapping:
az_ap(az_ap>90)=az_ap(az_ap>90)-360;

if do_pol
  % prevent wrapping:
  az_P(az_P>90)=az_P(az_P>90)-360;
  pa=azimuth(el_ap,az_ap,el_P,az_P);
  pa(pa>305)=pa(pa>305)-360;
  pa(pa>145)=pa(pa>145)-180;
  pa(pa<-90)=pa(pa<-90)+180;
end

function [k_hat alpha beta gamma]= calc_k_hat(az,el,p,dk,ind)
nch=length(ind.a);
ntot=length(p.r);
nsamp=length(az);
% Calculate k hat, the unit vector of the beam centroid as it exits the
% telescope:
lat=zeros(nch,nsamp);
lon=zeros(nch,nsamp);

for ii=1:nch
  [lat(ii,:) lon(ii,:)]=reckon(el,az,p.r(ind.a(ii)),90-p.theta(ind.a(ii))+dk);
end
k_hat=zeros([nch,3,nsamp]);

% longitude and azimuth are the same thing:
alpha=-lon;
beta=270+lat;
gamma=dk;
if length(gamma)==size(alpha,2)
  gamma=repmat(reshape(gamma,1,[]),size(alpha,1),1);
elseif length(gamma)==1
  gamma=repmat(gamma,size(alpha));
else
  error(['Don''t know what to do with oddly sized dk variable.']);
end
% gamma=repmat(dk,size(alpha));

k_hat(:,1,:)=cosd(-lon).*sind(90-lat);
k_hat(:,2,:)=sind(-lon).*sind(90-lat);
k_hat(:,3,:)=cosd(90-lat);

function P = calc_pol(P,ind,k,n_hat,km,az)
% Construct a polarization vector, [Ex Ey Ez];
nch = length(ind.e);
nsamp = length(az);

% kfull=zeros(nch,3,nsamp);
% kfull(ind.a,:,:)=k;
% kfull(ind.b,:,:)=k;
% k=kfull;
% 
% kfull=zeros(length(ind.e),3,nsamp);
% kfull(ind.a,:,:)=km;
% kfull(ind.b,:,:)=km;
% km=kfull;

% Reflect the polarization vector off of the mirror:
% The field perpendicular to the mirror is zero:
for ii=1:nch
  if length(az)==1
    Ps=nnormn(P(ii,:)');
    P(ii,:)=nnormn(Ps-2*dot(Ps,n_hat',1)'*n_hat');
  else
    Ps=squeeze(P(ii,:,:));
    Ps=nnormn(Ps);
    P(ii,:,:)=nnormn(Ps-2*repmat(dot(Ps,n_hat,1),3,1).*n_hat);
  end
end

function R = arbrot(cp,theta)
l = cp(:,1);
m = cp(:,2);
n = cp(:,3);
st = sin(theta);
ct = cos(theta);

% Second euler rotation: rotate normal vector by declination beta
% Rb=[cosd(beta) 0 sind(beta); ...
%     0          1 0 ; ...
%    -sind(beta) 0 cosd(beta)];
R=zeros(3,3,length(theta));
R(1,1,:) = l.*l.*(1-ct) + ct;
R(1,2,:) = m.*l.*(1-ct) - n.*st;
R(1,3,:) = n.*l.*(1-ct) + m.*st;
R(2,1,:) = l.*m.*(1-ct) + n.*st;
R(2,2,:) = m.*m.*(1-ct) + ct;
R(2,3,:) = n.*m.*(1-ct) - l.*st;
R(3,1,:) = l.*n.*(1-ct) - m.*st;
R(3,2,:) = m.*n.*(1-ct) + l.*st;
R(3,3,:) = n.*n.*(1-ct) + ct;

