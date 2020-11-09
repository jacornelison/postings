function mva=pointing_model(model,ide)
% Convert ide.az,.el,.dk to mva.az,.el using pointing model
% parameters model
%
% Note that due to offset between dk rot axis and radio
% point axis mva.az,el depends on the deck angle of observation
%
% model parameters:
%  1 flex_cos     antenna0.tracker.flexure(1)
%  2 flex_sin     antenna0.tracker.flexure(2)
%  3 az_tilt_ha   antenna0.tracker.tilts(1)
%  4 az_tilt_lat  antenna0.tracker.tilts(2)
%  5 el_tilt        antenna0.tracker.tilts(3)
%  6 collim_x     antenna0.tracker.fixedCollimation(1)
%  7 collim_y     antenna0.tracker.fixedCollimation(2)
%  8 collim_mag   antenna0.tracker.polarCollimation(1)
%  9 collim_dir   antenna0.tracker.polarCollimation(2)
% 10 az_zero      antenna0.tracker.encoder_off(1)
% 11 el_zero      antenna0.tracker.encoder_off(2)

if(~isfield(ide,'dk'))
  ide.dk=zeros(size(ide.az));
end

% All input params in degrees - convert to rad
d2r=pi/180;
model=model*d2r;
az=ide.az*d2r;
el=ide.el*d2r;
dk=ide.dk*d2r;

% Flexure - note QUaD is cos/sin versus sin/cos for SZA
el=el-model(1)*cos(el);
el=el-model(2)*sin(el);

% Az tilt

aztiltmode=2;

switch aztiltmode
  case 0
  case 1
    % Clem's version

    % astronomy az is clockwise from north
    % but want to work in frame with az anticlock from x as for matlab
    % sph2cart function etc. Also x-east, y-north, z-up seems most
    % natural.
    az=-az+pi/2;
    
    % Get normal vector to tilt plane
    c=cross([1,0,tan(model(3))],[0,1,tan(-model(4))]);
    % Find magnitude and dir of tilt
    phi=atan2(c(2),c(1));
    theta=atan(sqrt(c(1)^2+c(2)^2)/c(3));
    % Apply rotation
    [x,y,z]=sph2cart(az,el,ones(size(az)));
    [x,y,z]=rotaboutz(x,y,z,phi);   % rotate to x along tilt dir
    [x,y,z]=rotabouty(x,y,z,theta); % rotate by tilt angle
    [x,y,z]=rotaboutz(x,y,z,-phi);  % rotate back
    [az,el]=cart2sph(x,y,z);
    
    % Convert back to az clock from north
    az=-az+pi/2;

  case 2
    % Currently in the online code
    
    x = model(4);
    y = model(3);
    % Precompute trig terms.
    cos_x = cos(x);
    sin_x = sin(x);
    ycos_x = y * cos(x);
    cos_ycosx = cos(ycos_x);
    sin_ycosx = sin(ycos_x);
    % Compute the numerator and denominator of the atan2() that is used
    % to compute the azimuth of the tilt.
    top = sin_ycosx;
    bot = cos_ycosx * sin_x;
    % Compute the azimuth of the tilt.
    % Using -halfpi here seems to be arbitrary
    tilt_az = 0 - pi - atan2(top, bot);
    % Compute the magnitude of the tilt.
    tilt_mag = acos(cos_x * cos_ycosx);
    % Compute the direction between the azimuth of the source and the azimuth
    % of the axis around which the tilt is directed.
    w = tilt_az - pi/2 - az;
    % Precompute trig terms.
    sin_w = sin(w);
    cos_w = cos(w);
    sin_mag = sin(tilt_mag);
    cos_mag = cos(tilt_mag);
    % Compute the new target elevation.
    sin_el = sin(el) .* cos_mag - cos(el) .* sin_mag .* sin_w;
    el = asin(sin_el);
    cos_el = cos(el);
    % Compute the new target azimuth.
    top = cos_w .* cos_el;
    bot = -(cos_mag .* sin_w .* cos_el + sin_mag .* sin_el);
    az = tilt_az - atan2(top,bot);
end

% El tilt
% There is no way to do this with vector rotation as axes not perp.
% Below is taken from Tim/Martin
el = asin(sin(el)./cos(model(5)));
el(imag(el)~=0)=NaN;
az=az-asin(tan(model(5)).*sin(el)./cos(el));
az(imag(az)~=0)=NaN;

% Cross-el Collimation
% There is no way to do cross-el collimation with vector rotation
% Below is taken from Tim/Martin for the case of deck=90
az=az-asin(sin(model(6))./cos(el));
az(imag(az)~=0)=NaN;
el=asin(sin(el)./cos(model(6)));
el(imag(el)~=0)=NaN;

% El collimation
% Exactly the same as el zero point shift
el=el+model(7);

% mag/az collimation of third axis
deck=dk+model(9);
az=az-asin(sin(model(8)).*sin(deck)./cos(el));
el=asin(cos(model(8)).*sin(el) - sin(model(8)).*cos(deck).* ...
    sqrt(cos(model(8)).^2 - sin(el).^2 + cos(deck).^2 * sin(model(8)).^2) ./ ...
    (cos(model(8)).^2 + cos(deck).^2.*sin(model(8)).^2));

% Encoder zero points
az=az+model(10);
el=el+model(11);

% put phi into 0 to 2pi range
ind=az<0; az(ind)=az(ind)+2*pi;
ind=az>2*pi; az(ind)=az(ind)-2*pi;

% Convert back to deg
mva.az=az/d2r;
mva.el=el/d2r;

return
