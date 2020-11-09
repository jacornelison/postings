function d=invpointing_model(d,model,mirror,fs)
% d=invpointing_model(d,model)
%
% This is the algorithm in Cynthia's inverse_pointing.c
% It is an exact inverse of the original Clem pointing_model.m which
% John Kovac used to fit star data to derive model parameters for BICEP1.
% (The last statement is checked true for az/el)
% Note that it only handles az/el tilt and zero points.
%
% This algorithm checks out the same as Kiwon's

% get the raw values
if ~exist('fs','var')
  enc_mul= double(d.antenna0.tracker.encoder_mul(1,:))./360;
else
  enc_mul= double(d.antenna0.tracker.encoder_mul(fs.s(1),:))./360; 
end

if nargin<3
  mirror=0;
end

if mirror==1 
  if strcmp(get_experiment_name(), 'bicep3')
    if ~exist('fs','var')
      enc_sign= double(d.antenna0.tracker.encoder_sign(1,:));
    else
      enc_sign= double(d.antenna0.tracker.encoder_sign(fs.s(1),:));
    end    
    enc_mul(2)=enc_mul(2)*enc_sign(2);
  end
  enc_mul(3)=-enc_mul(3);
end

az=double(d.antenna0.pmac.fast_az_pos)/enc_mul(1);
el=double(d.antenna0.pmac.fast_el_pos)/enc_mul(2);
dk=double(d.antenna0.pmac.fast_dk_pos)/enc_mul(3);

% For az/el Cynthia subtracts the online offsets and then adds them
% again, but for dk it remains subtracted.
% We may wish to revisit this choice we made initially for BICEP1, and instead 
% allow an offline parameter for dk_zero in the same
% way AZ and EL zeros are included in the pointing model - JMK

samprate=length(d.antenna0.pmac.fast_az_pos)/size(d.antenna0.tracker.encoder_off,1);
dk_off=cvec(repmat((double(d.antenna0.tracker.encoder_off(:,3))/3.6e6)',samprate,1));
dk=dk-dk_off;

% (1) Encoder zeros
az=az-model.az_zero;
el=el-model.el_zero;

% convert to rad
az=az*pi/180;
el=el*pi/180;
dk=dk*pi/180;

el_tilt=model.el_tilt*pi/180;
az_tilt_lat=model.az_tilt_lat*pi/180;
az_tilt_ha=model.az_tilt_ha*pi/180;

% (2) El tilt: correct el first, then az
elraw=el;
el=sin(el).*cos(el_tilt);
el(el>+1) = +1;
el(el<-1) = -1;
el=asin(el);

a=tan(el).*tan(el_tilt);
a(a>+1)=+1;
a(a<-1)=-1;
az=az+asin(a);

dk_offset=sin(el_tilt)./cos(el);
dk_offset(dk_offset>+1)=+1;
dk_offset(dk_offset<-1)=-1;
dk=dk+asin(dk_offset);

% (3) Az tilt: correct dk first, then az, el
x = az_tilt_lat;
y = az_tilt_ha;
top = sin(y.*cos(x));
bot = cos(y.*cos(x)) .* sin(x);
tilt_az = 0 - pi - atan2(top,bot);
tilt_mag = acos(cos(x) .* cos(y*cos(x)));

azpt=az;
w=tilt_az - pi/2 - az;

top = sin(tilt_mag).*cos(w);
bot = cos(tilt_mag).*cos(el) - sin(el).*sin(tilt_mag).*sin(w);
dk_offset = atan2(top, bot);
dk = dk + dk_offset;
	  
top = cos(el).*cos(w);
bot = sin(el).*sin(tilt_mag) - cos(el).*cos(tilt_mag).*sin(w);
a = atan2(top, bot);
az = tilt_az - a;

% +/- 2pi until we're close to the original az value
margin=10*pi/180;
ind=az<azpt-margin; az(ind)=az(ind) + 2*pi;
ind=az>azpt+margin; az(ind)=az(ind) - 2*pi;
ind=az<azpt-margin; az(ind)=az(ind) + 2*pi;
ind=az>azpt+margin; az(ind)=az(ind) - 2*pi;

el = sin(el).*cos(tilt_mag) + cos(el).*sin(tilt_mag).*sin(w);
el(el>+1) = +1;
el(el<-1) = -1;
el=asin(el);
el(elraw>pi/2)=pi-el(elraw>pi/2);
az=az*180/pi;
el=el*180/pi;
dk=dk*180/pi;

d.pointing.hor.az=az;
d.pointing.hor.el=el;
d.pointing.hor.dk=dk;

return
