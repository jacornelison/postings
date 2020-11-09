function v=calc_tank_vol(h)
% Probe bottom is 39 mm above bottom of the tank
h=h+39;

for ii=1:length(h)
  v(ii)=calc_tank_body(h(ii));
end

function v=calc_tank_body(h)
% Tank Parameters (all in mm)
h1=302.41;
r1=1676; 
h2=172.11;
R2=840;
r2=210;
r3=1050;
h3=587.7;
ht=2*h1+2*h2+h3; % Total tank height


% Bottom sphere:
if h>h1
  v1=pi*(r1.*h1.^2-h1.^3/3);
else
  v1=pi*(r1.*h.^2-h.^3/3);
  v=v1*1E-6;
  return
end

% Partial toroid:
if h>(h1+h2)
  lt=(r2.^2-h2.^2).^(1/2);
  v2=pi*(-h2.^3/3+h2.*(r2.^2+R2.*(lt+R2))+r2.^2.*R2.*atan(h2/lt));
else
  h2=h1+h2-h;
  lt=(r2.^2-h2.^2).^(1/2);
  v2=566472974-pi*(-h2.^3/3+h2.*(r2.^2+R2.*(lt+R2))+r2.^2.*R2.*atan(h2/lt));
  v=(v1+v2)*1E-6;
  return
end

% Central Cylinder
if h>(h1+h2+h3)
  v3=pi.*r3.^2.*h3;
else
  h3=h-(h2+h1);
  v3=pi.*r3.^2.*h3;
  v=(v1+v2+v3)*1E-6;
  return
end

% Partial toroid:
if h>(h1+h2+h3+h2)
  lt=(r2.^2-h2.^2).^(1/2);
  v4=pi*(-h2.^3/3+h2.*(r2.^2+R2.*(lt+R2))+r2.^2.*R2.*atan(h2/lt));
else
  h4=h-(h3+h2+h1);
  lt=(r2.^2-h4.^2).^(1/2);
  v4=pi*(-h4.^3/3+h4.*(r2.^2+R2.*(lt+R2))+r2.^2.*R2.*atan(h4/lt));
  v=(v1+v2+v3+v4)*1E-6;
  return
end

% Top sphere:
h5=ht-h;
v5=452560821-pi*(r1.*h5.^2-h5.^3/3);
v=(v1+v2+v3+v4+v5)*1E-6;
