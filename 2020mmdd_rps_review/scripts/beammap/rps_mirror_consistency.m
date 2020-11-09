function rps_moon_consistency()

dirname = 'rps_data/2019_mirror/';
load([dirname 'b3rpsfiles'])
load([dirname 'fitdata'])

ch = 697;

% bs.r = p.r(ch);
% bs.theta = p.theta(ch);
% bs.chi = p.chi(ch);
% bs.chi_thetaref = p.chi_thetaref(ch);
% bs.drumangle = 0;


bs.r = 0;
bs.theta = 0;
bs.chi = 0;
bs.chi_thetaref = 0;
bs.drumangle = 0;

% load tods
%tods = rps_load_tods(sch,[1 4 7 12],[118 119 121 122],p,rpsopt,dirname,697, true);
ch=1991
tods = rps_load_tods(sch,[3 5 7 10],((([7 4 3 8]-1)*13)+1),p,rpsopt,dirname,1991, true);


% OG params, as fit by moon
%   rpsopt.mirror.height = 0.9540
%   rpsopt.mirror.tilt = 44.6870
%   rpsopt.mirror.roll = 0.1480

rpsopt.mirror.tilt = 44.6870
rpsopt.mirror.roll = 0.1480
%tilts = 44:0.1:45;
%rolls = 0:0.05:0.2;


hold off
for k = 1%:length(rolls)
%rpsopt.mirror.tilt = tilts(k);
%rpsopt.mirror.roll = rolls(k);

bparam = NaN(4,2);
for i = 1:4
   
[az, el, pa] = ...
    keck_beam_map_pointing(tods{i}.az,tods{i}.el,tods{i}.dk, rpsopt.mount, rpsopt.mirror, rpsopt.source, bs,'NoSource');

dx = 0.1;
dy = dx;
xbound = (max(az)-min(az))/2*0+15;
xbin = [-xbound:dx:xbound]+median(az);
ybound = (max(el)-min(el))/2*0+15;
ybin = [-ybound:dx:ybound]+median(el);

% Plot Beam Map
chind = tods{i}.ch==ch;

beam_map = grid_map(az,el, tods{i}.todcos(:,chind),xbin,ybin);

fig = figure(i)
h = imagesc(xbin, ybin, beam_map);

[m, ind] = max(tods{i}.todcos(:,chind));
bparam(i,1) = az(ind);
bparam(i,2) = el(ind);

end

%plot(bparam(:,1),bparam(:,2),'k')
%hold on
end


for i = 1:4
[r, theta, psi] = ...
    keck_beam_map_pointing(tods{i}.az,tods{i}.el,tods{i}.dk, rpsopt.mount, rpsopt.mirror, rpsopt.source, bs);

x = 2 * sind(r / 2) .* cosd(theta) * 180 / pi;
y = 2 * sind(r / 2) .* sind(theta) * 180 / pi;

dx = 0.1;
dy = dx;
xbound = 15;
xbin = [-xbound:dx:xbound];
ybin = xbin;
xbin = xbin+mean(x);
ybin = ybin+mean(y);

% Plot Beam Map
chind = tods{i}.ch==ch;
fig = figure(i)

beam_map = grid_map(x,y, tods{i}.todcos(:,chind),xbin,ybin);

h = imagesc(xbin, ybin, beam_map);


x0 = 2 * sind(bs.r / 2) .* cosd(bs.theta) * 180 / pi;
y0 = 2 * sind(bs.r / 2) .* sind(bs.theta) * 180 / pi;
[m, ind] = max(tods{i}.todcos(:,chind));
x(ind)
y(ind)


end


