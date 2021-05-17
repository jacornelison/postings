function sim_scan()

addpath(fullfile('Z:','pipeline','util'))
addpath(fullfile('Z:','dev'))
%% Init parameters
% Boresight
bs.r = 0;
bs.theta = 0;
bs.chi = 0;
bs.chi_thetaref = 0;
bs.drumangle = 0;
bs.expt = 'bicep3';

mount = [];
mount.aperture_offr = 0;
mount.aperture_offz = 1;
mount.dk_offy = 0;
mount.dk_offx = 0;
mount.el_tilt = 0;
mount.el_offx = 0;
mount.el_offz = 0;
mount.az_tilt_ha = 0;
mount.az_tilt_lat = 0;

% Mirror
mirror = [];
mirror.height = 2;
mirror.tilt = 45;
mirror.roll = 0;

% Pointing Model
pm = [];
pm.meas_name = {'170203 02:32:23'};
pm.mjd = 5.7787e+04;
pm.flex_cos = 0;
pm.flex_sin = 0;
pm.az_tilt_ha = 0*0.0046;
pm.az_tilt_lat = 0*0.0044;
pm.el_tilt = 0*-0.0286;
pm.collim_x = 0;
pm.collim_y = 0;
pm.collim_mag = 0;%-0.5802;
pm.collim_dir = 0;%143.9816;
pm.az_zero = 127;
pm.el_zero = -85;
pm.err_flex_cos = 0;
pm.err_flex_sin = 0;
pm.err_az_tilt_ha = 0.0614;
pm.err_az_tilt_lat = 0.0610;
pm.err_el_tilt = 0.2295;
pm.err_collim_x = 0;
pm.err_collim_y = 0;
pm.err_collim_mag = 0.0585;
pm.err_collim_dir = 5.8488;
pm.err_az_zero = 0.4323;
pm.err_el_zero = 0.0604;
pm.rms_res = 0.0050;
pm.numpts = 32;

source = [];
source.distance = 200;
source.azimuth = 180;
source.el = 2.5;
source.height = source.distance*tand(source.el);

% Check our new code against inv_pointing model
% Pointing
colormap('default');

de = 0.1;
se = 1;
scans = 10;
da = 0.1;
a0 = 130;
az0 = (-15:da:15) + a0;

md = {};
md2 = {};
for si = 1:scans
    [AZ, EL] = deal([]);
    for ie = 1:(se/de)
    el0 = ie*de+(si-1);
    [A,E] = meshgrid(az0,el0);
    A = reshape(A,[],1);
    E = reshape(E,[],1);
    
    if mod(ie,2) == 1
        A = flipud(A);        
    end
    AZ = [AZ; A];
    EL = [EL; E];    
    end
    DK = zeros(size(AZ))-1.5;
    
    md2{si}.az = AZ;
    md2{si}.el = EL;
    md2{si}.dk = DK;


T = 1:length(AZ);

d = [];
d.antenna0.tracker.encoder_mul = ones(length(AZ),3)*360;
d.antenna0.pmac.fast_az_pos = AZ;
d.antenna0.pmac.fast_el_pos = EL;
d.antenna0.pmac.fast_dk_pos = DK;
d.antenna0.tracker.encoder_off = repmat([0 0 0]*3.6e6,length(AZ),1);
d=invpointing_model(d,pm);
az_ideal = (d.pointing.hor.az);
el_ideal = d.pointing.hor.el;
dk_ideal = (d.pointing.hor.dk);

md{si}.az = az_ideal;
md{si}.el = el_ideal;
md{si}.dk = dk_ideal;
end


figure(1)
clf;
set(gcf,'Position',[500,50,1046,634])
for schind = 1:length(md)
    %[azapp, elapp, pa] = keck_beam_map_pointing(md{schind}.az,md{schind}.el,md{schind}.dk,rpsopt.mount,mirror,source,bs,'NoSource');
    [r, theta, psi] = keck_beam_map_pointing(md{schind}.az,md{schind}.el,md{schind}.dk,mount,mirror,source,bs);
    x = 2 * sind(r / 2) .* cosd(theta) * 180 / pi;
    y = 2 * sind(r / 2) .* sind(theta) * 180 / pi;
    
    subplot(2,3,1)
    hold on;
    plot(md2{schind}.az,md2{schind}.el,'Color',[schind,0,13]./13)
    grid on
    xlim([-20 20]+a0)
    ylim([-1 11])
    xlabel('az (^o)')
    ylabel('el (^o)')
    title('Raw Mount Coordinates')
    
    subplot(2,3,2)
    hold on;
    plot(md2{schind}.dk,md2{schind}.el,'Color',[schind,0,13]./13)
    grid on
    xlim([-3 3]-1.5)
    ylim([-1 11])
    xlabel('dk (^o)')
    ylabel('el (^o)')
    %title('Example RPS Scan: Raw Mount Coordinates')
    
    subplot(2,3,4)
    hold on;
    plot(md{schind}.az,md{schind}.el,'Color',[schind,0,13]./13)
    grid on
    xlim([-20 20]+a0-pm.az_zero)
    ylim([-1 11]-pm.el_zero)
    xlabel('az (^o)')
    ylabel('el (^o)')
    title('Ideal Horizontal Coordinates')
    
    subplot(2,3,5)
    hold on;
    plot(md{schind}.dk,md{schind}.el,'Color',[schind,0,13]./13)
    grid on
    xlim([-3 3]-3)
    ylim([2 8]+82)
    xlabel('dk (^o)')
    ylabel('el (^o)')
    %title('Example RPS Scan: Ideal Horizontal Coordinates')
    
    subplot(2,3,6)
    hold on;
    plot(x,y,'Color',[schind,0,13]./13)
    grid on
    xlim([-1 1]*20)
    ylim([-10 10])
    xlabel('x (^o)')
    ylabel('y (^o)')
    title('Boresight-Centered Coordinates')
    sgtitle('Example Scan')
end


%%

bs.expt = 'bicep3';
%pm = get_pointing_model(sch{1}.t1);

figure(2)
clf;
set(gcf,'Position',[500,50,1046,634])

for schind = 1:length(md)
    [azapp, elapp, pa] = beam_map_pointing_model(md2{schind}.az,md2{schind}.el,md2{schind}.dk,pm,mount,[],[],bs);
    [azapp, elapp, pa] = rectify_az(md2{schind}.az-a0,azapp,elapp,pa);
    [x, y, phi] = beam_map_pointing_model(md2{schind}.az,md2{schind}.el,md2{schind}.dk,pm,mount,mirror,source,bs);
%     x = 2 * sind(r / 2) .* cosd(theta) * 180 / pi;
%     y = 2 * sind(r / 2) .* sind(theta) * 180 / pi;
    subplot(2,3,1)
    hold on;
    plot(md2{schind}.az,md2{schind}.el,'Color',[schind,0,13]./13)
    grid on
    xlim([-20 20]+a0)
    ylim([-1 11])
    xlabel('az (^o)')
    ylabel('el (^o)')
    title('Raw Mount Coordinates')
    
    subplot(2,3,2)
    hold on;
    plot(md2{schind}.dk,md2{schind}.el,'Color',[schind,0,13]./13)
    grid on
    xlim([-3 3]-1.5)
    ylim([-1 11])
    xlabel('dk (^o)')
    ylabel('el (^o)')
    %title('Example RPS Scan: Raw Mount Coordinates')
    
    subplot(2,3,4)
    hold on;
    plot(azapp,elapp,'Color',[schind,0,13]./13)
    grid on
    xlim([-20 20]+a0-pm.az_zero)
    ylim([-1 11]-pm.el_zero)
    xlabel('az (^o)')
    ylabel('el (^o)')
    title('Ideal Horizontal Coordinates')
    
    subplot(2,3,5)
    hold on;
    plot(-1*(pa+90),elapp,'Color',[schind,0,13]./13)
    grid on
    xlim([-3 3]-3)
    ylim([-1 11]-pm.el_zero)
    xlabel('dk (^o)')
    ylabel('el (^o)')
    %title('Example RPS Scan: Ideal Horizontal Coordinates')
    
    subplot(2,3,6)
    hold on;
    plot(x,y,'Color',[schind,0,13]./13)
    grid on
    xlim([-1 1]*20)
    ylim([-1 1]*10)
    xlabel('x (^o)')
    ylabel('y (^o)')
    title('Boresight-Centered Coordinates')
    sgtitle('Example Scan')
end


%% Look at difference in pointing

load('../data/b3_array_info.mat')

pm.az_tilt_ha = 0*0.0046;
pm.az_tilt_lat = 0*0.0044;
pm.el_tilt = 0*-0.0286;

AZ = 130;
EL = 6;
DK = 0;

[x, y, phi] = beam_map_pointing_model(AZ,EL,DK,pm,mount,mirror,source,p);

pm.az_tilt_ha = 0.0046;
pm.az_tilt_lat = 0.0044;
pm.el_tilt = -0.0286;

[x2, y2, phi2] = beam_map_pointing_model(AZ,EL,DK,pm,mount,mirror,source,p);


figure(3)
clf; hold on;

subplot(1,3,1)
hist(x-x2)
grid on
ylabel('N')
xlabel('\Deltax'' (^o)')

subplot(1,3,2)
hist(y-y2)
grid on
xlabel('\Deltay'' (^o)')

subplot(1,3,3)
hist(phi-phi2)
grid on
xlabel('\Delta\phi'' (^o)')

sgtitle({'BICEP3 CMB-Derived Beam Center Residuals','None vs. Applied (New)'})



function [az, el, pa] = rectify_az(azpt,az,el,pa)

nsamp = length(az);
pntxy = [cosd(az) -sind(az) zeros(nsamp,1)];
pntyx = cross(pntxy,repmat([0,0,1],nsamp,1));
pnt = pntxy;
ort = pntyx;
for ind = 1:nsamp
    pnt(ind,:) = (rodmatrix(el(ind),pntyx(ind,:))*pnt(ind,:)')';
    ort(ind,:) = (rodmatrix(pa(ind)-90,-1*pnt(ind,:))*ort(ind,:)')';
end

az = wrapTo360(az);
margin=90;
for i = 1:5
    ind = (azpt-az)> margin; az(ind) = az(ind) + 180;
    ind = (azpt-az)<-margin; az(ind) = az(ind) - 180;
end

pntxy = [cosd(az) -sind(az) zeros(nsamp,1)];
pntyx = cross(pntxy,repmat([0,0,1],nsamp,1));
e3 = cross(pnt,pntyx);
for ind = 1:nsamp
    el(ind) = acosd(dot(pnt(ind,:),pntxy(ind,:)));
    pa(ind) = atan2(dot(ort(ind,:),pntyx(ind,:)),dot(ort(ind,:),e3(ind,:)))*180/pi;
end

