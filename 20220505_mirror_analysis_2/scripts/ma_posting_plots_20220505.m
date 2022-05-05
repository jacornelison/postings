function ma_posting_plots_20220505()

%% Compare GCP to Horizons data

% JPL Horizons data
% Either schedule works just fine.
if 1
    load('z:/dev/moon_analysis/timestream_moonsch_4_scan_1.mat')
    fname = 'z:/dev/moon_analysis/horizons_moonsch_4_scan_1_geo.txt';
    f = fopen(fname);
    horz_data_sch = textscan(f,'%f%s%s%f%f%f%f%s%s','delimiter',',','HeaderLines',60);
    fclose(f);
    index = 6;
elseif 1
    load('z:/dev/moon_analysis/timestream_moonsch_54_scan_10.mat')
    fname = 'z:/dev/moon_analysis/horizons_moonsch_54_scan_10.txt';
    f = fopen(fname);
    horz_data_sch = textscan(f,'%f%s%s%f%f','delimiter',',','HeaderLines',60);
    fclose(f);
    index = 4;
end

figdir = 'c:/Users/James/Documents/GitHub/postings/20220505_mirror_analysis_2/figs/';
% Define relevant variables
lat = double(d.antenna0.tracker.siteActual(1,2))/3.6e6; %mas to deg
alt = double(d.antenna0.tracker.siteActual(1,3))/1000; % mm to m
ttrack = double(d.antenna0.tracker.utc); % MJD
%ttrack = double(d.antenna0.time.utcslow);
%ttrack = double(d.array.frame.utc);
ddelt = double(d.antenna0.tracker.refraction(:,3))./3.6e6;
lst = double(d.antenna0.tracker.lst); % in ms
ddist = double(d.antenna0.tracker.equat_geoc(1,3))./1e6; %uAU to AU
ref = double(d.antenna0.tracker.refraction(:,1:2))./3.6e9; % uas to deg
[raxis, z, rcen, sitevel] = get_rcen(lat,alt);

hmjd = horz_data_sch{1}-2400000.5;
ra = horz_data_sch{index};
dec = horz_data_sch{index+1};

dra = double(d.antenna0.tracker.equat_geoc(:,1))./3.6e6;
ddec = double(d.antenna0.tracker.equat_geoc(:,2))./3.6e6;
dazg = wrapTo180(double(d.antenna0.tracker.horiz_geoc(:,1))./3.6e6);
delg = double(d.antenna0.tracker.horiz_geoc(:,2))./3.6e6;
dpag = double(d.antenna0.tracker.horiz_geoc(:,3))./3.6e6;
dazt = wrapTo180(double(d.antenna0.tracker.horiz_topo(:,1))./3.6e6);
delt = double(d.antenna0.tracker.horiz_topo(:,2))./3.6e6;
dpat = double(d.antenna0.tracker.horiz_topo(:,3))./3.6e6;
hra = interp1(hmjd,ra,ttrack);
hdec = interp1(hmjd,dec,ttrack);
t0 = (ttrack-ttrack(1))*84600;

% Compare equat_geoc to interpolated Horizons data

winscale = 1;
fig = figure(1);
fig.Position = [173 158 1175 450*winscale];
clf; hold on;
subplot(3,2,[1,3])
plot(t0,dra,t0,hra)
legend({'GCP-calculated','Matlab-Calc''d'},'Location','northwest')
grid on
ylabel({'RA','[Deg]'})

subplot(3,2,5)
plot(t0,(dra-hra)*3600)
grid on
xlabel('Time (s)')
ylabel({'\DeltaRA','[Arcseconds]'})


subplot(3,2,[2,4])
plot(t0,ddec,t0,hdec)
grid on
ylabel({'Dec','[Deg]'})
subplot(3,2,6)
plot(t0,(ddec-hdec)*3600)
grid on
xlabel('Time (s)')
ylabel({'\DeltaDec','[Arcseconds]'})

sgtitle('Geocentric Equatorial')

fname = fullfile(figdir,'gcp_vs_mat_equat_geoc.png');
saveas(fig,fname)

% Compare horiz_geoc to horizons
[hazg, helg, hpag] = equat_geoc2horiz_geoc(hra,hdec,lat,lst);

winscale = 1;
fig = figure(2);
fig.Position = [1950 0 1175 450*winscale];
clf; hold on;
subplot(3,3,[1,4])
plot(t0,dazg,t0,hazg)
grid on
title({'Az'})
ylabel('[Degrees]')
legend({'GCP-calculated','Matlab-Calc''d'},'Location','northeast')

subplot(3,3,[2,5])
plot(t0,delg,t0,helg)
grid on
title({'El'})
subplot(3,3,[3,6])
plot(t0,dpag,t0,hpag)
grid on
title({'PA'})



subplot(3,3,[7])
plot(t0,wrapTo180(dazg-hazg)*3600)
grid on
xlabel('Time (s)')
ylabel({'\DeltaAz [Arcsec]'})

subplot(3,3,[8])
plot(t0,(delg-helg)*3600)
grid on
xlabel('Time (s)')
ylabel({'\DeltaEl [Arcsec]'})

subplot(3,3,[9])
plot(t0,(dpag-hpag)*3600)
grid on
xlabel('Time (s)')
ylabel({'\DeltaPA [Arcsec]'})

sgtitle('Geocentric Horizontal')


fname = fullfile(figdir,'gcp_vs_mat_horiz_geoc.png');
saveas(fig,fname)

% Compare horiz_topo to horizons


[helt, dpx] = applyParallax(helg,ddist,rcen);
[helt, dref] = applyRefraction(helt,ref);
[hazt, helt] = applyDA(hazg,helt,sitevel);
hpat = hpag;

winscale = 1;
fig = figure(3);
fig.Position = [2900 0 1175 450*winscale];
clf; hold on;
subplot(3,3,[1,4])
plot(t0,dazt,t0,hazt)
grid on
title({'Az'})
ylabel('[Degrees]')
legend({'GCP-calculated','Matlab-Calc''d'},'Location','northeast')

subplot(3,3,[2,5])
plot(t0,delt,t0,helt)
grid on
title({'El'})
subplot(3,3,[3,6])
plot(t0,dpat,t0,hpat)
grid on
title({'PA'})



subplot(3,3,[7])
plot(t0,wrapTo180(dazt-hazt)*3600)
grid on
xlabel('Time (s)')
ylabel({'\DeltaAz [Arcsec]'})

subplot(3,3,[8])
plot(t0,(delt-helt)*3600)
grid on
xlabel('Time (s)')
ylabel({'\DeltaEl [Arcsec]'})

subplot(3,3,[9])
plot(t0,(dpat-hpat)*3600)
grid on
xlabel('Time (s)')
ylabel({'\DeltaPA [Arcsec]'})

sgtitle('Topocentric Horizontal')


fname = fullfile(figdir,'gcp_vs_mat_horiz_topo.png');
saveas(fig,fname)

%% Moon Temperature Profiles

load('z:/dev/moon_analysis/moon_temp_data.mat');

%% Fig 38 Reconstruction

fig = figure(4);
fig.Position = [1950 0 1296 585];
clf; hold on;
scatter(phase+180,lat,20,BT,'filled')
xlim([70, 360+90])
ylim([0,70])
grid on
colorbar()
xticks(70:20:450)
xlabel('Local Phase [Deg]')
ylabel('Lunar Latitude [Deg]')

fname = fullfile(figdir,'fig_38_recon.png');
saveas(fig,fname)


%% Phase/Temp/Convolution Maps

% MJD, , ,raapp,decapp,illum%,angdia,STO
fname = 'z:/dev/moon_analysis/illum_test_moon.txt';
f = fopen(fname);
hd_moon = textscan(f,'%f%s%s%f%f%f%f%f%f%f','delimiter',',','HeaderLines',60);
fclose(f);

% MJD, , ,raapp,decapp,illum%,angdia,STO
fname = 'z:/dev/moon_analysis/illum_test_sun.txt';
f = fopen(fname);
hd_sun = textscan(f,'%f%s%s%f%f','delimiter',',','HeaderLines',59);
fclose(f);
[stos, sto0s, illums, illum0s] = deal(NaN(size(hd_sun{1})));


lims = {[0,90], [100 300], [0 1]};
clab = {'\phi [Deg]','T [K]',''};
valname  = {'phase','temp','conv'};
valtitle = {'Local Phase','Brightness Temperature','Beam-Convolved Map'};

[xdiff, ydiff] = deal(size(hd_moon{1}));
winscale = 0.3;
fig = figure(1);
clf; hold on;
fig.Position = [2400 0 400 400];
tic
ppd = 150;
%ppd = 50;
for dayind = 1:length(hd_sun{1})
    sto0s(dayind) = hd_moon{10}(dayind);
    illum0s(dayind) = hd_moon{6}(dayind);
    ang_diam = hd_moon{7}(dayind)/3600;
    np_ang = hd_moon{end-2}(dayind);
    mpos = [hd_moon{4}(dayind), hd_moon{5}(dayind)];
    spos = [hd_sun{4}(dayind), hd_sun{5}(dayind)];
    
    utc = datestr(datenum(mjd2datestr(hd_moon{1}(dayind)),'yyyy-mmm-dd:HH:MM:SS'),'yyyy/mm/dd HH:MM:SS');
    [Az, El] = RADec2AzEl_no_c(mpos(1),mpos(2),utc,-89.9911,-44.65,2.8);
    %mpos = [Az El];
    [Az, El] = RADec2AzEl_no_c(spos(1),spos(2),utc,-89.9911,-44.65,2.8);
    %spos = [Az El];
    
    
    [cbin, xbin, ybin, STO, phasebin] = moon_illumination(mpos,spos,ang_diam,np_ang,ppd);
    
    ill = sum(length(find(reshape(cbin,[],1)>0.4)))./sum(length(find(~isnan(reshape(cbin,[],1)))))*100;
    np_v1 = [0;1;0]*abs(hd_moon{end-1}(dayind)/3600);
    np_v2 = [0;-1;0]*ang_diam/2*1.2;
    R = rodmatrix(np_ang,[0;0;-1]);
    np_v1 = R*np_v1;
    np_v2 = R*np_v2;
    
    cbin(isnan(cbin)) = 0;
    [X,Y] = meshgrid(xbin,ybin);
    A = gauss2d([1,0,0,0.16,0.16,0,0],X,Y,'std');
    
    M = conv2(cbin,A,'same');
    Mmax = max(max(M,[],1),[],2);
    ind = M'==Mmax;
    xdiff(dayind) = X(ind);
    ydiff(dayind) = Y(ind);%sqrt(X(ind)^2+Y(ind)^2);
    M = M./Mmax;
    
    if 1
    for valind = 1:3
        vals = {asind(phasebin),cbin,M};
        val = vals{valind};
        
        
        
        %plt = imagesc('XData',xbin+mpos(1),'YData',ybin+mpos(2),'CData',cbin');
        plt = imagesc('XData',xbin,'YData',ybin,'CData',val');
        %set(plt,'AlphaData',~isnan(val'))
        c = colorbar();
        c.Title.String = clab{valind};
        caxis(lims{valind})
        title({valtitle{valind},mjd2datestr(hd_moon{1}(dayind)-2400000.5)})
        plot([np_v1(1) np_v2(1)*0],[np_v1(2) np_v2(2)*0],'k')
        plot([np_v1(2) np_v2(1)*0],[-np_v1(1) np_v2(2)*0],'k')
        plot([np_v1(1)],[np_v1(2)],'kx')
        if valind==3
            plot(X(ind),Y(ind),'kx')
        end
        axis image
        grid on
        set(gca,'xdir','reverse')
        xlabel('Moon-Centered X [Deg]')
        ylabel('Moon-Centered Y [Deg]')
        xlim([-1 1]*winscale)
        ylim([-1 1]*winscale)
        stos(dayind) = STO;
        illums(dayind) = ill;
        
        fname = fullfile(figdir,[sprintf('moonmap_day_%02i_',dayind) valname{valind} '.png']);
        saveas(fig,fname)
    end
    end
end
toc


%% B3 Beams

fig = figure(4);
clf; hold on;
fig.Position = [2400 0 400 400];

[X,Y] = meshgrid(xbin,ybin);
A = gauss2d([1,0,0,0.16,0.16,0,0],X,Y,'std');
plt = imagesc('XData',xbin,'YData',ybin,'CData',A);
%set(plt,'AlphaData',~isnan(val'))
c = colorbar();
c.Title.String = '';
caxis(lims{valind})
title('B3 Beam')
%plot([np_v1(1) np_v2(1)*0],[np_v1(2) np_v2(2)*0],'k')
%plot([np_v1(2) np_v2(1)*0],[-np_v1(1) np_v2(2)*0],'k')
%plot([np_v1(1)],[np_v1(2)],'kx')
if valind==3
    %    plot(X(ind),Y(ind),'kx')
end
axis image
grid on
set(gca,'xdir','reverse')
xlabel('X [Deg]')
ylabel('Y [Deg]')
xlim([-1 1]*winscale)
ylim([-1 1]*winscale)

fname = fullfile(figdir,'b3beam.png');
saveas(fig,fname)

%% Az/El offset vs. day.

t = datenum(mjd2datestr(hd_moon{1}-2400000.5),'yyyy-mmm-dd:HH:MM:SS');
ind = [3,11, 57]; % Days we observed.
fig = figure(4);
clf; hold on;
fig.Position = [2400 0 800 400];
subplot(2,1,1)
hold on;
plot(t,xdiff)
plot(t(ind),xdiff(ind),'kx')
grid on
ylabel('Az offset [Deg]')
a = gca;
datetick('x','dd-mmm')
a.XTickLabel = [];

subplot(2,1,2)
hold on
plot(t,ydiff)
plot(t(ind),ydiff(ind),'kx')
grid on
ylabel('El offset [Deg]')
datetick('x','dd-mmm')
sgtitle('Convolved Beam Offsets')

fname = fullfile(figdir,'azeloffs_vs_time.png');
saveas(fig,fname)


% These functions are more or less copied straight from the GCP code:
function [az_geoc, el_geoc, pa] = equat_geoc2horiz_geoc(ra,dec,lat,lst)
% gcp/antenna/control/bicep/Pointing.cc : L254
%Hour angle

ha = (lst/3.6e6/24*360-ra);

% Geocentric Az
cos_el_sin_az = -cosd(dec).*sind(ha);
cos_el_cos_az = sind(dec).*cosd(lat)-cosd(dec).*sind(lat).*cosd(ha);
% az_geoc = (atan2(-cosd(dec).*sind(ha),...
%     sind(dec).*cosd(lat)-cosd(dec).*sind(lat).*cosd(ha))*180/pi);
az_geoc = atan2(cos_el_sin_az,cos_el_cos_az)*180/pi;

% Geocentric El
el_geoc = asind(sind(dec).*sind(lat)+cosd(dec).*cosd(lat).*cosd(ha));

% Parallactic Angle
pa = atan2(sind(ha)*cosd(lat),sind(lat).*cosd(dec)-sind(dec).*cosd(lat).*cosd(ha))*180/pi;

function [el, delta] = applyRefraction(el,refraction)
% gcp/antenna/control/bicep/Refraction.cc:L41;
% Refraction::applyRefraction
A = deg2rad(refraction(:,1));
B = deg2rad(refraction(:,2));

delta = rad2deg((A.*cosd(el).*sind(el).^3+B.*sind(el).*cosd(el).^3) ./...
    (sind(el).^4 + (A.*sind(el).^2+3.*B.*cosd(el).^2)));

el = el+delta;




function [el, delta] = applyParallax(el,dist,rcen)
% gcp/antenna/control/bicep/Site.cc:L111
% Site::applyParallax
delta = -atan2(cosd(el),dist/rcen-sind(el))*180/pi;
el = el + delta;


function [raxis, z, rcen, site_vel] = get_rcen(lat,alt)
% gcp/control/code/share/slalib/geoc.c
a0 = 6378140.0;
f = 1/298.257;
b = (1-f)^2;
au = 1.4959789e11; % meters per au
cp = cosd(lat);
sp = sind(lat);
c = 1./(sqrt(cp.^2+b.*sp.^2));
s = b.*c;
raxis = (a0.*c+alt).*cp./au;
z = (a0.*s+alt).*sp./au;

% gcp/control/code/unix/libunix_src/common/astrom.c:L520
rcen = sqrt(z.^2+raxis.^2);
site_vel = 2*pi*raxis*au/24/3600;


function [az, el] = applyDA(az,el,v)
% gcp/antenna/control/bicep/Site.cc:L128
%Site::applyDiurnalAberration
ct = cosd(el).*sind(az);
th = acosd(ct);
st = sind(th);

x = -cosd(el).*cosd(az);
psi = zeros(size(x));
ind = x~=0 | sind(el)~=0;
psi(ind) = atan2(x(ind),sind(el(ind)))*180/pi;
cp = cosd(psi);
sp = sind(psi);
th = th-v/3e8.*st;
ct = cosd(th);
st = sind(th);
el = asind(st.*cp);
az = wrapTo180(atan2(ct,-st.*sp)*180/pi);



