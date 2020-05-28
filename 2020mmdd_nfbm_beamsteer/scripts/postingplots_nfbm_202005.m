function postingplots_nfbm_202005()
% PLots for posting
load('holylfs/bicep3/beammaps/nfbm/nfbm_bicep3_mount_20170116.mat')

p = nfbm.p;
p = rmfield(p,'expt');


%% Beam Steer

%[p, pind] = get_array_info(nfbm.bmopt.datenum);
p = nfbm.p;
p = rmfield(p,'expt');

tiles = [14,9,4,15,5,16,11,6];
rows =  [3,5,4,6,3,3,4,2];
cols =  [4,3,3,2,6,6,6,4];
pol = 'A';

chind = [];
for tind = 1:length(tiles)
   chind(end+1) = find(p.tile==tiles(tind) & ...
       p.det_row==rows(tind) & p.det_col==cols(tind) & strcmp(p.pol,'A'));
%    chind(end+1) = find(p.tile==tiles(tind) & ...
%        p.det_row==rows(tind) & p.det_col==cols(tind) & strcmp(p.pol,'B'));
end

map = nfbm.quad_map(:,:,chind);
p0 = structcut(p,chind);

prx = 2. * sind(p0.r/2).*cosd(p0.theta)*180.0/pi;
pry = 2. * sind(p0.r/2).*sind(p0.theta)*180.0/pi;

figure(1)
imagesc(nfbm.x,nfbm.y,sum(map,3)/length(chind)); colorbar();


% Adjust distance until it looks good by eye.
distrange = 0:0.5:20;
phaserange = 0:10:360;
x = nfbm.x;
y = nfbm.y;

apt_pos = [0,0];
apt_diam = 0.53/0.0254;
circpos = [apt_pos-apt_diam/2, apt_diam, apt_diam];

% Manually setting the x/y moves
xmove = [-4.1,-4.474, -5.594, -2.794, -3.354, 0, -0.5, -1.114];
ymove = [4.1, 2.5, 0.7978, 4.718, 0.7978, 4.158, 1.918, -0.8822];

newmap = zeros(size(map));
%[xmove, ymove] = deal([]);
planedist = distrange(i);
%planedist=15.5;
%phase = phaserange(i);
phase = 0;

[X,Y] = meshgrid(x,y);

for ind = 1:length(chind)
    %xmove(ind) = planedist*tand(p.r(chind(ind)))*cosd(p.theta(chind(ind))+phase);
    %ymove(ind) = planedist*tand(p.r(chind(ind)))*sind(p.theta(chind(ind))+phase);
    Z = map(:,:,ind);
    newmap(:,:,ind) = interp2(X-xmove(ind),Y-ymove(ind),Z,X,Y);    
    %newmap(:,:,ind) = map(:,:,ind);
end
mask = ~isnan(newmap);
newmap(isnan(newmap)) = 0;

figure(2)
clf


for pltind = [1,2,3,4,5,6,7,8]
    if pltind<5
        subplot(3,3,pltind)
    else
        subplot(3,3,pltind+1)
    end
imagesc(x,y,newmap(:,:,pltind),[0,1.6]);% colorbar();
hold on
% circpos = [(apt_pos+[xmove(pltind), ymove(pltind)])-apt_diam/2, apt_diam, apt_diam];
circpos = [apt_pos-apt_diam/2, apt_diam, apt_diam];
rectangle('Position',circpos,'Curvature',[1 1])
title(p.tile(chind(pltind)))
set(gca,'YDir','normal')
grid on
end

% Cap the coadded beam to make a top hat
coaddmap = (sum(newmap,3)./sum(mask,3));
coaddmean = mean(mean(coaddmap,2),1);
coaddind = coaddmap>coaddmean/4;
coaddmap(coaddind) = coaddmean;

subplot(3,3,5)
imagesc(x,y,coaddmap,[0,1.6]);% colorbar();
hold on
circpos = [apt_pos-apt_diam/2, apt_diam, apt_diam];
rectangle('Position',circpos,'Curvature',[1 1])
set(gca,'YDir','normal')
title(sprintf('dist: %2.2f phase: %2.2f',planedist, phase))

pause(0.2)

end


%% Edge Taper
names = {'ideal','sharp','extended'};
tiles = [14,9,1];
rows =  [3,5,1];
cols =  [4,3,2];

chind = [];
for tind = 1:length(tiles)
   chind(end+1) = find(p.tile==tiles(tind) & ...
       p.det_row==rows(tind) & p.det_col==cols(tind) & strcmp(p.pol,'A'));
end

% Make maps
map = nfbm.quad_map(:,:,chind);
p0 = structcut(p,chind);

prx = 2. * sind(p0.r/2).*cosd(p0.theta)*180.0/pi;
pry = 2. * sind(p0.r/2).*sind(p0.theta)*180.0/pi;

% define the aperture and shift parameters
apt_pos = [0,0];
apt_diam = 19;
xmove = [-4.1, -4.474,-4.5];
ymove = [4.1, 2.5, -2.56];


% Shift the maps
x = nfbm.x;
y = nfbm.y;
[X,Y] = meshgrid(x,y);
newmap = zeros(size(map));

for ind = 1:length(chind)
    Z = map(:,:,ind);
    newmap(:,:,ind) = interp2(X-xmove(ind),Y-ymove(ind),Z,X,Y);    
end

% Define a margin for the edge
aptmargin = 0.1;
[X,Y] = meshgrid(x,y);
X = reshape(X,[],1);
Y = reshape(Y,[],1);
[edge_min, edge_max, edge_avg] = deal(zeros(size(chind)));

margin_diam = apt_diam*(1-aptmargin);
margind = inrange(sqrt(X.^2+Y.^2),margin_diam/2,apt_diam/2);

for ind = 1:length(chind)
marginmap = reshape(newmap(:,:,ind),[],1);
marginmap(~margind) = nan;
% Get the edge taper params

edge_min(ind) = nanmin(marginmap);
edge_max(ind) = nanmax(marginmap);
edge_avg(ind) = nanmean(marginmap);
end

log_edge_min = 10*log10(edge_min./nfbm.A(1,chind));
log_edge_max = 10*log10(edge_max./nfbm.A(1,chind));
log_edge_avg = 10*log10(edge_avg./nfbm.A(1,chind));

% Plot maps and margin overlay
figure(3)
clf
for pltind = 1:length(chind)
subplot(1,3,pltind)
imagesc(x,y,newmap(:,:,pltind),[0,1.6]);% colorbar();
hold on
circpos = [apt_pos-apt_diam/2, apt_diam, apt_diam];
rectangle('Position',circpos,'Curvature',[1 1])
circpos = [apt_pos-margin_diam/2, margin_diam, margin_diam];
rectangle('Position',circpos,'Curvature',[1 1])
xlabel('X (inches)')
ylabel('Y (inches)')
title([sprintf('min %2.2f / max %2.2f / avg %2.2f (dB)', ...
    log_edge_min(pltind), log_edge_max(pltind),log_edge_avg(pltind))])
set(gca,'YDir','normal')
grid on
axis image
end

%% Additional Metrics
% Artificilly make a sharp map
newmap(24:25,9:10,2) = 5;
for ind = 1:length(chind)
marginmap = reshape(newmap(:,:,ind),[],1);
marginmap(~margind) = nan;
% Get the edge taper params

edge_min(ind) = nanmin(marginmap);
edge_max(ind) = nanmax(marginmap);
edge_avg(ind) = nanmean(marginmap);
end

% Calculate additional metrics
maxmin_ratio = edge_max./edge_min;
minmax_diff = edge_max-edge_min;
frac_avg = (edge_avg-edge_min)./minmax_diff;

% Plot maps and margin overlay
figure(3)
clf
for pltind = 1:length(chind)
subplot(1,3,pltind)
imagesc(x,y,newmap(:,:,pltind),[0,1.6]);% colorbar();
hold on
circpos = [apt_pos-apt_diam/2, apt_diam, apt_diam];
rectangle('Position',circpos,'Curvature',[1 1])
circpos = [apt_pos-margin_diam/2, margin_diam, margin_diam];
rectangle('Position',circpos,'Curvature',[1 1])
xlabel('X (inches)')
ylabel('Y (inches)')
title([names{pltind} sprintf(': MMR %2.2f / fracavg %2.2f ', maxmin_ratio(pltind), frac_avg(pltind))])
set(gca,'YDir','normal')
grid on
axis image
end

