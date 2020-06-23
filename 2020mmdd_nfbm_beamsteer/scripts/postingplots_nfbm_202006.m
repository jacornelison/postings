function postingplots_nfbm_202005()
% PLots for posting
%%load('holylfs/bicep3/beammaps/nfbm/nfbm_bicep3_mount_20170116.mat')
close all
global pltcnt
pltcnt=1;
load('Z:\nfbm_bicep3_mount_20170116_bmsteer.mat')
%nfbm.x = nfbm.x*25.4;
%nfbm.y = nfbm.y*25.4;
p = nfbm.p;
p = rmfield(p,'expt');
pind = nfbm.ind;

inrange = @(A,B,C) (A > B) & (A<C);

% As-fit NFBM offset, scaling, and background.
%parms = [0.73341,1.4561,0.35424,16.864,0.3475];
parms = [17.9988   38.6928  535.6901   16.2122    0.5488];
xcen = parms(1);
ycen = parms(2);
planedist = parms(3);
phase = parms(4);
lims = [0,1];
loglims = [-12,0];
%% Coadding maps: Before and after shifting

tiles = [14,9,4,15,5,16,11,6];
rows =  [3,5,4,6,3,3,5,2];
cols =  [4,3,3,2,6,6,5,4];
pol = {'A','B'};

chind = [];
for polind = 1:2
for tind = 1:length(tiles)
   chind(end+1) = find(p.tile==tiles(tind) & ...
       p.det_row==rows(tind) & p.det_col==cols(tind) & strcmp(p.pol,pol{polind}));
end
end

%chind = find(p.tile'==11&ismember(1:2640',pind.rgl100));

map = nfbm.quad_map(:,:,chind);
maxes = nanmax(nanmax(map,[],2),[],1);
mins = nanmin(nanmin(map,[],2),[],1);

map = (map-repmat(mins,size(map(:,:,1))))./repmat(maxes-mins,size(map(:,:,1)));

p0 = structcut(p,chind);

prx = 2. * sind(p0.r/2).*cosd(p0.theta)*180.0/pi;
pry = 2. * sind(p0.r/2).*sind(p0.theta)*180.0/pi;


x = -1*nfbm.x;
y = nfbm.y;

apt_pos = [0,0];
apt_diam = 482.6;


%% Look at the masking model

% Poor fit
parms = [0   0  0   0    1];
[xmove, ymove] = get_xy_offsets(parms,p,chind);
model = reshape(nfbm_get_apt_model(parms,x,y,p,chind),length(x),length(y),[]);

[X,Y] = meshgrid(x,y);

newmap = map.*model;
for ind = 1:length(chind)
    newmap(:,:,ind) = interp2(X-xmove(ind),Y-ymove(ind),newmap(:,:,ind),X,Y);
end

coaddmap = newmap;
coaddmap(isnan(coaddmap)) = 0;
coaddmap = sum(coaddmap,3)./sum(coaddmap~=0,3);

plot_some_maps(map.*model,coaddmap,lims,x,y,tiles,rows,cols,'No Fit')
% Good fit
parms = [17.9988   38.6928  535.6901   16.2122    1];
[xmove, ymove] = get_xy_offsets(parms,p,chind);
model = reshape(nfbm_get_apt_model(parms,x,y,p,chind),length(x),length(y),[]);

[X,Y] = meshgrid(x,y);

newmap = map.*model;
for ind = 1:length(chind)
    newmap(:,:,ind) = interp2(X-xmove(ind),Y-ymove(ind),newmap(:,:,ind),X,Y);
end

coaddmap = newmap;
coaddmap(isnan(coaddmap)) = 0;
coaddmap = sum(coaddmap,3)./sum(coaddmap~=0,3);

plot_some_maps(map.*model,coaddmap,lims,x,y,tiles,rows,cols,'Good Fit')
%plot_some_maps(lognewmap,logcoaddmap,loglims,x,y,tiles,rows,cols,'After Shifting')

%% Shifting effects
% Make map with no offsets
coaddmap = sum(map,3)/length(chind);

plot_some_maps(map,coaddmap,lims,x,y,tiles,rows,cols,'No Shifting')
%plot_some_maps(logmap,logcoaddmap,loglims,x,y,tiles,rows,cols)

% Make map with offsets
[xmove, ymove] = get_xy_offsets(parms,p,chind);

[X,Y] = meshgrid(x,y);

newmap = map;
for ind = 1:length(chind)
    newmap(:,:,ind) = interp2(X-xmove(ind),Y-ymove(ind),newmap(:,:,ind),X,Y);
end

coaddmap = newmap;
coaddmap(isnan(coaddmap)) = 0;
coaddmap = sum(coaddmap,3)./sum(coaddmap~=0,3);

plot_some_maps(newmap,coaddmap,lims,x,y,tiles,rows,cols,'After Shifting')
%plot_some_maps(lognewmap,logcoaddmap,loglims,x,y,tiles,rows,cols,'After Shifting')


%% Figuring the aperture size

apt_diam = 0.52*1000;

% Plot maps and margin overlay
figure(pltcnt); pltcnt = pltcnt+1;
clf
set(gcf,'Position',[600,400,300,300])
subplot(1,1,1)
imagesc(x,y,10*log10(coaddmap),[-20,0]);
hold on
circpos = [apt_pos-apt_diam/2, apt_diam, apt_diam];
rectangle('Position',circpos,'Curvature',[1 1])
xlabel('X (mm)')
ylabel('Y (mm)')
set(gca,'YDir','normal')
grid on
axis image
title(sprintf('Aperture: %2.2fmm',apt_diam))
colormap(jet)
hp4 = get(subplot(1,1,1),'Position');
colorbar('Position', [hp4(1)+hp4(3)+0.01  hp4(2)  0.025  hp4(2)+hp4(3)*0.9])

% 
apt_diam = 19*25.4;

% Plot maps and margin overlay
figure(pltcnt); pltcnt = pltcnt+1;
clf
set(gcf,'Position',[600,400,300,300])
subplot(1,1,1)
%imagesc(x,y,coaddmap,[0,1]);% colorbar();
imagesc(x,y,10*log10(coaddmap),[-20,0]);
hold on
circpos = [apt_pos-apt_diam/2, apt_diam, apt_diam];
rectangle('Position',circpos,'Curvature',[1 1])
xlabel('X (mm)')
ylabel('Y (mm)')
set(gca,'YDir','normal')
grid on
axis image
title(sprintf('Aperture: %2.2fmm',apt_diam))
colormap(jet)
hp4 = get(subplot(1,1,1),'Position');
colorbar('Position', [hp4(1)+hp4(3)+0.01  hp4(2)  0.025  hp4(2)+hp4(3)*0.9])



%% Edge Taper Example
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
apt_diam = 19*25.4;
xmove = [-4.1, -4.474,-4.5]*25.4;
ymove = [4.1, 2.5, -2.56]*25.4;


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
figure(pltcnt); pltcnt = pltcnt+1;
clf
set(gcf,'Position',[50,50,900,300])
for pltind = 1:length(chind)
subplot(1,3,pltind)
imagesc(x,y,newmap(:,:,pltind),[0,1.6]);% colorbar();
hold on
circpos = [apt_pos-apt_diam/2, apt_diam, apt_diam];
rectangle('Position',circpos,'Curvature',[1 1])
circpos = [apt_pos-margin_diam/2, margin_diam, margin_diam];
rectangle('Position',circpos,'Curvature',[1 1])
xlabel('X (mm)')
ylabel('Y (mm)')
title([sprintf('min %2.2f / max %2.2f / avg %2.2f (dB)', ...
    log_edge_min(pltind), log_edge_max(pltind),log_edge_avg(pltind))])
set(gca,'YDir','normal')
grid on
axis image
end
colormap(jet)
hp4 = get(subplot(1,3,3),'Position');
colorbar('Position', [hp4(1)+hp4(3)+0.01  hp4(2)+0.08  0.025  hp4(2)+hp4(3)*2.5])

%% Beam Steer

[xmove, ymove] = get_xy_offsets(parms,p,1:2640);
params = nfbm.A;
%params(2,:) = params(2,:);
%params(3,:) = params(3,:);
bmsteer = [params(2,:)'-xmove, params(3,:)'-ymove];
%bmsteer = nfbm.beam_steer;

% Derotate the beam steer
xmovep = (((xmove-xcen))*cosd(phase)+((ymove-ycen))*sind(phase));
ymovep = (((xmove-xcen))*-sind(phase)+((ymove-ycen))*cosd(phase));
bmsteerp(:,1) = bmsteer(:,1)*cosd(phase)+bmsteer(:,2)*sind(phase);
bmsteerp(:,2) = bmsteer(:,1)*-1*sind(phase)+bmsteer(:,2)*cosd(phase);

figure(pltcnt); pltcnt = pltcnt+1;

scale = 0.35;
clf
set(gcf,'position',[500,0,600,600])
hold on
pol = {'A','B'};
clr = {'b','r'};
for polind = 1:2
    % ind = sqrt(sum(bmsteer.^2,2))<2 & (p.mce~=0) & ismember(1:2640,pind.rgl100)' & strcmp(p.pol,pol{polind});
    ind = sqrt(sum(bmsteer.^2,2))<20 & (p.mce~=0) & ismember(1:2640,pind.rgl100)' & strcmp(p.pol,pol{polind});
    quiver(ymovep(ind),xmovep(ind),bmsteerp(ind,2)*scale,bmsteerp(ind,1)*scale,0,'color',clr{polind});
end

limscale = 120;
xlim([-1,1]*limscale); ylim([-1,1]*limscale);
grid on
set(gca,'XDir','reverse')
ylabel('X (mm)')
xlabel('Y (mm)')
legend('Pol A','Pol B')
title(sprintf('Beam Steer (%2.2fx scaling)',scale))
% Per tile beam steer plots
keyboard()

%% Comparing sims to close detectors
%chind = 

function [xmove, ymove] = get_xy_offsets(parms,p,chind)
xcen = parms(1);
ycen = parms(2);
scaling = parms(3);
planedist = parms(3);
phase = parms(4);

% prx = 2. * sind(p.r(chind)/2).*cosd(p.theta(chind))*180.0/pi;
% pry = 2. * sind(p.r(chind)/2).*sind(p.theta(chind))*180.0/pi;
% 
% xmove = scaling*(prx*cosd(phase)-pry*sind(phase))+xcen;
% ymove = scaling*(prx*sind(phase)+pry*cosd(phase))+ycen;

xmove = planedist*tand(p.r(chind)).*cosd(p.theta(chind)+phase)+xcen;
ymove = planedist*tand(p.r(chind)).*sind(p.theta(chind)+phase)+ycen;


function plot_some_maps(map,coaddmap,lims,x,y,tiles,rows,cols,plttitle)
global pltcnt
figure(pltcnt); pltcnt = pltcnt+1;
set(gcf,'Position',[50,50,600,600])
clf

subplot(3,3,5)
imagesc(x,y,coaddmap,lims);% colorbar();
hold on
set(gca,'YDir','normal')
title('Coadd')
axis image
grid on
colormap(jet)
set(gca,'GridAlpha',0.5)

for pltind = 1:8
    if pltind<5
        cind = pltind;
    else
        cind = pltind+1;
    end
    
subplot(3,3,cind)
imagesc(x,y,map(:,:,pltind),lims);% colorbar();
hold on
title([sprintf('Tile %i Rows %i Col %i ',tiles(pltind),rows(pltind),cols(pltind)) 'Pol A'])
set(gca,'YDir','normal')


    switch cind
        case {1,4}
            ylabel('Y (mm)')
        case {3,6}
            %colorbar('east')
        case 7
            ylabel('Y (mm)')
            xlabel('X (mm)')
        case 8
            xlabel('X (mm)')
        case 9
            xlabel('X (mm)')
            %colorbar()
    end
axis image
grid on
set(gca,'GridAlpha',0.5)
colormap(jet)

hp4 = get(subplot(3,3,9),'Position');
colorbar('Position', [hp4(1)+hp4(3)+0.01  hp4(2)  0.025  hp4(2)+hp4(3)*3.2])
sgtitle(plttitle)
end





