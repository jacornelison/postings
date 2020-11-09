function rps_review_posting_plots()
global p p_ind

load('data/b3rpsfiles.mat')
load('data/fitdata_20201013.mat')
addpath('Z:\\b3reduc\b3_analysis\util\')
addpath('Z:\\b3reduc\b3_analysis\beammap\')
load('data/chflags_Jan2017.mat')

p_ind = pind_chflags;



prx = 2 * sind(p.r / 2) .* cosd(p.theta) * 180 / pi;
pry = 2 * sind(p.r / 2) .* sind(p.theta) * 180 / pi;

figdir = 'figs/';
SVPLT = true; % Save plots?
clr = get(groot,'DefaultAxesColorOrder');


% Cut Params
% These are to show what cuts I'm making to the data.
% these should be very loose cuts

% {Param, lower bound, upper bound, title}
param_array = {...
    {prx(fd.ch)-fd.bparam(:,1),-0.3,0.3,'x_{obs}-x_{fit}','dx'},...
    {pry(fd.ch)-fd.bparam(:,2),-0.3,0.3,'y_{obs}-y_{fit}','dy'},...
    {fd.aparam(:,2),-0.05,0.05,'Xpol leakage','xpol'},...
    {atand(tand(p.chi(fd.ch)+p.chi_thetaref(fd.ch)-fd.aparam(:,1))),-4,5,'\phi_{d,0}-\phi_{d,fit}','phid'},...
    {fd.stat,3,3,'Fit status','fit'},...
    {fd.aparam(:,7)./fd.data_rms,1000,8600,'Sig.-to-Noise','snr'}...
    };

fd_cut = structcut(fd,make_cut_index(fd,param_array));

%% Show cut thresholds
close all

for prmind = 1:length(param_array)
    arr = param_array{prmind};
    param = arr{1};
    lb = arr{2};
    ub = arr{3};
    ptitle = arr{4};
    
    fig = figure();
    set(fig,'Position',[900,300,840,420]);
    subplot(1,2,1)
    make_cut_hist(param,lb,ub,ptitle,true)
    
    subplot(1,2,2)
    if ~isnan(lb) & ~isnan(ub)
        make_cut_hist(param,lb,ub,ptitle,false)
        title('zoomed')
    else
        disp(['Skipping: ' ptitle])
    end
    saveas(fig,[figdir 'cutthresh_' arr{5} '.png'])
end


%% Cut channels
% For the Consistency checks, we don't want to include one of the two
% 46.25 boresight angles. Isolate via schedules.
dk = unique(fd_cut.dk);
spec_ind = true(size(fd_cut.ch,1),length(dk));
spec_ind(:,2) = ismember(fd_cut.sch,[7,8,9]);
fd_dk = make_dk_struct(fd_cut,param_array,unique(fd.dk),spec_ind);



%% Consistency check for stat uncert
% Uncomment to split the DK 46 data into their separate obs.
% This give you consistency without the mirror involved.
%
dk = unique(fd_cut.dk);
dk(end+1) = 46.25;
spec_ind = true(size(fd_cut.ch,1),length(dk));
spec_ind(:,2) = ismember(fd_cut.sch,[7,8,9]);
spec_ind(:,end) = ismember(fd_cut.sch,[13,14]);

% dk = [];
% dk(end+1) = 46.25;
% dk(end+1) = 46.25;
% spec_ind = true(size(fd_cut.ch,1),length(dk));
% spec_ind(:,1) = ismember(fd_cut.sch,[7,8,9]);
% spec_ind(:,2) = ismember(fd_cut.sch,[13,14]);


fd_dk = make_dk_struct(fd_cut,param_array,dk,spec_ind);

%% Get rid of weird tiles

cut_tiles = [1,4,8,20];
dk = unique(fd_cut.dk);
ind = ~ismember(p.tile(fd_cut.ch),cut_tiles);
spec_ind = repmat(ind,1,length(dk));
fd_dk = make_dk_struct(fd_cut,param_array,dk,spec_ind);


%% See how rms changes over averaging

parm = 'phip_corr';
parm = 'xpol';
fprintf(['\n\n\nParm: ' parm '\n'])
medsubstd = [];
medsubmed = [];
for dkind = 1:length(fd_dk.dk)
    tilemed = zeros(size(p.gcp));
    for tileind = 1:20
        i = intersect(find(p.tile==tileind),fd_dk.inda);
        tilemed(i) = nanmedian(fd_dk.(parm)(i,dkind));
    end
        
    medsubstd = nanstd(fd_dk.(parm)(fd_dk.inda,dkind)-tilemed(fd_dk.inda));
    medsubmed = nanmedian(fd_dk.(parm)(fd_dk.inda,dkind));%-tilemed(fd_dk.inda));
    disp(['DK ' num2str(fd_dk.dk(dkind)) ': ' num2str(medsubmed) ' -- ' num2str(medsubstd)])
end

phip_mean = nanmean(fd_dk.(parm),2);
    tilemed = zeros(size(p.gcp));
    for tileind = 1:20
        i = intersect(find(p.tile==tileind),fd_dk.inda);
        tilemed(i) = nanmedian(phip_mean(i));
    end

medsubstd = nanstd(phip_mean(fd_dk.inda)-tilemed(fd_dk.inda));
medsubmed = nanmedian(phip_mean(fd_dk.inda));%-tilemed(fd_dk.inda));
disp(['DK All: ' num2str(medsubmed) ' -- ' num2str(medsubstd)])

%% Consistency Check for sys uncert


plt_array = {...
    {'phip_corr',[-0.5 0.5],[0,160],'\phi_{pix}'},...
    {'xpol',[-0.02 0.02],[0,250],'\epsilon_{pix}'}
    };
close all
for pltind = 1:length(plt_array)
    
    make_diff_hist(fd_dk,plt_array{pltind},figdir)
    %make_scatter(fd_dk,plt_array{pltind},figdir)
end

%% Per-tile hists

%close all
xpall = nanmean(fd_dk.xpol,2);
figure(5)
clf
for i = 1:5
    for j = 1:4
        rcind = rowcol2ind(i,j,5,4);
        
        ind = ismember(1:2640,fd_dk.inda)&(p.tile==rcind)';
        
        subplot(5,4,rcind)
        V = xpall(ind);
        edges = -0.025:0.05/12:0.025;
        N = histc(V,edges);
        bar(edges,N,'histc')
        title(['Tile: ' num2str(rcind)])
        set(gca,'XTick',[-0.025:0.05/4:0.025])
        
    end
end



%% Plot averaged pixel estimates

p.expt='bicep3';
plt_array = {...
    {'phip_corr',[-2.0,0.5],[-0.25 0.25],'\phi_{pix}'},...
    {'xpol',[-0.025 0.025],[-0.01,0.01],'\epsilon_{pix}'}
    };
close all
for medind = 1:2
    if medind == 1
        pltsuffix = '';
        plttitle = ' Tile Plot';
        medfact = 0;
        
    else
        pltsuffix = '_medsub';
        plttitle = ' Tile Plot -- Median Subtracted';
        medfact = 1;
        
    end
    for pltind = 1:length(plt_array)
        arr = plt_array{pltind};
        parm = arr{1};
        plimx = arr{1+medind};
        pltlab = arr{4};
        tilemed = zeros(size(p.gcp));
        for tileind = 1:20
            i = intersect(find(p.tile==tileind),fd_dk.inda);
            tilemed(i) = nanmedian(nanmean(fd_dk.(parm)(i),2));
        end
                
        plot_tiles(nanmean(fd_dk.(parm),2)-tilemed*medfact,p,'pair','diff','clab',pltlab,'title',[pltlab plttitle],'clim',plimx)
        colormap jet
        %saveas(gcf,[figdir 'tile_plot_' parm pltsuffix '.png'])
        V = nanmean(fd_dk.(parm),2)-tilemed*medfact;
        disp([parm ', ' pltsuffix ': med ' num2str(nanmedian(V(fd_dk.inda))) ', std ' num2str(nanstd(V(fd_dk.inda)))])
        
        
        
    end
end


%% Load and map tods
load('data\ideal_tod.mat')
addpath('scripts\beammap')
addpath('scripts\util')


%% Show an example of the data
fig = figure(2);
clf
set(fig,'Position',[2200,100,700,500])
mrk = {'o','x'};
ch = [696,697];
pols = ['A','B'];
scaling = 20;
limx = [-195 195];
limy = [-0.1 1.3];
dpix = 0.1;
ybin = -20:dpix:20;
xbin = limx(1):dpix:limx(2);
[X,Y] = meshgrid(xbin,ybin);


for chi = 1:2
subplot(3,1,chi)
%[data,x,y,phi,az,el,dk,quad] = rps_concat_data(tods,ch(chi));
%plot(data/nanmax(data),'.')
az_plot = [];
el_plot = [];
A_plot = [];
for todind = 1:13
    chind = tods{todind}.ch==ch(chi);
    A = tods{todind}.todquad(:,chind);
    az = tods{todind}.az;
    el = tods{todind}.el;
    
    az = (az-nanmedian(az)-0.43)*scaling;
    el = (el-nanmedian(el))*scaling+tods{todind}.rot;
    
    az_plot = [az_plot; az];
    el_plot = [el_plot; el];
    A_plot = [A_plot; A];
end

map = griddata(el_plot,az_plot,A_plot/nanmax(A_plot),X,Y);
imagesc(xbin,ybin,map,[0,1])

%plot(az_plot,A_plot/nanmax(A_plot),'.')
grid on
xlim(limx)
%ylim([-5 5])
a = gca;
a.XTick = [-180:30:180]+15;
a.XTickLabel = [];
a.YTick = [];
a.GridAlpha = 1;
title([pols(chi) ' Detector Rasterset Data'])
end
%ylabel('Source Amplitude (Normalized)','FontSize',12)



subplot(3,1,3)
plot(-100,-100,'Color',clr(1,:));
hold on
plot(-100,-100,'Color',clr(2,:));
for chi = 1:2
chind = find(fd.sch==1 & fd.ch==ch(chi) & fd.row==10);
modcurve = fd.bparam(chind,6:18);
plot(fd.rot(chind,:),modcurve/max(modcurve),['k' mrk{chi}])
grid on
hold on
%plot(fd.rot(chind,:),fd.bparam(chind,6:18))
rot = -180:180;
plot(rot,rps_get_mod_model([fd.phi_d(chind)+90,0,0,0,1/2],rot),'Color',clr(chi,:))

end
xlim(limx)
ylim(limy)
legend('A Detector','B Detector')
xlabel('Source Angle (^\circ)','FontSize',12)
ylabel('Source Amplitude (Normalized)','FontSize',10)
set(gca, 'XTick', [-180:60:180]);
saveas(fig,[figdir 'modcurve.png'])
%%%%%%%%%%%%%%%%%%%
%%% Extra Plots %%%
%%%%%%%%%%%%%%%%%%%
%% Plot tiles of Phi P

close all
% Q / U are just the weighted means of the cos / sin of the angles
% A/B pols are weighted by one.
% Cross pol contribution is weighted by xpol leakage.
for dkind = 1:length(dk)
    
    plot_tiles(fd_dk.phip_corr(:,dkind),p,'pair','diff','clab','\Delta\phi_{pix}','clim',[-2.5 1])
    colormap jet
    
end


%% Plot tiles with median subtraction

close all

for dkind = 1:length(dk)
    tilemed = zeros(size(p.gcp));
    for tileind = 1:20
        i = intersect(find(p.tile==tileind),inda);
        tilemed(i) = nanmedian(fd_dk.phip_corr(i,dkind));
    end
    
    
    plot_tiles(fd_dk.phip_corr(:,dkind)-tilemed,p,'pair','diff','clab','\Delta\phi_{pix}','title','Median Subtracted','clim',[-0.25 0.25])
    colormap jet
    
    
end





%% Per-dl xpol
close all
for dkind = 1:length(fd_dk.dk)
    
    plot_tiles(fd_dk.xpol(:,dkind),p,'pair','diff','clab','\epsilon_P','title','Median Subtracted','clim',[-0.03 0.03])
    colormap jet
end


%% Averaged xpol-leakage

close all
%hist(nanmean(fd_dk.xpol,2),100)
plot_tiles(nanmean(fd_dk.xpol,2),p,'pair','diff','clab','\epsilon_Q','title','Median Subtracted','clim',[-0.03 0.03])
colormap jet



%% Make single raster

close all
fig = figure(1)
clf

todind = 1;
chind = tods{todind}.ch==697;
az = tods{todind}.az;
el = tods{todind}.el;
A = tods{todind}.todquad(:,chind)/nanmax(tods{todind}.todquad(:,chind));

plot3(az,el,A)
xlabel('Azimuth (^\circ)')
ylabel('Elevation (^\circ)')
zlabel('Amplitude')
title('RPS Raster')
grid on
saveas(fig,[figdir 'raster.png'])



%%%%%%%%%%%%%%%%%%%%%%%%
%%% Myriad functions %%%
%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%
function ind = make_cut_index(fd,param_array)
% Makes cut list. If true, channels PASS the cuts.

% Cut bad rasters
% cut_sch = [08 06 07 07 07 07 08 09 09 05 06 06 04 04 05 05 05 09 11];
% cut_row = [11 14 18 02 11 17 17 04 11 11 13 04 17 09 07 10 12 15 02];

%cut_sch = [05 06 06 09];
%cut_row = [18 06 13 11];

cut_sch = [05 06 06 08 09 09];
cut_row = [18 06 13 11 11 15];


cut_raster = true(size(fd.ch));
for i = 1:length(cut_sch)
    cut_raster = cut_raster & ~(fd.sch==cut_sch(i) & fd.row==cut_row(i));
end



arr = param_array{1};
ind = true(size(arr{1}));% & cut_raster;
for prmind = 1:length(param_array)
    arr = param_array{prmind};
    ind = ind & inrange(arr{1},arr{2},arr{3});
    
end

ind = ind&cut_raster;

%%%%%
function make_cut_hist(param,lb,ub,ptitle,pltraw)

if pltraw
    plim = [nanmin(param),nanmax(param)];
    
else
    plim = [lb-abs(lb), ub+ub];
end

edges = plim(1):diff(plim)/100:plim(2);

N = histc(param,edges);
mxn = nanmax(N);


bar(edges,N,'histc');
hold on
plot([1 1]*lb,[0 1]*2*mxn,'r')
plot([1 1]*ub,[0 1]*2*mxn,'r')
grid on
ylim([-0.1, 1.1*mxn])
xlim(plim)
title(ptitle)

%%%%%
function make_diff_hist(fd_dk,arr,figdir)
parm = arr{1};
plimx = arr{2};
plimy = arr{3};
xlab = arr{4};
clr = get(groot,'DefaultAxesColorOrder');
count = 1;
m = [];
s = [];
x = [1 1];
y = [-1 1].*plimy*2;
ndiffs = sum(1:(length(fd_dk.dk)-1));
nrows = ceil(sqrt(ndiffs));
if length(fd_dk.dk)<=2
    ncols = 1;
else
    ncols = floor(sqrt(ndiffs));
end
endrow = get_endrows(ndiffs,nrows,ncols);

dk = unique(fd_dk.dk);
dk_count = ones(size(dk));
dk_names = {};
for dkind = 1:length(fd_dk.dk)
    ind = dk==fd_dk.dk(dkind);
    if dk_count(ind)==1
        dk_names{end+1} = '';
    else
        dk_names{end+1} = ['\_' num2str(dk_count(ind))];
    end
    dk_count = dk_count+double(dk==fd_dk.dk(dkind));
end



fig = figure();
set(fig,'Position',[2200 0 285*ncols 240*nrows]);
for i = 1:length(fd_dk.dk)
for j = (i+1):length(fd_dk.dk)
    if i ~= j
        dkname1 = ['DK ' num2str(floor(fd_dk.dk(i))) dk_names{i}];
        dkname2 = ['DK ' num2str(floor(fd_dk.dk(j))) dk_names{j}];
        subplot(nrows,ncols,count)
        V = fd_dk.(parm)(fd_dk.inda,i)-fd_dk.(parm)(fd_dk.inda,j);
        m(end+1) = nanmedian(V);
        s(end+1) = nanstd(V);
        edges = plimx(1):abs(diff(plimx))/25:plimx(2);
        N = histc(V,edges);
        plt = bar(edges,N,'histc');
        set(plt,'EdgeColor','k','FaceColor',clr(1,:));
        title([dkname1 ' - ' dkname2])
        hold on
        plt = plot(x*m(end),y,'Color',clr(2,:));
        set(plt,'LineWidth',1.5)
        plt = plot(x*m(end)+s(end),y,'Color',clr(5,:));
        set(plt,'LineWidth',1.5)
        plt = plot(x*m(end)-s(end),y,'Color',clr(5,:));
        set(plt,'LineWidth',1.5)
        grid on
        xlim(plimx)
        ylim(plimy)
        set(gca,'XTick',plimx(1):abs(diff(plimx))/4:plimx(2))
        axis square
        
        % Are we in the first column?
        [r,c] = ind2rowcol(count,nrows,ncols);
        if c==1
            ylabel('N','FontSize',14)
        end
        
        % Are we at the bottom most row for that column?
        if r==endrow(c)
            xlabel(['\Delta' xlab],'FontSize',14)
        end
        
        count=count+1;
    end
end
end
disp([parm 'max diff: ' num2str(nanmax(abs(m)))])
disp([parm 'mean std: ' num2str(nanstd(m))])
disp([parm 'min std: ' num2str(nanmin(s))])
saveas(gcf,[figdir 'diffhist_' parm '_sys.png'])

%%%%%
function make_scatter(fd_dk,arr,figdir)
parm = arr{1};
plimx = arr{2};
plimy = arr{3};
xlab = arr{4};
clr = get(groot,'DefaultAxesColorOrder');
count = 1;
mn = nanmin(nanmin(fd_dk.(parm),[],2),[],1)*1.1;
mx = nanmax(nanmax(fd_dk.(parm),[],2),[],1)*1.1;
ndiffs = sum(1:(length(fd_dk.dk)-1));
nrows = ceil(sqrt(ndiffs));
if length(fd_dk.dk)<=2
    ncols = 1;
else
    ncols = floor(sqrt(ndiffs));
end
fig = figure();
set(fig,'Position',[500 0 285*ncols 240*nrows]);
for i = 1:length(fd_dk.dk)
for j = (i+1):length(fd_dk.dk)
    if i ~= j
        subplot(nrows,ncols,count)
        V1 = fd_dk.(parm)(fd_dk.inda,i);
        V2 = fd_dk.(parm)(fd_dk.inda,j);
        plt = plot(V1,V2,'.');
        title(['DK ' num2str(floor(fd_dk.dk(i))) ' Vs. DK ' num2str(floor(fd_dk.dk(j)))])
        hold on
        grid on
        plt = plot([-1000 1000],[-1000 1000],'k--');
        xlim([mn mx])
        ylim([mn mx])
%         set(plt,'LineWidth',1.5)
%         xlim(plimx)
%         ylim(plimy)
%         set(gca,'XTick',plimx(1):abs(diff(plimx))/4:plimx(2))
        axis square
        
        
        if count == 2*nrows-1 | count == 2*nrows
            xlabel(xlab)
        end
        
        count=count+1;
    end
end
end
saveas(gcf,[figdir 'scatter_' parm '.png'])


%%%%%
function make_consist_plots()
a = 1;


%%%%%
function fd_dk = make_dk_struct(fd_cut,param_array,dk,spec_ind)
global p p_ind;




if ~exist('spec_ind','var')% | isempty(spec_ind)
    spec_ind = true(size(fd_cut.ch,1),length(dk));
end
%dk = unique(fd.dk);

%dk = dk([1 2 4]);
fd_dk = [];
[phid,eps] = deal(NaN(size(p.gcp,1),length(dk)));
for dkind = 1:length(dk)
    
    % Collect the channels per DK.
    for chind = 1:length(p.gcp)
        ind = find(fd_cut.dk==dk(dkind) & fd_cut.ch==chind & spec_ind(:,dkind));
        if ~isempty(ind)
            phid(chind,dkind) = nanmean(fd_cut.phi_d(ind));
            eps(chind,dkind) = nanmean(fd_cut.aparam(ind,2));
        end
    end
    
end

fd_dk.phid = phid;
fd_dk.eps = eps;
fd_dk.dk = dk';
% Q / U are just the weighted means of the cos / sin of the angles
% A/B pols are weighted by one.
% Cross pol contribution is weighted by xpol leakage.

inda = p_ind.rgl100a;
indb = p_ind.rgl100b;

dkind = 1:length(dk);
phia = fd_dk.phid(inda,dkind);
phib = fd_dk.phid(indb,dkind);
epsa = fd_dk.eps(inda,dkind);
epsb = fd_dk.eps(indb,dkind);

Q = (cosd(2*phia)-cosd(2*phib)-epsa.*cosd(2*phia)+epsb.*cosd(2*phib))./(2+epsa+epsb);
U = (sind(2*phia)-sind(2*phib)-epsa.*sind(2*phia)+epsb.*sind(2*phib))./(2+epsa+epsb);
phip = atan2(U,Q)/2*180/pi;
xpol = 1-sqrt(Q.^2+U.^2);

[pltphi, plteps] = deal(NaN(size(phid)));
pltphi(inda,:) = atand(tand(phip-repmat((p.mce(inda)==0)*90,1,length(dk))));
pltphi(indb,:) = 0;
plteps(inda,:) = xpol;
plteps(indb,:) = 0;


fd_dk.phip = NaN(size(phid));
fd_dk.phip(inda,dkind) = phip;
fd_dk.phip(indb,dkind) = 0;
fd_dk.phip_corr(:,dkind) = pltphi;
fd_dk.xpol(:,dkind) = plteps;
fd_dk.inda = inda;
fd_dk.indb = indb;

function endrow = get_endrows(ndiffs,nrows,ncols)
V = get_rcvec(nrows,ncols);

Vind = V>ndiffs;
endrow = ones(1,ncols)*nrows;
for colind=1:ncols
    if Vind(nrows,colind)==1
        endrow(colind)=nrows-1;
    end
end



function V = get_rcvec(nrows,ncols)
    V = reshape(1:nrows*ncols,ncols,nrows)';

function i = rowcol2ind(r,c,nrows,ncols)
    V = reshape(1:nrows*ncols,ncols,nrows)';
    i = V(r,c);


function [r,c] = ind2rowcol(ind,nrows,ncols)
    V = reshape(1:nrows*ncols,ncols,nrows)';
    [r,c] = find(V==ind);
    
function index = inrange(A,B,C)

index = (A >= B) & (A <= C);

