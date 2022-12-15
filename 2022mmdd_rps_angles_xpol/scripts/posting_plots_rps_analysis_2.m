function posting_plots_rps_analysis_2()

clc; clear all;
load('z:/dev/rps/rps_beam_fits_type5_notilt_cut.mat')
fd_notilt = fd;
%load('z:/dev/rps/rps_beam_fits_type11_rerun_cut.mat')
load('z:/dev/rps/rps_beam_fits_type5_rerun_mirrorperrast_cut.mat')
fd_mpr = fd;
load('z:/dev/rps/rps_beam_fits_type5_rerun_cut.mat')
fd_type5 = fd;


%load('z:/dev/rps/rps_beam_fits_rerun_all_cut.mat')
figdir = 'c:/Users/James/Documents/GitHub/postings/2022mmdd_rps_pointing/figs/';
load('z:/dev/rps/fpu_data_obs.mat')
load('z:/dev/rps/pm.mat')
load('z:/dev/rps/source_fit_data.mat')
inrange = @(A,B,C) B <= A & A <= C;
outrange = @(A,B,C) A <= B | C <= A;
addpath('z:/dev')
addpath('z:/pipeline/beammap')
addpath('z:/pipeline/util')
addpath('z:/dev/rps')


prx = 2*sind(p.r/2).*cosd(p.theta)*180/pi;
pry = 2*sind(p.r/2).*sind(p.theta)*180/pi;

clc
% Moon-derived mirror params -- Use only these and nothing else!
mirror = struct();
mirror.height = 1.4592;
mirror.tilt= 44.8870;
mirror.roll = -0.075;
rpsopt.mirror = mirror;

rpsopt.source.distance = 195.5;
% Fit for the source params given our mirror info:
source = rps_fit_source(fd,rpsopt,p,'');
rpsopt.source = source;


% With new mirror and source parameters, update the pointing.
[fd.x,fd.y,phi] = beam_map_pointing_model(fd.az_cen,fd.el_cen,fd.dk_cen,model,'bicep3',mirror,source,[]);
fd.x = reshape(fd.x,size(fd.ch));
fd.y = reshape(fd.y,size(fd.ch));
fd.resx = reshape(prx(fd.ch),size(fd.ch))-fd.x;
fd.resy = reshape(pry(fd.ch),size(fd.ch))-fd.y;
[resth, resr] = cart2pol(fd.resx,fd.resy);
[fd.resx_rot, fd.resy_rot] = pol2cart(resth-fd.dk_cen*pi/180,resr);

% Make a median subtraction array.
fd.phi_medsub = fd.phi;
for chind = 1:2640
    ind = fd.ch==chind;
    if ~isempty(find(ind))
        fd.phi_medsub(ind) = fd.phi_medsub(ind)-nanmedian(fd.phi(ind));
    end
end

%%
fd_type5 = fd;
load('z:/dev/rps/rps_beam_fits_type5_rerun_mirrorperrast.mat')
fd_mpr = fd;
cutind = false(size(fdrast.ch));
for fdind = 1:length(fdtype5.ch)
    ind = fdrast.schnum == fdtype5.schnum(fdind) & fdrast.rowind == fd_type5.rowind(fdind) & fdrast.ch == fdtype5.ch(fdind);
    
    cutind(ind) = true;
end

fd_mpr = structcut(fdrast,cutind);
fd = fd_mpr;
save('z:/dev/rps/rps_beam_fits_type5_rerun_mirrorperrast_cut.mat','fd','scheds','titles','dks')

%%
% Order the a/b phis based on schedule
clc
[phis, xpols, phis_err, xs, ys, phi_pair, poleff_pair, phi_pair_err] = sort_fd_by_scheds(fd_type5,scheds,p,p_ind);

[phis_mpr, xpols_mpr, phis_err_mpr, xs_mpr, ys_mpr, phi_pair_mpr, poleff_pair_mpr, phi_pair_err_mpr] = sort_fd_by_scheds(fd_mpr,scheds,p,p_ind);

[phis_notilt, xpols_notilt, phis_err_notilt, xs_notilt, ys_notilt, phi_pair_notilt, poleff_pair_notilt, phi_pair_err_notilt] = sort_fd_by_scheds(fd_notilt,scheds,p,p_ind);

% moar cuts
cutind = true(size(phi_pair)) & abs(phi_pair)>5 & abs(phi_pair_notilt)>5 & abs(phi_pair_mpr)>5 ;
phi_pair(cutind) = NaN;
phi_pair_notilt(cutind) = NaN;
phi_pair_mpr(cutind) = NaN;


%% Compare pol angles with/without tiltmeter correction

clc
cm = colormap('lines');
V1 = phi_pair;
V2 = phi_pair_notilt;
ind = isnan(V1+V2);
V1(ind) = NaN;
V2(ind) = NaN;

%
edges = (-1:0.05/2:1);

fig = figure(1);
fig.Position(3:4) = [750 400];
clf;
t = tiledlayout(1,2);
nexttile(1)


phinm = nanmedian(V2,1);
ind = phinm==0 | sum(~isnan(V2),1)<5;
phinm(ind) = NaN;
v2 = reshape(V2-phinm,[],1);
v2(v2==0) = NaN;

N = histc(v2,edges);
b = bar(edges,N,'histc');
b.FaceColor = cm(1,:);
M = nanmean(v2);
S = nanstd(v2);
L = length(find(~isnan(v2)));
title({'No Tilt Meter Correction',...
    sprintf('M: %0.3f, STD: %0.3f, N: %0i',M,S,L)})
xlim([-1 1]*0.6)
ylim([-0.1 1550])
xlabel('\phi_{pair} - median(\phi_{pair}) [Deg]')
ylabel('N');
grid on
pbaspect([1 1 1])

nexttile(2)
phinm = nanmedian(V1,1);
ind = phinm==0 | sum(~isnan(V1),1)<5;
phinm(ind) = NaN;
v1 = reshape(V1-phinm,[],1);
v1(v1==0) = NaN;
N = histc(v1,edges);
b = bar(edges,N,'histc');
b.FaceColor = cm(1,:);
M = nanmean(v1);
S = nanstd(v1);
L = length(find(~isnan(v1)));
title({'With Tilt Meter Correction',...
    sprintf('M: %0.3f, STD: %0.3f, N: %0i',M,S,L)})
pbaspect([1 1 1])
xlim([-1 1]*0.6)
ylim([-0.1 1550])
ax = gca;
ax.YTickLabel = {};
xlabel('\phi_{pair} - median(\phi_{pair}) [Deg]')
grid on
t.TileSpacing = 'none';
t.Padding = 'loose';

%saveas(fig,'C:\Users\James\Documents\GitHub\postings\2022_spie_plots\figs\tilt_corr_hist.png')
%exportgraphics(fig,'C:\Users\James\Documents\GitHub\postings\2022_spie_plots\figs\tilt_corr_hist.pdf','ContentType','vector')

%% Compare pol angles with/without perraster mirror correction

clc
cm = colormap('lines');
V1 = phi_pair;
V2 = phi_pair_mpr;
ind = isnan(V1+V2);
V1(ind) = NaN;
V2(ind) = NaN;

%
edges = (-1:0.05/2:1);

fig = figure(2);
fig.Position(3:4) = [750 400];
clf;
t = tiledlayout(1,2);
nexttile(2)


phinm = nanmedian(V2,1);
ind = phinm==0 | sum(~isnan(V2),1)<5;
phinm(ind) = NaN;
v2 = reshape(V2-phinm,[],1);
v2(v2==0) = NaN;

N = histc(v2,edges);
b = bar(edges,N,'histc');
b.FaceColor = cm(1,:);
M = nanmean(v2);
S = nanstd(v2);
L = length(find(~isnan(v2)));
title({'Per-Rasterset Mirror Fit',...
    sprintf('M: %0.3f, STD: %0.3f, N: %0i',M,S,L)})
xlim([-1 1]*0.6)
ylim([-0.1 1550])
ax = gca;
ax.YTickLabel = {};
xlabel('\phi_{pair} - median(\phi_{pair}) [Deg]')

grid on
pbaspect([1 1 1])

nexttile(1)
phinm = nanmedian(V1,1);
ind = phinm==0 | sum(~isnan(V1),1)<5;
phinm(ind) = NaN;
v1 = reshape(V1-phinm,[],1);
v1(v1==0) = NaN;
N = histc(v1,edges);
b = bar(edges,N,'histc');
b.FaceColor = cm(1,:);
M = nanmean(v1);
S = nanstd(v1);
L = length(find(~isnan(v1)));
title({'Overall Mirror Fit',...
    sprintf('M: %0.3f, STD: %0.3f, N: %0i',M,S,L)})
pbaspect([1 1 1])
xlim([-1 1]*0.6)
ylim([-0.1 1550])
ylabel('N');
xlabel('\phi_{pair} - median(\phi_{pair}) [Deg]')
grid on
t.TileSpacing = 'none';
t.Padding = 'loose';

%saveas(fig,'C:\Users\James\Documents\GitHub\postings\2022_spie_plots\figs\tilt_corr_hist.png')
%exportgraphics(fig,'C:\Users\James\Documents\GitHub\postings\2022_spie_plots\figs\tilt_corr_hist.pdf','ContentType','vector')

%% Plot these against each other

cm = colormap('lines');
%V1 = reshape(phi_pair,[],1);
%V2 = reshape(phi_pair_mpr,[],1);
V1 = nanmean(phi_pair,1);
V2 = nanmean(phi_pair_mpr,1);
ind = isnan(V1+V2);
V1(ind) = NaN;
V2(ind) = NaN;

%
edges = (-1:0.05/2:1);

fig = figure(3);
fig.Position(3:4) = [750 400];
clf;
t = tiledlayout(1,2);
nexttile(1)
scatter(V1,V2,16,cm(1,:),'filled')
grid on
pbaspect([1 1 1])

nexttile(2)
edges = (-1:0.05:1)*0.1;
N = histc(V1-V2,edges);
bar(edges,N,'histc')
grid on
xlim([min(edges) max(edges)])
pbaspect([1 1 1])

%% Plot as a function of DK
clc
ind0 = inrange(nanmean(phis,1),-5,5);
ind90 = inrange(nanmean(phis,1),-5+90,5+90);
cm = colormap('lines');
V1 = phis;
V2 = phis_mpr;
DK = repmat(dks',1,2640);
ind = isnan(V1+V2);
V1(ind) = NaN;
V2(ind) = NaN;

%
edges = (-1:0.05/2:1);

fig = figure(4);
fig.Position(3:4) = [750 400];
clf; hold on;
t = tiledlayout(2,1);
nexttile(1)
hold on;
scatter(DK(:,ind0),V1(:,ind0),16,cm(1,:),'filled')
scatter(DK(:,ind90),V1(:,ind90)-90,16,cm(2,:),'filled')
grid on

nexttile(2)
hold on;
scatter(DK(:,ind0),V2(:,ind0),16,cm(1,:),'filled')
scatter(DK(:,ind90),V2(:,ind90)-90,16,cm(2,:),'filled')
grid on

%pbaspect([1 1 1])



%%

function [phis, xpols, phis_err, xs, ys, phi_pair, poleff_pair, phi_pair_err] = sort_fd_by_scheds(fd,scheds,p,p_ind)

len = length(scheds);

inda = p_ind.a;
indb = p_ind.b;

[phis, xpols, phis_err, xs, ys] = deal(NaN(len,2640));
for schind = 1:len
    for chind = 1:2640
        ci = find(fd.ch==chind & ismember(fd.schnum,scheds{schind}));

        if ~isempty(ci)

            ch_chi = atand(tand(p.theta(chind)+p.chi(chind)));

            val = reshape(fd.phi(ci),[],1);
            xp = reshape(fd.xpol(ci),[],1);
            err = reshape(fd.phi_err(ci),[],1);

            %err(isnan(val))=1e10;
            xs(schind,chind) = mean(fd.x(ci));
            ys(schind,chind) = mean(fd.y(ci));
            phis(schind,chind) = mean(val);%wmean(val,1./err,1);
            xpols(schind,chind) = mean(xp);%wmean(xp,1./err,1);
            phis_err(schind,chind) = nanmin(err);

        end
    end
end


% Calculate pair-diff angles
% Loop over channels to account for MCE0;
[phi_pair, poleff_pair, phi_pair_err] = deal(NaN(len,2640));
for schind = 1:len
    for chind = 1:length(inda)

        pha = phis(schind,inda(chind));
        phb = phis(schind,indb(chind));
        ea = xpols(schind,inda(chind));
        eb = xpols(schind,indb(chind));

        phi_pair_err(schind,inda(chind)) = max(phis_err(schind,[inda(chind) indb(chind)]));
        if p.mce(inda(chind))~=0

            [phi_pair(schind,inda(chind)), poleff_pair(schind,inda(chind))] = calc_pair_diff_pol(pha,phb,ea,eb);
        else

            [phi_pair(schind,inda(chind)), poleff_pair(schind,inda(chind))] = calc_pair_diff_pol(phb,pha,eb,ea);
        end

        if 0%~inrange(phi_pair(schind,inda(chind))+2.5,-2.5,2.5)
            phi_pair(schind,inda(chind)) = NaN;
            poleff_pair(schind,inda(chind)) = NaN;
            phi_pair_err(schind,inda(chind)) = 1e10;
        end
    end
    %ind = abs(phi_pair(schind,:)-nanmean(phi_pair(schind,:)))<1;
end



