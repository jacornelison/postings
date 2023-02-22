function posting_plots_mirror_diff_pol()

%%
set(groot,'defaulttextinterpreter','latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(groot,'defaultAxesFontSize',12)


%%

clc
clear all
close all
postdir = fullfile('C:','Users','James','Documents','GitHub','postings','20230220_diff_reflection_polarization','');
datadir = fullfile(postdir,'data','');
figdir = fullfile(postdir,'figs','');
%load(fullfile(datadir,'b3rpsfiles.mat'))
load('z:/dev/rps/fpu_data_obs.mat')
addpath('z:/dev/')
addpath('z:/dev/rps/')
addpath('z:/dev/diff_polarization/')
n = 2200;
k = 2200;
dk = -180:45:180;

[phipair,phia,phib,Ainc,epspair] = calc_mirror_diff_pol(p,n,k,dk);

fig = figure(1);
fig.Position(3:4) = [900 700];

%%
scaling = 0.000002;
names = {'phipair','phia','phib','epspair'};
vals = {phipair,phia,phib-90,epspair};
labels = {'\phi_{pair} [Deg]','\phi_{a} [Deg]','\phi_{b}-90 [Deg]','\epsilon_{pair}'};
scales = [1e-10,4e-2,4e-2,4e-7];

for valind = 1:length(vals)

for dkind = 1:length(dk)
    figure(1); clf;
    plot_tiles(vals{valind}(:,dkind),p,'fig',1,'clim',[-1,1]*scales(valind),'title',['dk ' num2str(dk(dkind))],'pair','mean','clab',labels{valind});
    %plot_tiles(vals{valind}(:,dkind),p,'fig',1,'title',['dk ' num2str(dk(dkind))],'pair','mean','clab',labels{valind});
    colormap jet
    fname = sprintf('tileplot_%s_dk_%i',names{valind},dk(dkind));
    fname = fullfile(figdir,fname);
    saveas(1,fname,'png')
end


plot_tiles(mean(vals{valind},2),p,'fig',1,'clim',[-1,1]*scales(valind),'title','mean','pair','mean','clab',labels{valind});
    colormap jet
    fname = sprintf('tileplot_%s_dk_mean',names{valind});
    fname = fullfile(figdir,fname);
    saveas(1,fname,'png')

    
    plot_tiles(mean(vals{valind}(:,1:4),2),p,'fig',1,'clim',[-1,1]*scales(valind),'title','mean 0-135','pair','mean','clab',labels{valind});
    colormap jet
    fname = sprintf('tileplot_%s_dk_mean_rps',names{valind});
    fname = fullfile(figdir,fname);
    saveas(1,fname,'png')

end
%% Scatter plot of phia/b/pair vs DK
clc
fig = figure(1);
fig.Position(3:4) = [900 800];
clf;

p_dummy = struct;
p_dummy.r = [0; 7; 15];
p_dummy.theta = zeros(size(p_dummy.r));
p_dummy.drumangle = zeros(size(p_dummy.r));
p_dummy.chi = zeros(size(p_dummy.r));
p_dummy.chi_thetaref=zeros(size(p_dummy.r));


[phipair,phia,phib,Ainc,epspair] = calc_mirror_diff_pol(p_dummy,n,k,dk);

V0 = {phia, phib-90, phipair};
ylims = {[-1 1]*1, [-1 1]*1, [-1 1]*0.5};
yttls = {'\phi_0','\phi_{90}','\phi_{pair}'};
for valind = 1:length(V0)
    K = repmat(dk,length(p_dummy.r),1);
    K = reshape(K,[],1);
    clr = repmat(p_dummy.r,1,length(dk));
    clr = reshape(clr,[],1);
    V = reshape(V0{valind},[],1);
    subplot(3,1,valind)
    scatter(K,V,14,clr,'filled')
    grid on
    colormap('turbo')
    ylabel([yttls{valind} ' [Degrees]'])
    %ylim(ylims{valind})

end
xlabel('DK angle [Degrees]')

%% Plot the reflection coefficients

clc
fig = figure(11);
fig.Position(3:4) = [450 400];
clf; hold on;

Nair = 1;
Nal = 2200*(1+sqrt(-1));

theta_i = 0:65;
Atx = sqrt(1-(Nair/Nal*sind(theta_i)).^2);
Rs = abs((Nair*cosd(theta_i)-Nal*Atx)./(Nair*cosd(theta_i)+Nal*Atx)).^2;
Rp = abs((Nair*Atx-Nal*cosd(theta_i))./(Nair*Atx+Nal*cosd(theta_i))).^2;

plot(theta_i,Rs,'b')
plot(theta_i,Rp,'r')
grid on
title('Reflectivity of Aluminum @ 90GHz, 240 Kelvin')
xlabel('\theta_{inc} [Degrees]')
ylabel('Reflectivity')
legend({'S (Perpendicular)','P (Parallel)'},'location','southwest')
fname = 'al_reflectivity.png';
saveas(fig,fullfile(figdir,fname),'png')

%% Plot mirror diff pol
fig = figure(1);
fig.Position(3:4) = [900, 700*2/3];
clf; hold on;
set(groot,'defaultAxesFontSize',14)
clc
rs = [0 7 14];
DK = -180:5:180;
phi_rot = NaN(length(DK),2);
phi0 = {0 90};
ylabs = {'$\phi$','$\phi$','$\phi_{pair}$'};
ttls = {'Input of 0 degrees','Input of 90 Degrees','Pair-Diff Angle'};
for phiind = 1:length(phi0)
    subplot(2,1,phiind)
    hold on
for rind = rs
    if phiind ~=3
    phi_in = zeros(size(DK))+phi0{phiind};
    r = zeros(size(DK))+rind;
    theta = zeros(size(DK));
    mtilt = zeros(size(DK))+45;
    mroll = zeros(size(DK));
    
    V = mirror_diff_pol_calc(phi_in,r,theta,2200,2200,DK,mtilt,mroll);
    %V = mirror_diff_pol_correction(V,r,theta,2400,2400,DK,mtilt,mroll);
    phi_rot(:,phiind) = V;
    else
        V = (phi_rot(:,1)+phi_rot(:,2)-90)/2;
    end
    plot(DK,V)
    
end
legend({'r = 0deg' 'r = 7deg' 'r = 14deg'})
grid on
title(ttls{phiind})
    ylabel([ylabs{phiind} ' [Degrees]'])

end
xlabel('DK [Degrees]')
fname = 'phi_bias.png';
saveas(fig,fullfile(figdir,fname),'png')


%% Plot Correction
fig = figure(1);
fig.Position(3:4) = [900, 700*2/3];
clf; hold on;
set(groot,'defaultAxesFontSize',14)
clc
rs = [0 7 14];
DK = -180:5:180;
phi_rot = NaN(length(DK),2);
phi0 = {0 90};
ylabs = {'$\phi$','$\phi$','$\phi_{pair}$'};
ttls = {'Input of 0 degrees','Input of 90 Degrees','Pair-Diff Angle'};
for phiind = 1:length(phi0)
    subplot(2,1,phiind)
    hold on
for rind = rs
    if phiind ~=3
    phi_in = zeros(size(DK))+phi0{phiind};
    r = zeros(size(DK))+rind;
    theta = zeros(size(DK));
    mtilt = zeros(size(DK))+45;
    mroll = zeros(size(DK));
    
    V = mirror_diff_pol_calc(phi_in,r,theta,2200,2200,DK,mtilt,mroll);
    V = mirror_diff_pol_correction(V,r,theta,2200,2200,DK,mtilt,mroll);
    phi_rot(:,phiind) = V;
    else
        V = (phi_rot(:,1)+phi_rot(:,2)-90)/2;
    end
    plot(DK,V)
    
end
legend({'r = 0deg' 'r = 7deg' 'r = 14deg'})
grid on
title(ttls{phiind})
    ylabel([ylabs{phiind} ' [Degrees]'])

end
xlabel('DK [Degrees]')
fname = 'phi_correction.png';
saveas(fig,fullfile(figdir,fname),'png')



%% Check for bias in fits
fig = figure(1);
fig.Position(3:4) = [900, 700*2/3];
clf; hold on;
set(groot,'defaultAxesFontSize',14)

clc
DK = -180:5:180;
parmout_a = NaN(size(DK,1),5);
for dkind = 1:length(DK)
rot = -180:30:180;
Z = zeros(size(rot));
r = Z+14;
theta = Z;
k = Z-DK(dkind);
mtilt = Z+45;
mroll = Z+0;
V = mirror_diff_pol_calc(rot,r,theta,2200,2200,k,mtilt,mroll);
V = reshape(V,size(rot));
rot_biased = V;

parmin = [0, 0, 0, 0, 1];
A = rps_get_mod_model(parmin,rot_biased);

mxfev = 100000;
mxiter = 100000;
options = optimset('TolFun',1e-10,'MaxIter',mxiter,'MaxFunEvals',mxfev,'Display','off');
lb = [-20 -0.5 -10 -10 1e6];
ub = [20 0.5 10 10 1e6];
guess = [0 0 0 0 1];
parmout_a(dkind,:) = lsqcurvefit(@rps_get_mod_model,guess,rot,A,lb,ub,options);

end

fig = figure(1);
clf;

plot(DK,parmout_a(:,1))

