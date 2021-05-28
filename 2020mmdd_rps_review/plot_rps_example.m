function plot_rps_example()

%% Init things
% Load and map tods
load(fullfile('b3_inst_paper_data','rps_ideal_tod.mat'))
load(fullfile('b3_inst_paper_data','rps_fitdata_20201013.mat'))

% Add pipeline paths
%addpath(fullfile('Z:','pipeline','beammap'))
%addpath(fullfile('Z:','pipeline','util'))

% Init misc variables
clr = get(groot,'DefaultAxesColorOrder'); %get the default colors
mrk = {'o','x'}; % Change the marker per plot
ch = [696,697]; % Known good channels
pols = ['A','B']; % Polarization
scaling = 20; % Scaling for the beam maps
limx = [-195 195]; % xlim range
limy = [-0.1 1.3]; % ylim range

dpix = 0.1; %size of the maps 
ybin = -20:dpix:20;
xbin = limx(1):dpix:limx(2);

%% Raster plots
fig = figure(2);
clf;
set(fig,'Position',[0*2200,100,700,500])

[X,Y] = meshgrid(xbin,ybin);

% Make the maps
for chi = 1:2
    subplot(4,1,chi)
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
    imagesc('XData',xbin,'YData',ybin,'CData',map,[0,1])
    grid on
    xlim(limx)
    
    a = gca;
    a.XTick = (-180:30:180)+15;
    a.XTickLabel = [];
    a.YTick = [];
    a.GridAlpha = 1;
    title([pols(chi) ' Detector Rasterset Data'])
    axis image
end

% Plot the rasters
clr1 = clr([1,3],:);
subplot(4,1,3:4)
plot(-100,-100,'Color',clr1(1,:));
hold on
plot(-100,-100,'Color',clr1(2,:));
for chi = 1:2
    chind = find(fd.sch==1 & fd.ch==ch(chi) & fd.row==10);
    modcurve = fd.bparam(chind,6:18);
    plot(fd.rot(chind,:),modcurve/max(modcurve),['k' mrk{chi}])
    grid on
    hold on
    rot = -180:180;
    plot(rot,rps_get_mod_model([fd.phi_d(chind)+90,0,0,0,1/2],rot),'Color',clr1(chi,:))
    
end
xlim(limx)
ylim(limy)
legend('A Detector','B Detector')
xlabel('Source Angle (^\circ)','FontSize',12)
ylabel('Source Amplitude (Normalized)','FontSize',10)
set(gca, 'XTick', [-180:60:180]);
%saveas(fig,[figdir 'modcurve.png'])

