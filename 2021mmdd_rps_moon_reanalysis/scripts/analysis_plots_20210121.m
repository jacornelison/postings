function analysis_plots_20210121()


%%
clr = {[0, 0.4470, 0.7410],...
    [0.8500, 0.3250, 0.0980],...
    [0.4660, 0.6740, 0.1880]*0.8,...
    [0.9290, 0.6940, 0.1250]};

%load('../data/fit_data_temp.mat')
load('../data/fit_data_final.mat')
load('../data/cmb_derived_centers.mat')
%load('../data/moon_beamfit_model.mat')
close all

%res_beam = reshape(fd.data-model,[],2);
%resx = res_beam(:,1);
%resy = res_beam(:,2);

%%
% Moon Obs Coverage

figure()
set(gcf,'Position',[1000 100 600 600])
clf; hold on;
for i = 1:3
    plot(md{i}.x,md{i}.y,'Color',clr{i})
end
plot(brx,bry,'kx')
grid on
xlabel('x (^o)')
ylabel('y (^o)')
title('Moon Observation Coverage Jan 2017')
legend({'DK = 0','DK = 45','DK = 90','CMB-derived pnts'})
axis square


% Looking at a specific detector
figure();
set(gcf,'Position',[100,400,1200,400])
clf; hold on;
for schind = 1:3
    subplot(1,3,schind)
    chind = md{schind}.ch==696;%fit_ch(690);
    if ~isempty(find(chind))
        px = prx(md{schind}.ch(chind));
        py = pry(md{schind}.ch(chind));
        x = md{schind}.x;
        y = md{schind}.y;
        az = md{schind}.az;
        el = md{schind}.el;
        maskind = sqrt((x-px).^2+(y-py).^2) < 1;
        plot3(x(maskind),y(maskind),fbstruct{schind}.fb(maskind,chind),'Color',clr{schind})
        %plot3(el(maskind),az(maskind),fbstruct{schind}.fb(maskind,chind),'Color',clr{schind})
        grid on
        xlabel('x (^o)')
        ylabel('y (^o)')
    end
end

%% Look at residuals per Pol

plts = {{'rgl100a'};{'rgl100b'};{'rgl100a','rgl100b'}};
pltlabs = {'a','b','both'};
for pltind = 1:length(plts)
    scale = 10;
    fig = figure(3);
    fig.Position = [700 400 500 500];
    clf; hold on;
    
    pxy = reshape(fd.data,[],2);
    
    pols = plts{pltind};
    clr = {'b','r'};
    for polind = 1:length(pols)
        chind = ismember(fd.fit_ch,p_ind.(pols{polind}));
        quiver(pxy(chind,1),pxy(chind,2),fd.resx(chind)*scale,fd.resy(chind)*scale,0,'Color',clr{polind})
        
    end
    grid on
    %legend({'DK=0','DK=45','DK=90'})
    legend(plts{pltind})
    title('Beam Center Best-Fit Residuals x10')
    xlabel('x (^o)')
    ylabel('y (^o)')
    figname = ['../figs/' 'moonfit_residuals_quiver_zoomed_pol' pltlabs{pltind}];
    saveas(gcf,figname,'png')
end

%% Look at residuals per DK

clr = {[0, 0.4470, 0.7410],...
    [0.8500, 0.3250, 0.0980],...
    [0.4660, 0.6740, 0.1880]*0.8,...
    [0.9290, 0.6940, 0.1250]};


plts = {[1];[2];[3];,[1,2,3]};
pltlabs = {'0','45','90','all'};
legs = {{'0'},{'45'},{'90'},{'0','45','90'}};
for pltind = 1:length(plts)
    scale = 10;
    fig = figure(3);
    fig.Position = [700 400 500 500];
    clf; hold on;
    
    pxy = reshape(fd.data,[],2);
    
    dks = plts{pltind};
    
    for dkind = 1:length(dks)
        chind = ismember(fd.sch,dks(dkind)) & sqrt(fd.resx.^2+fd.resy.^2)'<1;
        quiver(pxy(chind,1),pxy(chind,2),fd.resx(chind)*scale,fd.resy(chind)*scale,0,'Color',clr{dkind})
        
    end
    grid on
    %legend({'DK=0','DK=45','DK=90'})
    legend(legs{pltind})
    title('Beam Center Best-Fit Residuals x10')
    xlabel('x (^o)')
    ylabel('y (^o)')
    xlim([-15 15])
    ylim([-15 15])
    figname = ['../figs/' 'moonfit_residuals_quiver_zoomed_dk' pltlabs{pltind}];
    saveas(gcf,figname,'png')
end


%% Look at Tile 11
clr = {[0, 0.4470, 0.7410],...
    [0.8500, 0.3250, 0.0980],...
    [0.4660, 0.6740, 0.1880]*0.8,...
    [0.9290, 0.6940, 0.1250]};

plts = {[1];[2];[3];,[1,2,3]};
pltlabs = {'0','45','90','all'};
legs = {{'0'},{'45'},{'90'},{'0','45','90'}};
for pltind = 1:length(plts)
    scale = 5;
    fig = figure(3);
    fig.Position = [700 400 500 500];
    clf; hold on;
    
    pxy = reshape(fd.data,[],2);
    
    dks = plts{pltind};
    
    for dkind = 1:length(dks)
        chind = ismember(fd.sch,dks(dkind)) & p.tile(fd.fit_ch)'==11 & sqrt(fd.resx.^2+fd.resy.^2)'<1;
        quiver(pxy(chind,1),pxy(chind,2),fd.resx(chind)*scale,fd.resy(chind)*scale,0,'Color',clr{dkind})
        
    end
    grid on
    legend(legs{pltind})
    title('Beam Center Best-Fit Residuals x5')
    xlabel('x (^o)')
    ylabel('y (^o)')
    xlim([-2.5 2.5])
    ylim([-2.5 2.5])
    figname = ['../figs/' 'moonfit_residuals_quiver_zoomed_dk' pltlabs{pltind} '_tile11'];
    saveas(gcf,figname,'png')
end







%% Scatter plots with marginal histograms per DK.
% DOES NOT WORK WITH MATAB 2009!


if str2num(datestr(version('-date'),'yyyy'))>2017
    
    plts = {[1];[2];[3];,[1,2,3]};
    pltlabs = {'0','45','90','all'};
    
    for pltind = 1:length(plts)
        chind = ismember(fd.sch,plts{pltind}) & sqrt(fd.resx.^2+fd.resy.^2)'<0.15;
        lims = [0, 15];
        fig = figure(4);
        fig.Position = [700 400 500 500];
        clf; hold on;
        scaling = 0.15;
        h = scatterhist(fd.resx(chind),fd.resy(chind),'group',p.pol(fd.fit_ch(chind)),...
            'Kernel','off','location','southeast','Color','br',...
            'Marker','xo','MarkerSize',[5,3]);
        
        legend({'Pol A','Pol B'})
        title('Beam Center Best-Fit Residuals')
        xlabel('\Delta x (^o)')
        ylabel('\Delta y (^o)')
        xlim([-1 1]*scaling)
        ylim([-1 1]*scaling)
        grid on
        
        pause(0.1)
        
        set(h(3),...
            'Position',h(3).Position-[0.0 0 0 0],...
            'XLim',h(3).XLim,...
            'YLim',lims);
        
        ax3 = axes('Position', h(3).Position,...
            'Color', 'none',...'XColor', 'none',...
            'XAxisLocation', 'bottom',...
            'YAxisLocation', 'right',...
            'XDir','reverse',...
            'YTickLabel','',...
            'XLim', h(3).YLim,...
            'YLim', h(3).XLim);
        
        xlabel(ax3, 'N');
        grid on
        
        
        set(h(2),...
            'XLim',h(2).XLim,...
            'YLim',lims);
        ax2 = axes('Position', h(2).Position,...
            'Color', 'none',...
            'XAxisLocation', 'bottom',...
            'YAxisLocation', 'left',...
            'YDir','normal',...
            'XTickLabel','',...
            'XLim', h(2).XLim,...
            'YLim', h(2).YLim);
        ylabel(ax2, 'N');
        grid on
        
        figname = ['../figs/' 'moonfit_residuals_scatterhist_dk' pltlabs{pltind}];
        saveas(gcf,figname,'png')
        
    end
    
end

%% Look at RPS fits


%% RPS Obs Coverage

load('../data/rps_fit_data_temp.mat')
load('../data/b3rpsfiles_2017.mat')
prx = 2 * sind(p.r / 2) .* cosd(p.theta) * 180 / pi;
pry = 2 * sind(p.r / 2) .* sind(p.theta) * 180 / pi;

addpath('z:/pipeline/util')
map = colormap('lines');
figure(3)
set(gcf,'Position',[1000 100 600 600])
clf; hold on;
h = [];
for i = 1:length(md)
    [x,y,phi] = beam_map_pointing_model(md{i}.az,md{i}.el,md{i}.dk,...
        rpsopt.pm,'bicep3',rpsopt.mirror,rpsopt.source,[]);
    
    h(end+1) = plot(x,y,'Color',map(i,:));
end
h(end+1) = plot(prx,pry,'kx');
grid on
xlabel('x (^o)')
ylabel('y (^o)')
title('Moon Observation Coverage Jan 2017')
legend(h([1,end]),{'RPS Scan','CMB-derived pnts'})
axis square
figname = ['../figs/rps_coverage_2017'];
    saveas(gcf,figname,'png')


%% 2018 RPS Residuals
load('../data/rps_fit_data_final.mat')

clr = {[0, 0.4470, 0.7410],...
    [0.8500, 0.3250, 0.0980],...
    [0.4660, 0.6740, 0.1880]*0.8,...
    [0.9290, 0.6940, 0.1250]};

plts = {[1]};%;[2];[3];,[1,2,3]};
pltlabs = {'0'};%,'45','90','all'};
legs = {{'0'}};%,{'45'},{'90'},{'0','45','90'}};
for pltind = 1:length(plts)
    scale = 10;
    fig = figure(3);
    fig.Position = [700 400 500 500];
    clf; hold on;
    
    pxy = reshape(fd.data,[],2);
    
    dks = plts{pltind};
    
    for dkind = 1:length(dks)
        chind = true(size(fd.sch));%ismember(fd.sch,dks(dkind)) & sqrt(fd.resx.^2+fd.resy.^2)'<1;
        quiver(pxy(chind,1),pxy(chind,2),fd.resx(chind)*scale,fd.resy(chind)*scale,0,'Color',clr{dkind})
        
    end
    grid on
    %legend({'DK=0','DK=45','DK=90'})
    legend(legs{pltind})
    title(['Beam Center Best-Fit Residuals x' num2str(scale)])
    xlabel('x (^o)')
    ylabel('y (^o)')
    xlim([-15 15])
    ylim([-15 15])
    figname = ['../figs/' 'rpsfit_residuals_quiver_zoomed_dk' pltlabs{pltind}];
    saveas(gcf,figname,'png')
end

    
%% 2018 RPS Obs Coverage

load('../data/rps18_fit_data_temp.mat')
load('../data/b3rpsfiles_2018.mat')
rpsopt.pm = fd.rpsopt.pm;
prx = 2 * sind(p.r / 2) .* cosd(p.theta) * 180 / pi;
pry = 2 * sind(p.r / 2) .* sind(p.theta) * 180 / pi;

addpath('z:/pipeline/util')
map = colormap('lines');

plts = {[1];[2];[3];[4]};
pltlabs = {'0','45','90','135'};
dk = unique(fd.dk);

for pltind = 1:length(plts)
    dks = dk(plts{pltind});
    figure(3)
set(gcf,'Position',[1000 100 600 600])
clf; hold on;
    
h = [];
count = 1;
for i = 1:length(md)
    if md{i}.dk0 == dks;
    [x,y,phi] = beam_map_pointing_model(md{i}.az,md{i}.el,md{i}.dk,...
        rpsopt.pm,'bicep3',rpsopt.mirror,rpsopt.source,[]);
    
    h(end+1) = plot(x,y,'Color',map(count,:));
    if count == size(map,1)
        count = 1;
    else
        count = count+1;
    end
    
    end
end
h(end+1) = plot(prx,pry,'kx');
grid on
xlabel('x (^o)')
ylabel('y (^o)')
title('RPS Observation Coverage Jan 2018')
legend(h([1,end]),{'RPS Scan','CMB-derived pnts'})
axis square
xlim([-1 1]*20)
ylim([-1 1]*20)
figname = ['../figs/rps_coverage_2018_dk_' pltlabs{pltind}];
    saveas(gcf,figname,'png')

end

%% 2018 RPS residuals in map coords.
load('../data/rps18_fit_data_final.mat')

clr = {[0, 0.4470, 0.7410],...
    [0.8500, 0.3250, 0.0980],...
    [0.4660, 0.6740, 0.1880]*0.8,...
    [0.9290, 0.6940, 0.1250]};

plts = {[1];[2];[3];[4];[1,2,3,4]};
pltlabs = {'0','45','90','135','all'};
legs = {{'0'},{'45'},{'90'},{'135'},{'0','45','90','135'}};
dk = unique(fd.dk);
for pltind = 1:length(plts)
    scale = 10;
    fig = figure(3);
    fig.Position = [700 100 500 500];
    clf; hold on;
    
    pxy = reshape(fd.data,[],2);
    mod = reshape(fd.model,[],2);
    dks = dk(plts{pltind});
    
    for dkind = 1:length(dks)
        chind = ismember(fd.dk,dks(dkind));
        
        [x_mirr, y_mirr] = get_mirror_coords(fd.dk(chind),pxy(chind,1),pxy(chind,2),1.5,rpsopt.mount,rpsopt.mirror);
        [xmod_mirr, ymod_mirr] = get_mirror_coords(fd.dk(chind),mod(chind,1),mod(chind,2),1.5,rpsopt.mount,rpsopt.mirror);
        resxm = x_mirr-xmod_mirr;
        resym = y_mirr-ymod_mirr;
        quiver(x_mirr,y_mirr,resxm*scale,resym*scale,0,'Color',clr{dkind})
        
    end
    grid on
    %legend({'DK=0','DK=45','DK=90'})
    legend(legs{pltind})
    title({['2018 RPS Beam Center Best-Fit Residuals x' num2str(scale)],'Mirror Coordinates'})
    xlabel('x (m)')
    ylabel('y (m)')
    xlim([-1 1]*0.4)
    ylim([-1 1.70]*0.4)
    figname = ['../figs/' 'rpsfit18_residuals_quiver_mirrorcoords_dk' pltlabs{pltind}];
    saveas(gcf,figname,'png')
end

%% Look at residuals of residuals

fd1 = load('../data/fit_data_final.mat')
fd2 = load('../data/rps18_fit_data_final.mat')


fig = figure(3);
fig.Position = [700 100 500 500];
clf; hold on;
    


