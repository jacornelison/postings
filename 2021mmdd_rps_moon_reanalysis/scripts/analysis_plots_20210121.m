function analysis_plots_20210121()


%%
clr = {[0, 0.4470, 0.7410],...
    [0.8500, 0.3250, 0.0980],...
    [0.4660, 0.6740, 0.1880]*0.8,...
    [0.9290, 0.6940, 0.1250]};

load('../data/fit_data_temp.mat')
load('../data/fit_data_final.mat')
load('../data/cmb_derived_centers.mat')

close all


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
%%
% Look at residuals
scale = 10;
fig = figure(3);
fig.Position = [700 400 500 500];
clf; hold on;

pxy = reshape(fd.data,[],2);
%count = 0;
pols = {'rgl100a','rgl100b'};
clr = {'b','r'};
for polind = 1:2
    chind = ismember(fd.fit_ch,p_ind.(pols{polind}));
    
    %plot((1:length(chind))+count,fd.fit_ch(chind),'color',clr{schind})
    %count = count+length(chind);
    quiver(pxy(chind,1),pxy(chind,2),fd.resx(chind)*scale,fd.resy(chind)*scale,0,'Color',clr{polind})
    %plot(res(chind,1),res(chind,2),'.','color',clr{schind})
end
grid on
%legend({'DK=0','DK=45','DK=90'})
legend({'Pol A','Pol B'})
title('Beam Center Best-Fit Residuals x10')
xlabel('x (^o)')
ylabel('y (^o)')
figname = ['../figs/' 'moonfit_residuals_quiver_zoomed'];
saveas(gcf,figname,'png')



%%
lims = [0, 15];
fig = figure(4);
fig.Position = [700 400 500 500];
clf; hold on;
scaling = 0.15;
h = scatterhist(fd.resx,fd.resy,'group',p.pol(fd.fit_ch),...
    'Kernel','on','location','southeast','Color','br');

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

figname = ['../figs/' 'moonfit_residuals_scatterhist'];
saveas(gcf,figname,'png')



