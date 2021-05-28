function posting_plots_mirror_diff_pol()


load('../data/b3rpsfiles.mat')
n = 2400;
k = 2400;
dk = -180:45:135;

[phipair,phia,phib,Ainc] = calc_mirror_diff_pol(p,n,k,dk);

%%
scaling = 0.000002;
for dkind = 1:length(dk)
    figure(1); clf;
    plot_tiles(phipair(:,dkind),p,'fig',1,'clim',[-1,1]*scaling,'title',['dk ' num2str(dk(dkind))],'pair','mean','clab','\phi_{pair}');
    colormap jet
    fname = sprintf('../figs/tileplot_phipair_dk_%i',dk(dkind));
    saveas(1,fname,'png')
end

plot_tiles(mean(phipair,2),p,'fig',1,'clim',[-1,1]*scaling,'title',['mean'],'pair','mean','clab','\phi_{pair}');
    colormap jet
    fname = sprintf('../figs/tileplot_phipair_dk_mean');
    saveas(1,fname,'png')

    
    plot_tiles(mean(phipair(:,1:4),2),p,'fig',1,'clim',[-1,1]*scaling,'title',['mean 0-135'],'pair','mean','clab','\phi_{pair}');
    colormap jet
    fname = sprintf('../figs/tileplot_phipair_dk_mean_rps');
    saveas(1,fname,'png')


 %%
    scaling = 0.005;
for dkind = 1:length(dk)
    figure(1); clf;
    plot_tiles(phia(:,dkind),p,'fig',1,'clim',[-1,1]*scaling,'title',['dk ' num2str(dk(dkind))],'pair','mean','clab','\phi_a');
    colormap jet
    fname = sprintf('../figs/tileplot_phia_dk_%i',dk(dkind));
    saveas(1,fname,'png')
end

plot_tiles(mean(phia,2),p,'fig',1,'clim',[-1,1]*scaling,'title',['mean'],'pair','mean','clab','\phi_a');
    colormap jet
    fname = sprintf('../figs/tileplot_phia_dk_mean');
    saveas(1,fname,'png')

    
    plot_tiles(mean(phia(:,1:4),2),p,'fig',1,'clim',[-1,1]*scaling,'title',['mean 0-135'],'pair','mean','clab','\phi_a');
    colormap jet
    fname = sprintf('../figs/tileplot_phia_dk_mean_rps');
    saveas(1,fname,'png')

scaling = 0.005;
for dkind = 1:length(dk)
    figure(1); clf;
    plot_tiles(phib(:,dkind)-90,p,'fig',1,'clim',[-1,1]*scaling,'title',['dk ' num2str(dk(dkind))],'pair','mean','clab','\phi_b-90');
    colormap jet
    fname = sprintf('../figs/tileplot_phib_dk_%i',dk(dkind));
    saveas(1,fname,'png')
end

plot_tiles(mean(phib,2)-90,p,'fig',1,'clim',[-1,1]*scaling,'title',['mean'],'pair','mean','clab','\phi_b-90');
    colormap jet
    fname = sprintf('../figs/tileplot_phib_dk_mean');
    saveas(1,fname,'png')

    
    plot_tiles(mean(phib(:,1:4),2)-90,p,'fig',1,'clim',[-1,1]*scaling,'title',['mean 0-135'],'pair','mean','clab','\phi_b-90');
    colormap jet
    fname = sprintf('../figs/tileplot_phib_dk_mean_rps');
    saveas(1,fname,'png')
