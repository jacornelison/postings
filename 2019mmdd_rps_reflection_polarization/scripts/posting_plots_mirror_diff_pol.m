function posting_plots_mirror_diff_pol()


load('../data/b3rpsfiles.mat')
n = 2000;
k = 2000;
dk = -180:45:135;
scaling = 0.02;
[phipair,phia,phib,Ainc] = calc_mirror_diff_pol(p,n,k,dk);

for dkind = 1:length(dk)
    figure(1); clf;
    plot_tiles(phipair(:,dkind),p,'fig',1,'clim',[-1,1]*scaling,'title',['dk ' num2str(dk(dkind))],'pair','mean');
    colormap jet
    fname = sprintf('../figs/tileplot_phipair_dk_%i',dk(dkind));
    saveas(1,fname,'png')
end

plot_tiles(mean(phipair,2),p,'fig',1,'clim',[-1,1]*scaling,'title',['mean'],'pair','mean');
    colormap jet
    fname = sprintf('../figs/tileplot_phipair_dk_mean');
    saveas(1,fname,'png')

    
    plot_tiles(mean(phipair(:,1:4),2),p,'fig',1,'clim',[-1,1]*scaling,'title',['mean 0-135'],'pair','mean');
    colormap jet
    fname = sprintf('../figs/tileplot_phipair_dk_mean_rps');
    saveas(1,fname,'png')

    
