clc
datadir = 'postings/20191230_bicep3_opteff/data/';
m1 = load(sprintf('%s/B3_SP2018_20181217_run10_OE/inter_results.mat',datadir));
m2 = load(sprintf('%s/B3_SP2019_20191228_run11_OE/inter_results.mat',datadir));
mk = {[0,0,0.7],[0.7,0,0],[0,0.7,0.2],[0.5,0,0.5]};
p = get_array_info(20191228);
kboltz = 1.38e-23;
dnu = 25e9;

load('dev/paper_data/rps_workspace_for_plotting.mat')
fd_uniform = load('dev/rps2022_data_refit_uniform');
fd_uniform = fd_uniform.fd;
[fd_uniform, phis_uniform, phi_pair_uniform, xpols_uniform, poleffs_uniform,n1s_uniform,n2s_uniform,amps_uniform] = get_pol_params_per_obs(fd_uniform,p,scheds);


figdir = 'postings/2025mmdd_mod_curve_weighting/figs/';
%%

idx = ismember(p.tile,[3,4,5,7,18]);
fig = figure(2);
hold off
v1 = nanmedian(poleffs_uniform,1);
scatter(m2.PSat{2}(idx,1)*1e12,v1(idx),8,p.tile(idx))
xlim([50,250])
%xlim([87,110])
% scatter(m2.dPdT*1e12,nanmedian(poleffs_uniform,1),8,p.tile)
% xlim([0.05,0.15])
%ylim([-0.015,0.025])
grid on
ylim([0.95,1.025])
%%
fname = 'psat2019';
fig = jfigure(1,[600,550]);
clc
plot_tiles(m2.PSat{1}*1e12,p,'fig',fig,'clim',[90,220],'pair','mean','title','2019 PSat','clab','[pW]');
colormap(flipud(colormap('parula')))
fname = sprintf('tile_plots_%s.png',fname);
exportgraphics(fig,fullfile(figdir,fname),"Resolution",printres)

%%

idx = ismember(p.tile,[3,4,5,7,18]);
fig = figure(3);
hold off
v1 = nanmedian(poleffs_uniform,1);
scatter(m2.PSat{2}(idx,1)*1e12,m2.dPdT(idx,1)*1e12,8,p.tile(idx))
xlim([50,250])
%xlim([87,110])
% scatter(m2.dPdT*1e12,nanmedian(poleffs_uniform,1),8,p.tile)
% xlim([0.05,0.15])
%ylim([-0.015,0.025])
grid on
%ylim([0.95,1.025])

%%

fig = jfigure(124523,[800,500],[2,1]);
clf
ax1 = nexttile;
ax2 = nexttile;
hold(ax1,'on')
hold(ax2,'on')

for mcidx = 1:1000%length(fd_uniform.ch)
   if ~ismember(p.tile(fd_uniform.ch),11)
      continue 
   end
   if ismember(fd_uniform.ch(mcidx),ind0)
    ax = ax1;
   else
    ax = ax2;
   end
   
   A = fd_uniform.bparam{mcidx}(6:end);
   res = A-rps_get_mod_model(fd_uniform.parm_refit{mcidx},fd_uniform.rot{mcidx}+fd_uniform.phi_s(mcidx));
   plot(ax,fd_uniform.rot{mcidx}+fd_uniform.dk_cen(mcidx),res/fd_uniform.Amp(mcidx),'Color',[1 1 1]*0.35)
   
    
end
grid on


