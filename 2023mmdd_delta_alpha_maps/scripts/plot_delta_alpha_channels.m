function plot_delta_alpha_channels

addpath('z:/pipeline')
addpath('z:/pipeline/util')
addpath('z:/pipeline/beammap')
addpath('z:/dev/sims')

%%
result_ttl = {'Legacy','Reanalysis','BK'};
subset_ttl = {'Lower','Middle','Upper','All'};
sernames = {'6607','6608','6609','6606';...
    '6626','6627','6628','6621';...
    '','','','3553'};
daughters = {'f','g','h'};
yrs = {'2016','2017','2018'};
figdir = '~/postings/2023mmdd_delta_alpha_maps/figs/';
for resind = 1:length(result_ttl)
    for subind = 4%1:length(subset_ttl)
        ser = sernames{resind,subind};
        for yrind = 1:length(yrs)
            daught = daughters{yrind};
    fname = sprintf('maps/%s/0018_%s_filtp3_weight3_gs_dp1100_jack01.mat',ser,daught);
    if ~exist(fname,'file')
        continue
    end

    load(fname,'coaddopt')
    ttlname = sprintf('%s | %s | %s',result_ttl{resind},subset_ttl{subind},yrs{yrind});
    fig = figure(1);
    fig.Position(3:4) = [900 800];
    clf;
    plot_tiles(ismember(coaddopt.ind.e,coaddopt.ind.rgl),coaddopt.p,...
        'pair','mean','clim',[0 1],'title',ttlname,'fig',1);
    
    fname = sprintf('plot_tiles_%s_%s.png',ser,daught);
    saveas(fig,fullfile(figdir,fname),'png');

        end
    end
end








