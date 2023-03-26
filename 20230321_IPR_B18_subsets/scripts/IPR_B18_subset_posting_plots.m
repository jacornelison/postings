function IPR_B18_subset_posting_plots()

%% Add relevant paths / default settings
addpath('z:/pipeline')
addpath('z:/pipeline/util')
addpath('z:/pipeline/beammap')
addpath('z:/dev/sims')
set(groot,'defaulttextinterpreter','latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(groot,'defaultAxesFontSize',12)

figdir = fullfile('C:','Users','James','Documents','GitHub','postings','20230321_IPR_B18_subsets','figs','');

cmlines = colormap('lines');
%%

load('z:/dev/sims/6622_fgh_global_pol_fits_bins_2_10_offdiag_0_matrix_cross_normal_repsim_6614xxx8.mat')
ps_sub = ps;
aps_sub = polstruct2aps(ps_sub);

load('z:/dev/sims/6600_fgh_global_pol_fits_bins_2_10_offdiag_0_matrix_cross_normal_repsim_6614xxx8.mat');
aps = polstruct2aps(ps);

fig = figure(1);
fig.Position(3:4) = [560 700];
clf; hold on;

rc2idx = reshape(1:12,2,[])';
t = tiledlayout(6,2);
t.Padding = 'tight';
t.TileSpacing = 'tight';
clc
ylabs = {'$T_{8}\times T_{2}$',...
    '$T_{8}\times E_{2}$',...
    '$E_{8}\times E_{2}$',...
    '$B_{8}\times B_{2}$',...
    '$T_{8}\times B_{2}$',...
    '$E_{8}\times B_{2}$',...
};



for specind = [1:6]
    
    
    D = squeeze(mean(aps(2).Cs_l(:,specind,1:10),3));
    D_sub = squeeze(mean(aps_sub(2).Cs_l(:,specind,1:10),3));
    
    nexttile(rc2idx(specind,1))
    hold on;
    plot(D,'Color',cmlines(1,:))
    plot(D_sub,'Color',cmlines(2,:))
    plot(D,'.','Color',cmlines(1,:))
    plot(D_sub,'.','Color',cmlines(2,:))
    grid on
    if specind ==1
        title({'Spectra','Blue: B18, Red: B18-subset'})
    end
    ylabel(ylabs{specind})

    nexttile(rc2idx(specind,2))
    hold on;
    plot((D-D_sub)./D,'.','Color',cmlines(1,:))
    plot((D-D_sub)./D,'Color',cmlines(1,:))
    grid on
    if specind ==1
        title({'fractional difference','$(D-D_{sub})/D$'})
    end
    %ylim([-1 1]*0.1)
end



