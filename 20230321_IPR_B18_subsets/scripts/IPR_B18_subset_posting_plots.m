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

%% Grab the pol angles for a given sernum

polopt = struct;
polopt.offdiag = 0;

fitopt = struct;
fitopt.signame = [2 3 5 6 7 8];
fitopt.sername = '6622';
fitopt.daughter = 'gh';
%fitopt.purename = '';
fitopt.iscross = true;
fitopt.covtype = 'normal';
fitopt.polopt = polopt;
%fitopt.estimator = 'linear';
fitopt.bpcm_simset = '6614/xxx8_fgh_filtp3_weight3_gs_dp1100_jack0_matrix_cm_overfreq.mat';
fitopt.usebins = 2:10;

clear accumulate_pol_rot_fits
tic;
accumulate_pol_rot_fits(fitopt)
toc


%%
clc
matname = {'','_matrix'};
for matind = 1:2
    psname = sprintf('z:/dev/sims/6622_fgh_global_pol_fits_bins_2_10_offdiag_0%s_cross_normal_repsim_6614xxx8.mat',matname{matind});
    load(psname)
    ps_sub = ps;
    aps_sub = polstruct2aps(ps_sub);

    psname = sprintf('z:/dev/sims/6614_fgh_global_pol_fits_bins_2_10_offdiag_0%s_cross_normal_repsim_6614xxx8.mat',matname{matind});
    load(psname)
    aps = polstruct2aps(ps);

    for apsind = 1:6
        sigind = ps{apsind}.signame;
        ylabs = {...
            sprintf('$T_{8}\\times T_{%i}$',sigind),...
            sprintf('$T_{8}\\times E_{%i}$',sigind),...
            sprintf('$E_{8}\\times E_{%i}$',sigind),...
            sprintf('$B_{8}\\times B_{%i}$',sigind),...
            sprintf('$T_{8}\\times B_{%i}$',sigind),...
            sprintf('$E_{8}\\times B_{%i}$',sigind),...
            };

        fig = figure(1);
        fig.Position(3:4) = [560 700];
        clf; hold on;
        rc2idx = reshape(1:12,2,[])';
        t = tiledlayout(6,2);
        %t.Padding = 'tight';
        t.TileSpacing = 'tight';

        for specind = [1:6]

            D = squeeze(mean(aps(apsind).Cs_l(:,specind,1:10),3));
            D_sub = squeeze(mean(aps_sub(apsind).Cs_l(:,specind,1:10),3));
            Ddiff = (D_sub-D)./D;
            %Ddiff(abs(Ddiff)>1.5) = NaN;

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
            %ylabel({ylabs{specind},'$\ell(\ell+1)C_{\ell}\,[\muK^2]$'})
            ylabel({ylabs{specind},'$D_{\ell}\,[\mu K^2]$'})
            xlim([0 18])


            nexttile(rc2idx(specind,2))
            hold on;
            plot(Ddiff,'.','Color',cmlines(1,:))
            plot(Ddiff,'Color',cmlines(1,:))
            grid on
            if specind ==1
                title({'fractional difference','$(D_sub-D})/D$'})
            elseif ismember(specind,[5 6])
                ylim([-1 1])
            end
            xlim([0 18])

        end
        figname = sprintf('specplot_sig%i%s',sigind,matname{matind});
        saveas(fig,fullfile(figdir,figname),'png')
        disp(figname)
    end
end

%% Overplot BB-powers

sigtitle = {'Unlensed LCDM','Gaussian Dust','Lensed LCDM','Sign-Flip Noise','Lensed-LCDM+Noise','Lensed-LCDM+Noise+Dust'};
psname = sprintf('z:/dev/sims/6622_fgh_global_pol_fits_bins_2_10_offdiag_0%s_cross_normal_repsim_6614xxx8.mat',matname{2});
    load(psname)
    ps_sub = ps;
    aps_sub = polstruct2aps(ps_sub);

fig = figure(2);
fig.Position(3:4) = [900 400];
clf; hold on;

clear z;
for apsind = 1:6
    D = real(log10(squeeze(mean(aps_sub(apsind).Cs_l(:,4,:),3))));
    plot(D,'.','Color',cmlines(apsind,:))
    z(apsind) = plot(D,'Color',cmlines(apsind,:));
end
grid on
legend(z,sigtitle,'Location','northeastoutside')
ylim([-6 0])

%% Grab means/std's from simsets


sigind = 1;
clc
for sigind = 1:6

fprintf('\n\n%i:\n',sigind)
load('z:/dev/sims/6600_fgh_global_pol_fits_bins_2_10_offdiag_0_matrix_cross_normal_repsim_6614xxx8.mat');
M = mean(ps{sigind}.alpha(1,:));
S1 = std(ps{sigind}.alpha(1,:));
fprintf('6600 | 2-10 | M: %1.4f | S: %1.4f\n',M,S1)

load('z:/dev/sims/6622_fgh_global_pol_fits_bins_2_10_offdiag_0_matrix_cross_normal_repsim_6614xxx8.mat');
M = mean(ps{sigind}.alpha(1,:));
S3 = std(ps{sigind}.alpha(1,:));
fprintf('6622 | 2-10 | M: %1.4f | S: %1.4f\n',M,S3)

fprintf('%%-diff | 2-10 | %1.4f\n',((S3-S1)/S1))

load('z:/dev/sims/6600_fgh_global_pol_fits_bins_2_15_offdiag_0_matrix_cross_normal_repsim_6614xxx8.mat');
M = mean(ps{sigind}.alpha(1,:));
S2 = std(ps{sigind}.alpha(1,:));
fprintf('6600 | 2-15 | M: %1.4f | S: %1.4f\n',M,S2)

load('z:/dev/sims/6622_fgh_global_pol_fits_bins_2_15_offdiag_0_matrix_cross_normal_repsim_6614xxx8.mat');
M = mean(ps{sigind}.alpha(1,:));
S4 = std(ps{sigind}.alpha(1,:));
fprintf('6622 | 2-15 | M: %1.4f | S: %1.4f\n',M,S4)

fprintf('%%-diff | 2-15 | %1.4f\n',((S4-S2)/S2))

end

    






