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
fitopt.sername = '6614';
fitopt.daughter = 'fgh';
fitopt.purename = '';
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


%% Plot power spectra

clc
matname = {'','_matrix'};
nsims = 1:30;
usebins = 2:16;
for matind = 1:2
    psname = sprintf('z:/dev/sims/6622_gh_global_pol_fits_bins_2_10_offdiag_0%s_cross_normal_repsim_6614xxx8.mat',matname{matind});
    %psname = sprintf('z:/dev/sims/6614_gh_global_pol_fits_bins_2_10_offdiag_0%s_cross_normal_repsim_6614xxx8.mat',matname{matind});
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

        for specind = 1:6
            D = squeeze(mean(aps(apsind).Cs_l(usebins,specind,nsims),3));
            D_sub = squeeze(mean(aps_sub(apsind).Cs_l(usebins,specind,nsims),3));
            Ddiff = (D_sub-D)./D;
            %Ddiff(abs(Ddiff)>1.5) = NaN;

            nexttile(rc2idx(specind,1))
            hold on;
            plot(usebins,D,'Color',cmlines(1,:))
            plot(usebins,D_sub,'Color',cmlines(2,:))
            plot(usebins,D,'.','Color',cmlines(1,:))
            plot(usebins,D_sub,'.','Color',cmlines(2,:))
            grid on
            if specind ==1
                title({'Spectra','Blue: Original, Red: 17+18 Subset'})
            end
            %ylabel({ylabs{specind},'$\ell(\ell+1)C_{\ell}\,[\muK^2]$'})
            ylabel({ylabs{specind},'$D_{\ell}\,[\mu K^2]$'})
            xlim([0 17])


            nexttile(rc2idx(specind,2))
            hold on;
            plot(usebins,Ddiff,'.','Color',cmlines(1,:))
            plot(usebins,Ddiff,'Color',cmlines(1,:))
            grid on
            if specind ==1
                title({'fractional difference','$(D_{sub}-D)/D$'},'Interpreter','Latex')
            end
            
            % mess with the y axes
            ysc = 1.2;
            ymin = max([min(Ddiff)*ysc,-1]);
            ymax = min([max(Ddiff)*ysc,1]);
            if ymin>0 & ymax>0
                ymin = 0;
            elseif ymin<0 & ymax<0
                ymax = 0;
            end
            ylim([ymin ymax])

            xlim([0 17])

        end
        figname = sprintf('specplot_sig%i%s',sigind,matname{matind});
        saveas(fig,fullfile(figdir,figname),'png')
        disp(figname)
    end
end

%% Overplot BB-powers



sigtitle = {'Unlensed LCDM','Gaussian Dust','Lensed LCDM','Sign-Flip Noise','Lensed-LCDM+Noise','Lensed-LCDM+Noise+Dust'};
psname = sprintf('z:/dev/sims/6622_gh_global_pol_fits_bins_2_10_offdiag_0%s_cross_normal_repsim_6614xxx8.mat',matname{2});
load(psname)
ps_sub = ps;
aps_sub = polstruct2aps(ps_sub);

usebins = 2:16;
fig = figure(2);
fig.Position(3:4) = [670 400];
clf; hold on;
apsloop = [2,3,4,6];
clear z;
for apsind = apsloop
    D = real(log10(squeeze(mean(aps_sub(apsind).Cs_l(usebins,4,:),3))));
    plot(usebins,D,'.','Color',cmlines(apsind,:),'MarkerSize',14)
    z(apsind) = plot(usebins,D,'Color',cmlines(apsind,:),'LineWidth',1);
end
grid on
%legend(z(apsloop),sigtitle(apsloop),'Location','northeastoutside')
legend(z(apsloop),sigtitle(apsloop),'Location','northwest')
ylim([-4 0])
xlim([0 17])
xlabel('$\ell$-bins')
ylabel('log10($D^{BB}_\ell$) [log10($\mu K^2$)]')
title('BB Power spectra, 17+18 subset')
figname = 'spectra_compare_17+18subset.png';
exportgraphics(fig,fullfile(figdir,figname),'Resolution',600)


%% Grab means/std's from simsets


sigind = 1;
clc
A = NaN(6,6);
usesims = 1:30;

for sigind = 1:6

fprintf('\n\n%i:\n',sigind)
load('z:/dev/sims/6614_fgh_global_pol_fits_bins_2_10_offdiag_0_matrix_cross_normal_repsim_6614xxx8.mat');
M1 = mean(ps{sigind}.alpha(1,usesims))/0.87;
S1 = std(ps{sigind}.alpha(1,usesims))/0.87;
fprintf('6614 | 2-10 | M: %1.4f | S: %1.4f\n',M,S1)


load('z:/dev/sims/6622_gh_global_pol_fits_bins_2_10_offdiag_0_matrix_cross_normal_repsim_6614xxx8.mat');
M2 = mean(ps{sigind}.alpha(1,usesims))/0.87;
S2 = std(ps{sigind}.alpha(1,usesims))/0.87;
fprintf('6622 | 2-10 | M: %1.4f | S: %1.4f\n',M,S2)
fprintf('%%-diff | 2-10 | %1.4f\n',((S2-S1)/S1))

load('z:/dev/sims/6614_fgh_global_pol_fits_bins_2_15_offdiag_0_matrix_cross_normal_repsim_6614xxx8.mat');
M3 = mean(ps{sigind}.alpha(1,usesims))/0.87;
S3 = std(ps{sigind}.alpha(1,usesims))/0.87;
fprintf('6614 | 2-15 | M: %1.4f | S: %1.4f\n',M,S3)

load('z:/dev/sims/6622_gh_global_pol_fits_bins_2_15_offdiag_0_matrix_cross_normal_repsim_6614xxx8.mat');
M4= mean(ps{sigind}.alpha(1,usesims))/0.87;
S4 = std(ps{sigind}.alpha(1,usesims))/0.87;
fprintf('6622 | 2-15 | M: %1.4f | S: %1.4f\n',M,S4)

fprintf('%%-diff | 2-15 | %1.4f\n',((S4-S3)/S3))

A(sigind,:) = [S1 S2 ((S2-S1)/S1), S3, S4, ((S4-S3)/S3)];
if sigind==6
    B = [M1 S1/sqrt(30) M3 S3/sqrt(30); M2 S2/sqrt(30) M4 S4/sqrt(30);];% NaN ((S2-S1)/S1) NaN ((S4-S3)/S3)];
end

end

%%
hd = {};
rl = {'L-LCDM','Noise','L-LCDM+N+Dust'};
simple_html_table(A([3 4 6],:),hd,rl)    

%%

rl = {'B18','17+18 Subset'};%,'Frac-Diff'};
simple_html_table(B,{},rl)





