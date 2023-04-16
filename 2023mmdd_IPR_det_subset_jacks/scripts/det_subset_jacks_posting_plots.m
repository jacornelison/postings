function get_sim_global_rotation_B18()

%% Add relevant paths
addpath('z:/pipeline')
addpath('z:/pipeline/util')
addpath('z:/pipeline/beammap')
addpath('z:/dev/sims')
set(groot,'defaulttextinterpreter','latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(groot,'defaultAxesFontSize',12)

%% Initialize relevant variables
clear all
close all
% per sigtype
signame = [2,3,5,6,7,8];
sigtitle = {'Unlensed LCDM','Gaussian Dust','Lensed LCDM','Sign-Flip Noise','Lensed-LCDM+Noise','Lensed-LCDM+Noise+Dust'};
% per purification
%purename = {'','pureB_','matrix_'};
%puretitle = {'No purification','Kendrick Pure-B','Matrix Purification'};
purename = {'','_matrix'};
puretitle = {'No purification','Matrix Purification'};

% Per sernum
%sername = {'6600','6614','6614'};
sername = {'6622','6600','6614'};
daughter = {'fgh','fgh','fgh'};
% Per fixed/free axes
axname = {'_free','_fixed'};

%figdir = fullfile('C:','Users','James','Documents','GitHub','postings','20221213_birefringence_fits_to_B2018','figs','');
%figdir = fullfile('C:','Users','James','Documents','GitHub','postings','20221213_birefringence_pol_rot_sims','figs','');
%figdir = fullfile('C:','Users','James','Documents','GitHub','postings','20221213_birefringence_data_subset_fits','figs','');
%figdir = fullfile('C:','Users','James','Documents','GitHub','postings','20230214_IPR_B18_mat_pure','figs','');
%figdir = fullfile('C:','Users','James','Documents','GitHub','postings','20230321_IPR_B18_subsets','figs','');
%figdir = fullfile('C:','Users','James','Documents','GitHub','postings','20230321_IPR_high_ell','figs','');
%figdir = fullfile('C:','Users','James','Documents','GitHub','postings','20230321_IPR_linest','figs','');
figdir = fullfile('C:','Users','James','Documents','GitHub','postings','2023mmdd_IPR_det_subset_jacks','figs','');

final = load('z:/dev/sims/final/3553/real_supfac');

% per fittype
specnames = {'EB','TB','EB+TB'};
crossname = {'','_cross'};

cm = colormap('lines');
close all

%% Grab the pol angles for a given sernum

polopt = struct;
polopt.offdiag = 0;

fitopt = struct;
fitopt.signame = [2 3 5 6 7 8];
fitopt.sername = '6607';
fitopt.daughter = 'h';
%fitopt.purename = '';
%fitopt.iscross = true;
fitopt.covtype = 'normal';
fitopt.polopt = polopt;
fitopt.estimator = 'linear';
fitopt.bpcm_simset = '6614/xxx8_fgh_filtp3_weight3_gs_dp1100_jack0_matrix_cm_overfreq.mat';
fitopt.usebins = 11:15;

clear accumulate_pol_rot_fits
tic;
accumulate_pol_rot_fits(fitopt)
toc



%% Plot the APS now
%

clc
fig = figure(1);
fig.Position(3:4) = [600 600];
clf;

figname = 'blank_6x6.png';
saveas(fig,fullfile(figdir,figname))
%ylims = [0.25 5 0.4]; % For Posting 1
%ylims = [0.4 5 0.01]; % for posting 2
ylims = [0.4 5 0.4]; % for posting 3
ylims = [0.5 10 2];

usebins = 2:10;

for crossind = 1:2
    for serind = 2:3%1:length(sername)
        load(sprintf('z:/dev/sims/%s_%s_global_pol_fits%s_bins_%i_%i.mat',...
            sername{serind},daughter{serind},crossname{crossind},min(usebins),max(usebins)))

        for sigind = 1:length(signame)
            for pureind = 1:2%1:length(purename)
                % Organize the APS

                nsims = length(ps{1,1}.alpha)/length(specnames);

                for fitspecind = 0%:length(specnames)
                    if fitspecind == 0
                        apsloop = 1:length(specnames):length(ps{1,1}.alpha);
                    else
                        apsloop = fitspecind:length(specnames):length(ps{1,1}.alpha);
                    end

                    count = 1;
                    ell = ps{1,1}.l;
                    [EB, TB, BB] = deal(NaN(17,1,nsims));
                    [EBexp, TBexp, BBexp] = deal(NaN(17,1,nsims));
                    for apsind = apsloop
                        EB(:,:,count) = ps{sigind,pureind}.rrot(apsind).real(:,6);
                        TB(:,:,count) = ps{sigind,pureind}.rrot(apsind).real(:,5);
                        BB(:,:,count) = ps{sigind,pureind}.rrot(apsind).real(:,4);

                        EBexp(:,:,count) = ps{sigind,pureind}.rrot(apsind).expv(:,6);
                        TBexp(:,:,count) = ps{sigind,pureind}.rrot(apsind).expv(:,5);
                        BBexp(:,:,count) = ps{sigind,pureind}.rrot(apsind).expv(:,4);

                        count = count+1;
                    end

                    ell = ell(usebins);

                    EB = squeeze(EB(usebins,:,:));
                    TB = squeeze(TB(usebins,:,:));
                    BB = squeeze(BB(usebins,:,:));

                    EBexp = squeeze(EBexp(usebins,:,:));
                    TBexp = squeeze(TBexp(usebins,:,:));
                    BBexp = squeeze(BBexp(usebins,:,:));


                    for axind = 1%1:2
                        clf;
                        subplot(3,1,1)
                        hold on;

                        z1 = plot(repmat(ell,1,size(EB,2)),EB,'Color',[1 1 1]*0.7);
                        z2 = plot(ell,nanmean(EB,2),'k');
                        legs = {'Sim Rlz','Mean Sims'};
                        if fitspecind ~= 0
                            legs(3:4) = {'Expv','Mean Expv'};
                            z3 = plot(repmat(ell,1,size(EB,2)),EBexp,'Color',[0 0.7 1]*0.7);
                            z4 = plot(ell,nanmean(EBexp,2),'b');
                            legend([z1(1), z2(1), z3(1), z4(1)] ,legs,'Location','southwest')
                        else
                            legend([z1(1), z2(1)] ,legs,'Location','southwest')
                        end
                        grid on
                        xlabel('$\ell$','interpreter','latex')
                        ylabel('$D^{EB}_{\ell}[\mu K^2]$','interpreter','latex')
                        if axind == 2
                            ylim([-1 1]*ylims(1))
                        end

                        subplot(3,1,2)
                        hold on;
                        plot(repmat(ell,1,size(TB,2)),TB,'Color',[1 1 1]*0.7)
                        plot(ell,nanmean(TB,2),'k')
                        if fitspecind ~=0
                            plot(repmat(ell,1,size(TB,2)),TBexp,'Color',[0 0.7 1]*0.7)
                            plot(ell,nanmean(TBexp,2),'b')
                        end
                        grid on
                        xlabel('$\ell$','interpreter','latex')
                        ylabel('$D^{TB}_{\ell}[\mu K^2]$','interpreter','latex')
                        if axind == 2
                            ylim([-1 1]*ylims(2))
                        end

                        subplot(3,1,3)
                        hold on;
                        plot(repmat(ell,1,size(BB,2)),BB,'Color',[1 1 1]*0.7)
                        plot(ell,nanmean(BB,2),'k')
                        %plot(repmat(ell,1,size(TB,2)),BBexp,'Color',[0 0.7 1]*0.7)
                        %plot(ell,nanmean(BBexp,2),'b')
                        grid on
                        xlabel('$\ell$','interpreter','latex')
                        ylabel('$D^{BB}_{\ell}[\mu K^2]$','interpreter','latex')
                        if axind == 2
                            ylim([0 1]*ylims(3))
                        end

                        sgtitle(sprintf('%s, %s',sigtitle{sigind},puretitle{pureind}))
                        if fitspecind ==0
                            figname = sprintf('aps_%s_%s_%s%01i%s%s_bins_%i_%i.png',...
                                sername{serind},daughter{serind},purename{pureind},...
                                signame(sigind),axname{axind},crossname{crossind},...
                                min(usebins),max(usebins));
                        else
                            figname = sprintf('aps_%s_%s_%s%01i%s_%s_bins_%i_%i.png',...
                                sername{serind},daughter{serind},purename{pureind},...
                                signame(sigind),axname{axind},specnames{fitspecind},...
                                min(usebins),max(usebins));
                        end
                        saveas(fig,fullfile(figdir,figname))
                    end
                end
            end
        end
    end
end

%% Plot histograms of Alpha

if 1
fig = figure(1);
fig.Position(3:4) = [1100 400];
clf;

figname = 'blank_6x4.png';
saveas(fig,fullfile(figdir,figname))
end
scaling = {[1 1 1],[0.865 0.803 0.848]};
scalename = {'','_corr'};

crossname = {'','_cross'};

%offs = [0,0,0,0,-0.8,-0.8,-0.8,-0.8,0.5,-0.5,0.25,1,0];
set(groot,'defaultAxesFontSize',11)

% Loops we wanna do.
crossloop = 1:2;
scaleloop = 1;%1:2;
serloop = 3;%1:length(sername);
pureloop = 1:2;%1:length(purename);

bins = {2:10};%{2:10, 2:15};
covnames = {'normal','normal_repsim_6614xxx8','legacy'};%{'normal','legacy'};
covnames = covnames(2);
repnames = {''};%{'','_alloff'};
pseudonames = {''};%,'_pseudo'};
offdiagnames = 0;%{2,102};
% Put all of our stupid bullshit loops and options into one loop.
combos = product(crossloop, scaleloop, serloop, pureloop, bins,...
    covnames,repnames,pseudonames,offdiagnames);

clc
for combidx = 1:length(combos)

    [crossind, scaleind, serind, pureind,...
        usebins,covtypename,repname,pseudoname,offdiagname] = deal(combos{combidx}{:});
    

    fname = sprintf('%s_%s_global_pol_fits_bins_%i_%i_offdiag_%i%s%s%s%s%s',sername{serind},daughter{serind},...
            min(usebins),max(usebins),offdiagname,purename{pureind},crossname{crossind},...
            ['_' covtypename],repname);
    disp(fname)
    load(fullfile('z:/dev/sims/',fname))

    for sigind = 1:size(ps,1)
        figname = sprintf('%s_%s_sig_%i_alpha_hist_bins_%i_%i_offdiag_%i%s%s%s%s%s',sername{serind},daughter{serind},...
            ps{sigind}.signame,min(usebins),max(usebins),offdiagname,purename{pureind},crossname{crossind},...
            ['_' covtypename],repname);
        
        do_alpha_hist_plots(ps{sigind},figname,figdir)
    end
end



%% Make a blank image
fig = figure(1);
fig.Position(3:4) = [600 400];
clf;

figname = 'blank.png';
saveas(fig,fullfile(figdir,figname))

%% Make a blank image
fig = figure(1);
fig.Position(3:4) = [600 600];
clf;



%% Look at the distributions of data subsets

clc
close all
sername = {'6600','6603','6604','6605'};
sigind = 1;
pureind = 1;
for serind = 1:length(sername)
    apsname = sprintf('z:/dev/sims/aps/%s/xxx%i_h_filtp3_weight3_gs_dp1100_jack01_%scm_overfreq.mat',sername{serind},signame(sigind),purename{pureind});
    load(apsname,'coaddopt')
    plot_tiles(double(ismember(1:2640,coaddopt.ind.rgl100)),coaddopt.p,'title','FPU coverage');
    %length(coaddopt.ind.rgl100a)
    fname = sprintf('%s_tile_plot.png',sername{serind});
    saveas(serind,fullfile(figdir,fname))
end

%% Histogram of upper-middle-lower angle data splits

load('z:/dev/sims/phi_dummy_bicep3_20170101_struct.mat')

dp = k.diffphi;
q = quantile(dp,[0,0.3333,0.6666,1]);
idx = logical(k.hasmeas);

fig = figure(101);
fig.Position(3:4) = [440 430];
clf; hold on;
t = tiledlayout(1,1);
t.TileSpacing = 'compact';
t.Padding = 'compact';
nexttile()
hold on

edges = -1.5:3.5/71:0.6;
N = histc(dp(idx),edges);
bar(edges,N,'histc')

yedges = [0, max(N)*1.2];
clr = {'b','r','g'};
for pltind=1:3
    plot([1 1]*q(pltind),yedges,'--','Color',clr{pltind})
    plot([1 1]*q(pltind+1),yedges,'--','Color',clr{pltind})
end
text(-1.35,max(N)*1.1,'Lower')
text(-0.88,max(N)*1.1,'Middle')
text(-0.45,max(N)*1.1,'Upper')
grid on
xlim([min(edges), max(edges)])
ylim(yedges)
xlabel('$\phi_{pair,RPS}\;-\;\phi_{pair,obs}$ [Deg]')
ylabel('N')
title('Histogram of as-measured angle vs. fiducial')
pbaspect([1 1 1])
fname = 'hist_diffphi.png';
saveas(fig,fullfile(figdir,fname))


%% Plot the actual mean vs. expected (with error bars?)
clc

Aexp = [-1.069 -0.796 -0.778 -0.468];


cm = colormap('lines');
mk = {'-','--',':'};
signums = [2,5, 7, 8];
% Loop over signals
clear z
for sigind = 1:size(Amean,2)
    fig = figure(1);
    fig.Position(3:4) = [600 500]*0.75;
    clf; hold on;

    % Loop over fit-types
    for fitind = 1:size(Amean,3)
        Dexp = Aexp;
        Dsim = squeeze(Amean(:,sigind,fitind))';
        Dstd = squeeze(Astd(:,sigind,fitind))';
        z(fitind) = errorbar(Dexp,Dsim,Dstd/2,'Color',cm(fitind,:),'LineStyle',mk{fitind});
        plot(Dexp,Dsim,'.','Color',cm(fitind,:),'MarkerSize',12)
        grid on

    end

    lims = [-1.2 -0.3];
    plot(lims,lims,'k--')

    xlim(lims)
    ylim(lims)

    xlabel('Expected Angle [Degrees]')
    ylabel('Actual Angle [Degrees]')
    legend(z,{'EB','TB','EB+TB'},'location','northwest')

    text(-1.125,-1.05,'Lower','Rotation',90)
    text(-0.85,-0.79,'Middle','Rotation',90)
    text(-0.7,-1.05,'B2018+RPS','Rotation',90)
    text(-0.4,-.7,'Upper','Rotation',90)


    fname = sprintf('subset_mean_compare_type_%i',signums(sigind));
    %saveas(fig,fullfile(figdir,fname),'png')

end



%% Now plot the differences between each subset
% and compare to the expected
clc

Aexp = [-1.069 -0.796 -0.778 -0.468];

ind1 = [1 2 1];
ind2 = [2 4 4];
cm = colormap('lines');
mk = {'-','--',':'};
signums = [2,5, 7, 8];
% Loop over signals
clear z
for sigind = 1:size(Amean,2)
    fig = figure(1);
    fig.Position(3:4) = [600 500]*0.75;
    clf; hold on;

    % Loop over fit-types
    for fitind = 1:size(Amean,3)
        Dexp = Aexp(ind2)-Aexp(ind1);
        Dsim = squeeze(Amean(ind2,sigind,fitind)-Amean(ind1,sigind,fitind))';
        Dstd = sqrt(squeeze(Astd(ind2,sigind,fitind)).^2+squeeze(Astd(ind1,sigind,fitind)).^2)';
        z(fitind) = errorbar(Dexp,Dsim,Dstd/2,'Color',cm(fitind,:),'LineStyle',mk{fitind});
        plot(Dexp,Dsim,'.','Color',cm(fitind,:),'MarkerSize',12)
        grid on

    end

    lims = [0.15 0.65];
    plot(lims,lims,'k--')

    xlim(lims)
    ylim(lims)

    xlabel('Expected Difference [Degrees]')
    ylabel('Actual Difference [Degrees]')
    legend(z,{'EB','TB','EB+TB'},'location','northwest')

    text(0.25,0.3,'Mid - Low','Rotation',90)
    text(0.32,0.4,'Upper - Mid','Rotation',90)
    text(0.6,0.35,'Upper - Low','Rotation',90)


    fname = sprintf('subset_diff_compare_type_%i',signums(sigind));
    %saveas(fig,fullfile(figdir,fname),'png')

end


%end of main function
end

function printcol(x,cstart,rend)
txt = '<td class=''';
if cstart
    txt = [txt 'cstart'];
end
if rend
    txt = [txt ' rend'];
end
if ischar(x)
    txt = [txt '''>%s</td>\n'];
else
    txt = [txt '''>%0.4f</td>\n'];
end
fprintf(txt,x);

end

function do_covmat_plots(ps,fname, figdir)
if ~exist('covfig','var')
    covfig = randi([0,1000],1,1);
end
fig = figure(covfig);
fig.Position(3:4) = [1100 400];
clf;
t = tiledlayout(1,3);
t.TileSpacing = 'compact';
t.Padding = 'compact';


end


function do_alpha_hist_plots(ps,fname,figdir,doebtb)
% ps should be ps{idx}
specnames = {'EB-Only','TB-Only','EB+TB'};

if ~exist("doebtb",'var')
    doebtb=true;
end

if ~doebtb
    specmax = 2;
else
    specmax = 3;
end

histfig = figure(1);
histfig.Position(3:4) = [367*specmax 400];
clf;
t = tiledlayout(1,specmax);
t.TileSpacing = 'compact';
t.Padding = 'compact';

a = ps.alpha;
x = size(a);
for specind = 1:specmax
    edges = (-1:0.075:1)*0.5;
    N = histc(a(specind,:),edges);
    
    nexttile(specind)
    bar(edges,N,'histc')
    xlim([min(edges) max(edges)])
    grid on
    Z = a(specind,:);
    title({sprintf('%s-fit',specnames{specind}), ...
        sprintf('Mean: %0.4f STD:%0.4f EOM:%0.4f',mean(Z),std(Z),std(Z)./sqrt(length(Z)))})
    xlabel('$\alpha$ [Degrees]')
    pbaspect([1,1,1])
end

disp(fname)
saveas(histfig,fullfile(figdir,fname),'png');

end