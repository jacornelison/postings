function check_QU_map_rot

%% 
restoredefaultpath
addpath('z:/pipeline')
addpath('z:/pipeline/util')
addpath('z:/pipeline/beammap')
addpath('z:/dev/sims')
addpath('z:/dev/')
set(groot,'defaulttextinterpreter','latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(groot,'defaultAxesFontSize',12)
figdir = fullfile('C:','Users','jcornelison','Documents','GitHub','postings','20240305_IPR_rotation_scaling','figs','');


%%
% Loading the maps takes forever to run, so don't rerun if you've already fit params.
Ain = [0.5 -0.5 0.25 1.0];
Aout = NaN(10,4,3); % Q, U, Q+U

%%
clc
if 0
    for simind = 1:10

        map2 = load(sprintf('z:/dev/sims/maps/6600_%03i2_h_filtp3_weight3_gs_dp1100_jack01.mat',simind));
        [MAP2, ~]=make_map(map2.ac,map2.m,map2.coaddopt);
        for angind = 1:4
            if 0 % A check to make sure we can rotate things properly
                MAP1 = rotqumaps(MAP2,Ain(angind));
            else
                map1 = load(['z:/dev/sims/maps/661' num2str(angind-1) sprintf('_%03i2_h_filtp3_weight3_gs_dp1100_jack01.mat',simind)]);
                [MAP1, ~]=make_map(map1.ac,map1.m,map1.coaddopt);
            end

            mxfev = 100000;
            mxiter = 100000;
            options = optimset('TolFun',1e-10,'MaxIter',mxiter,'MaxFunEvals',mxfev,'Display','off');

            % Only use common pixels
            idx = isnan(MAP1.Q) | isnan(MAP1.U) | isnan(MAP2.Q) | isnan(MAP2.U);
            MAP1.Q(idx) = 0; MAP1.U(idx) = 0; MAP2.Q(idx) = 0; MAP2.U(idx) = 0;

            for fitmode = 1:3
                %fitfunc = @(x,y) y.Q*cosd(2*x)-y.U*sind(2*x);
                fitfunc = @(x,y) rotmodel(x,y,fitmode);
                switch fitmode
                    case 1
                        parm = lsqcurvefit(fitfunc,0,MAP1,reshape(MAP2.Q,[],1),[-10],[10],options);
                    case 2
                        parm = lsqcurvefit(fitfunc,0,MAP1,reshape(MAP2.U,[],1),[-10],[10],options);
                    case 3
                        parm = lsqcurvefit(fitfunc,0,MAP1,[reshape(MAP2.Q,[],1);reshape(MAP2.U,[],1)],[-10],[10],options);
                end
                Aout(simind,angind,fitmode) = parm;
            end


            fprintf('ang %i,sim %i\n',angind,simind)
        end
    end
    %
    save(fullfile(figdir,'..','scripts','map_angle_fits'),'Aout','Ain')
else
    load(fullfile(figdir,'..','scripts','map_angle_fits'))
end

%%
clr = colormap('lines');
fig = figure(1); clf;
fig.Position(3:4) = [900 600]*1.1;
t = tiledlayout(2,3);
t.TileSpacing = 'compact';
t.Padding = 'compact';

plt_names = {'Fit Q Maps only','Fit U Maps Only','Fit Q+U Maps'};
for fitmode = 1:3
    nexttile(); hold on
    Ain_rep = repmat(Ain,10,1);
    Aout_rep = squeeze(Aout(:,:,fitmode));
    [C,S] = polyfit(reshape(Ain_rep,1,[]),reshape(Aout_rep,1,[]),1);
    fiterr = (inv(S.R)*inv(S.R)')*S.normr^2/S.df;


    %plot(Ain,Aout,'.','MarkerSize',14,color=clr(1,:))
    M = nanmedian(Aout_rep,1);
    S = nanstd(Aout_rep,[],1);
    errorbar(Ain,M,S,'.','MarkerSize',14,'Color',clr(1,:),'CapSize',10)
    %plot(Ain,Aout,color=clr(1,:))

    plot([-1 1]*1.2,[-1 1]*1.2,'k--')
    plot([-1,1]*1.2,polyval(C,[-1 1]*1.2),'k')
    xlim([-0.7 1]*1.2)
    ylim([-0.7 1]*1.2)
    grid on
    pbaspect([1 1 1])
    if fitmode == 3
    xlabel('Sim Input Angle [deg]')
    end
    if fitmode == 1
    ylabel('Map Output angle [deg]')
    end
    title(plt_names{fitmode})
    text(-0.75,0.4, ...
        { ...
        sprintf('Slope: %0.3f $\\pm$ %0.1E',C(1),fiterr(1,1)), ...
        sprintf('x0: %0.1f $\\pm$ %0.1E',C(2),fiterr(2,2))...
        },'Interpreter','latex')

end


% Cross Check: Input Angle Vs. Output Angle
% Independent of Number of ell-bins used to fit, so only use standard 9.
ang_in = [0.5 -0.5 0.25 1.0];
sername1 = {'6610','6611','6612','6613'};
ebscaling = 0.87;
lims = [-0.7 1]*1.2;
cmlines = colormap('lines');

[ang_out, s_out] = deal(NaN(1,4));
[ang_out_all] = deal(NaN(10,4));
for serind = 1:4
    load(sprintf('z:/dev/sims/angle_fits/%s_h_jack01_global_pol_fits_bins_2_10_offdiag_0_matrix_normal_repsim_6614xxx8.mat',sername1{serind}))
    ang_out(serind) = nanmean(ps{1}.alpha(1,:),2);
    s_out(serind) = nanstd(ps{1}.alpha(1,:),[],2);
    ang_out_all(:,serind) = ps{1}.alpha(1,:)';
end

[C,S] = polyfit(reshape(repmat(ang_in,10,1),1,[]),reshape(ang_out_all,1,[]),1);
fiterr = (inv(S.R)*inv(S.R)')*S.normr^2/S.df;

%fig = jfigure(401369,[350 350]);
nexttile(4); hold on;

plot([-10 10],[-10 10],'k--');
plot([-1,1]*1.2,polyval(C,[-1 1]*1.2),'k');
%z(2) = plot(ang_in,ang_out,'.-','Color',cmlines(1,:),'MarkerSize',14);
errorbar(ang_in,ang_out,s_out,'.','MarkerSize',14,'Color',clr(1,:),'CapSize',10);
%z(3) = plot(ang_in,ang_out./ebscaling,'.-','Color',cmlines(4,:),'MarkerSize',14);
%legend(z(2:3),{'Uncorrected','Corrected'},'Location','northwest')
    text(-0.75,0.4, ...
        { ...
        sprintf('Slope: %0.3f $\\pm$ %0.1E',C(1),fiterr(1,1)), ...
        sprintf('x0: %0.1f $\\pm$ %0.1E',C(2),fiterr(2,2))...
        },'Interpreter','latex')

xlim(lims)
ylim(lims)
grid on
pbaspect([1 1 1])
xlabel('Sim Input Angle [Degrees]')
ylabel('Spectra Output Angle [Degrees]')
title('Fit EB Spectra')


% Now do TB
ang_in = [0.5 -0.5 0.25 1.0];
sername1 = {'6610','6611','6612','6613'};
ebscaling = 0.87;
lims = [-0.7 1]*1.2;
cmlines = colormap('lines');

[ang_out, s_out] = deal(NaN(1,4));
[ang_out_all] = deal(NaN(10,4));
for serind = 1:4
    load(sprintf('z:/dev/sims/angle_fits/%s_h_jack01_global_pol_fits_bins_2_10_offdiag_0_matrix_normal_repsim_6614xxx8.mat',sername1{serind}))
    ang_out(serind) = nanmean(ps{1}.alpha(2,:),2);
    s_out(serind) = nanstd(ps{1}.alpha(2,:),[],2);
    ang_out_all(:,serind) = ps{1}.alpha(2,:)';
end

[C,S] = polyfit(reshape(repmat(ang_in,10,1),1,[]),reshape(ang_out_all,1,[]),1);
fiterr = (inv(S.R)*inv(S.R)')*S.normr^2/S.df;
nexttile(5); hold on;

plot([-10 10],[-10 10],'k--');
plot([-1,1]*1.2,polyval(C,[-1 1]*1.2),'k');

errorbar(ang_in,ang_out,s_out,'.','MarkerSize',14,'Color',clr(1,:),'CapSize',10);

    text(-0.75,0.4, ...
        { ...
        sprintf('Slope: %0.3f $\\pm$ %0.1E',C(1),fiterr(1,1)), ...
        sprintf('x0: %0.1f $\\pm$ %0.1E',C(2),fiterr(2,2))...
        },'Interpreter','latex')

xlim(lims)
ylim(lims)
grid on
pbaspect([1 1 1])
xlabel('Sim Input Angle [Degrees]')
%ylabel('Output Angle fit from sims [Degrees]')
title('Fit TB Spectra')

sgtitle({'Input Vs. Output Angle','Sims: LCDM-only, N=10'},'Interpreter','latex')
exportgraphics(fig,fullfile(figdir,'Aout_vs_Ain.png'),"Resolution",300)

function model = rotmodel(ang,mapin,fitmode)
mapout = rotqumaps(mapin,ang);

switch fitmode
    case 1
        model=reshape(mapout.Q,[],1);
    case 2
        model=reshape(mapout.U,[],1);
    case 3
        model=[reshape(mapout.Q,[],1); reshape(mapout.U,[],1)];
end

function [Q_out,U_out]=do_qu_rotation(Q_in,U_in,rotang)

[th,r]=cart2pol(Q_in,U_in);
[Q_out,U_out]=pol2cart(th+2*rotang,r);


