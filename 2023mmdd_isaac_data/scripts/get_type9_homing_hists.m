function get_type9_homing_hists()

addpath('z:/dev/rps/')
addpath('z:/pipeline/util')
p22 = load('z:/dev/rps/fpu_data_obs.mat'); % pointing from 2022
p18 = load('z:/dev/rps/fpu_data_B18.mat'); % RGL's from B18
p = p22.p;
p_ind = p18.p_ind;

% shortcuts for 0/90 orientiations
inda = p_ind.a;
indb = p_ind.b;

ind0 = [inda(ismember(inda,find(p.mce~=0))) ...
    indb(ismember(indb,find(p.mce==0)))];
ind90 = [indb(ismember(indb,find(p.mce~=0))) ...
    inda(ismember(inda,find(p.mce==0)))];

load('z:/dev/rps/sch_type9.mat')
fd_type9 = load('z:/dev/rps/rps_beam_fits_type9_21feb_rerun.mat');
fd_type9 = fd_type9.fd;
fd_type9 = rps_cut_fitdata(fd_type9,p,p_ind,false);
gitdir = fullfile('C:','Users','James','Documents','');
figdir = fullfile(gitdir,'GitHub','postings','2023mmdd_isaac_data','figs');

%% Refit the 10-inc measurements and 
scheds = [1:34 50:68];

fig = figure(4923); clf;
fig = figure(4924); clf;
phis = NaN(1,length(fd_type9.ch));
for chind = 1:length(fd_type9.ch)
    if ~ismember(fd_type9.schnum(chind),scheds)
        continue
    end

    ch = fd_type9.ch(chind);
    mxfev = 100000;
    mxiter = 100000;
    options = optimset('TolFun',1e-10,'MaxIter',mxiter,'MaxFunEvals',mxfev,'Display','off');
    lb = [-20 -0.5 -10 -10 0];
    ub = [339.999 0.5 10 10 1e6];
    if ch == 696
        angguess = 0;
    else
        angguess = 90;
    end
    R = fd_type9.bparam{chind}(6:end);
    a = fd_type9.rot{chind};
    % Downsample 10-inc measurements
    if ismember(fd_type9.schnum(chind),49:58)
        fig = figure(4923);
        hold on;
        R = R(1:3:end);
        a = a(1:3:end);
        plot(a,R)
    else
        fig = figure(4924);
        hold on;
        plot(a,R)
    end

    guess = [angguess 0 0 0 1];
    parm = lsqcurvefit(@rps_get_mod_model,guess,a+fd_type9.phi_s(chind),R,lb,ub,options);
    phis(chind) = atand(tand(parm(1)));
    %fprintf('%i: %2.2f\n',fd_type9.schnum(chind),parm(1))

end

%% Plot the histograms

inds = {...
...ismember(fd_type9.schnum,[49:58]) & ismember(fd_type9.ch,ind0),...
~ismember(fd_type9.schnum,[59:68]) & ismember(fd_type9.ch,ind0),...
ismember(fd_type9.schnum,[59:68]) & ismember(fd_type9.ch,ind0),...
...ismember(fd_type9.schnum,[49:58]) & ismember(fd_type9.ch,ind90),...
~ismember(fd_type9.schnum,[59:68]) & ismember(fd_type9.ch,ind90),...
ismember(fd_type9.schnum,[59:68]) & ismember(fd_type9.ch,ind90),...
    };

edgescell = {linspace(-2.1,-1.6,10),...
linspace(-2.1,-1.6,10),...
linspace(86.8,87.4,10),...
linspace(86.8,87.4,10),...
};

fig = figure(7254);
fig.Position(3:4) = [560 650];
clf;
t = tiledlayout(2,2);
t.TileSpacing = 'compact';
t.Padding = 'compact';

for valind = 1:length(inds)
    nexttile()
    idx = inds{valind};
    V = phis(idx);
    edges = edgescell{valind};
    N = histc(V,edges);
    bar(edges,N,'histc')
    M = nanmean(V);
    S = nanstd(V);
    L = length(V(~isnan(V)));
    ttlstats = sprintf('M: %2.3f | S: %2.3f | N: %i | EOM: %2.3f',M,S,L,S/sqrt(L));
    grid on
    xlim([edges(1) edges(end)])
    pbaspect([1 1 1])
    if valind == 1
        title({'Homed Every Rasterset',ttlstats})
        ylabel('Pol 0')
    elseif valind == 2
        title({'No Homing',ttlstats})
    elseif valind == 3
        title({ttlstats})
        ylabel('Pol 90')
        xlabel('Pol Angle [Deg]')
    elseif valind == 4
        title({ttlstats})
        xlabel('Pol Angle [Deg]')
    end

end

sgtitle({'Type 9 RPS Schedules on BICEP3','with and without homing'})

fname = 'type9_homing_hists';
saveas(fig,fullfile(figdir,fname),'png')