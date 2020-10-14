function plot_cut_params()

load('data/b3rpsfiles.mat')
load('data/fitdata_20201013.mat')
addpath('Z:\\b3reduc\b3_analysis\util\')

prx = 2 * sind(p.r / 2) .* cosd(p.theta) * 180 / pi;
pry = 2 * sind(p.r / 2) .* sind(p.theta) * 180 / pi;

%%

% {Param, lower bound, upper bound, title}
param_array = {...
    {prx(fd.ch)-fd.bparam(:,1),-0.3,0.3,'x_{obs}-x_{fit}'},... 
    {pry(fd.ch)-fd.bparam(:,2),-0.3,0.3,'y_{obs}-y_{fit}'},...
    {fd.aparam(:,2),-0.05,0.05,'Xpol'},...
    {atand(tand(p.chi(fd.ch)+p.chi_thetaref(fd.ch)-fd.aparam(:,1))),-4,4,'\phi_{d,0}-\phi_{d,fit}'}...
    };




close all

for prmind = length(param_array)
    arr = param_array{prmind};
    param = arr{1};
    lb = arr{2};
    ub = arr{3};
    ptitle = arr{4};
    
    fig = figure();
    hist(param,length(param)/100)
    grid on
    title(['Raw: ' ptitle])
    
    if ~isnan(lb) & ~isnan(ub)
        make_cut_hist(param,lb,ub,ptitle)
    else
        disp(['Skipping: ' ptitle])
    end
end


function ind = make_cut_index(fd,param_array)
% Makes cut list. If true, channels PASS the cuts.

arr = param_array{1};
ind = true(size(arr{1}));
for prmind = 1:length(param_array)
    arr = param_array{prmind};
    ind = ind & inrange(arr{1},arr{2},arr{3});
    
end

% Toss chans without pairs:
% Do it by DK since they're separated in time and one channel may be
% picked up later in the schedule.

pair_index = zeros(size(fd.ch));
schlist = unique(fd.sch);
for schind = 1:length(schlist)
    
    ind = ind & ismember(fd.ch(fitind),chind);
    
    
end



function make_cut_hist(param,lb,ub,ptitle)
plim = [lb-abs(lb), ub+ub];
edges = plim(1):diff(plim)/100:plim(2);

N = histc(param,edges);
mxn = nanmax(N);

fig = figure();
set(fig,'Position',[900,300,420,420]);
bar(edges,N,'histc');
hold on
plot([1 1]*lb,[0 1]*2*mxn,'r')
plot([1 1]*ub,[0 1]*2*mxn,'r')
grid on
ylim([-0.1, 1.1*mxn])
xlim(plim)
title(ptitle)


function index = inrange(A,B,C)

index = (A >= B) & (A <= C);
