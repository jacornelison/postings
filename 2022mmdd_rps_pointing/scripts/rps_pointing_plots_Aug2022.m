function rps_pointing_plots_Aug2022()

figdir = 'c:/Users/James/Documents/GitHub/postings/2022mmdd_rps_pointing/figs/';
inrange = @(A,B,C) B <= A & A <= C;
outrange = @(A,B,C) A <= B | C <= A;

clc
% Moon-derived mirror params -- Use only these and nothing else!
mirror = struct();
mirror.height = 1.4592;
mirror.tilt= 44.8870;
mirror.roll = -0.0750;
rpsopt.mirror = mirror;

rpsopt.source.distance = 195.5;
% Fit for the source params given our mirror info:
source = rps_fit_source(fd,rpsopt,p,'')
rpsopt.source = source;


%% With new mirror and source parameters, update the pointing.
[fd.x,fd.y,phi] = beam_map_pointing_model(fd.az_cen,fd.el_cen,fd.dk_cen,model,'bicep3',mirror,source,[]);
fd.x = reshape(fd.x,size(fd.ch));
fd.y = reshape(fd.y,size(fd.ch));
fd.resx = reshape(prx(fd.ch),size(fd.ch))-fd.x;
fd.resy = reshape(pry(fd.ch),size(fd.ch))-fd.y;
[resth, resr] = cart2pol(fd.resx,fd.resy);
[fd.resx_rot, fd.resy_rot] = pol2cart(resth-fd.dk_cen*pi/180,resr);


% Organize the x/y pointing by obs
[xps, yps] = deal(NaN(len,2640));
for schind = 1:len
    for chind = 1:2640
        ci = find(fd.ch==chind & ismember(fd.schnum,scheds{schind}));

        if ~isempty(ci)

            %err(isnan(val))=1e10;
            xps(schind,chind) = mean(fd.x(ci));
            yps(schind,chind) = mean(fd.y(ci));
            
        end
    end
end

%% Fit for a rotation and scaling

fpu = fit_fpu_angle_and_scaling_from_xy(fdsch,p)

%%

winscale = 1.5;
scaling = 10;
fig = figure(1);
fig.Position(3:4) = [500*winscale 450*winscale];
clf;
cm = colormap('turbo');
clridx = floor(linspace(1,size(cm,1),3*19));

% Things dealing with projection
xlims = {[-1 1]*15 [-1 1]*15 [-1 1]*0.5};
ylims = {[-1 1]*15 [-1 1]*15 [-0.5 0.8]};
projlabels = {' [Degrees]','_m [Degrees]','_m [Meters]'};
projnames = {'','_mirror','_mirror'};

% Things dealing with fits
fittype = {'overall','perdk'};%,'persch'};


casename = {'none','ang','scale','both'};
parms = [fpu.angle,fpu.scaling,fpu.xtrans,fpu.ytrans];
for caseind = 1:4
    switch caseind
        case 1
            params = [0,1,0,0];
            x = xps;
            y = yps;
            caselabel = sprintf('FPU Rot: %0.3f^o, Scaling: %0.3f',params(1),params(2));
        case 2
            params = parms.*[1,0,0,0];
            params(2) = 1;
            x = params(2).*(xps.*cosd(params(1))-yps.*sind(params(1)));
            y = params(2).*(xps.*sind(params(1))+yps.*cosd(params(1)));
            caselabel = sprintf('FPU Rot: %0.3f^o, Scaling: %0.3f',params(1),params(2));
        case 3
            params = parms.*[0,1,0,0];
            x = params(2).*(xps.*cosd(params(1))-yps.*sind(params(1)));
            y = params(2).*(xps.*sind(params(1))+yps.*cosd(params(1)));
            caselabel = sprintf('FPU Rot: %0.3f^o, Scaling: %0.3f',params(1),params(2));
        case 4
            params = parms;
            x = params(2).*(xps.*cosd(params(1))-yps.*sind(params(1)));
            y = params(2).*(xps.*sind(params(1))+yps.*cosd(params(1)));
            caselabel = sprintf('FPU Rot: %0.3f^o, Scaling: %0.3f',params(1),params(2));
    end

projind = 1;
corrind = 2;
winscale = 1;
scaling = 10;
fig = figure(1);
fig.Position(3:4) = [500*winscale 450*winscale];
clf;
ind = 1:length(scheds);
quiver(prx,pry,(prx-nanmean(x(ind,:),1)')*scaling,(pry-nanmean(y(ind,:),1)')*scaling,0)
grid on;
xlim(xlims{projind})
ylim(ylims{projind})
xlabel(sprintf('X%s',projlabels{projind}))
ylabel(sprintf('Y%s',projlabels{projind}))
    title({sprintf('Beam Center Residuals, x%i', scaling),...
        sprintf('All-DK average, Tilt: %1.2f  Roll: %1.3f',mirror.tilt,mirror.roll),...
        caselabel...
        })
    figname = fullfile(figdir,sprintf('quiver_mean_fit_%s.png',casename{caseind}));
    saveas(fig,figname)

end