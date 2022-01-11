function posting_plots()


%% Make histograms per detector with (or without) tilt (or homing) correction
close all

prefix = fullfile('C:','Users','James','Documents','GitHub','postings','2022mmdd_rps_update','');
paths = {fullfile(prefix,'data','2021_pole_notilt_params','params',''),fullfile(prefix,'data','2021_pole_params','params','')};
pols = {'A','B'};
titles = {'No Tilt Correction', 'Tilt Corrected'};
names = {'no_tilt','with_tilt'};
homenames = {'','_cut'};
resolution = 10;
lims = {[-2.5, -1.5],...
    [86.5 87.5]};

for homeind = 1:2
    
    for pltind = 1:length(paths)
        
        dirname = paths{pltind};
        
        param_files = dir(fullfile(dirname,'param_09*'));
        for chind = 1:2
            
            fname = ['type9_hist_pol_' pols{chind} '_' names{pltind} homenames{homeind}];
            det_angles = [];
            
            for fileind = 1:length(param_files)
                load(fullfile(dirname,param_files(fileind).name))
                det_angles(end+1) = atand(tand(i_param{chind}.phi_d));
            end
            
            if homeind == 2
                ind = [10 12 14 16:18];
                det_angles(ind) = NaN;
            end
            
            fig = figure(chind);
            clf; hold on;
            edges = lims{chind}(1):diff(lims{chind})/(resolution-1):lims{chind}(2);
            N = histc(det_angles,edges);
            bar(edges,N,'histc')
            %hist(det_angles,10)
            grid on
            xlabel('Pol Angle [Degrees]')
            ylabel('N')
            ylim([0,8])
            text(mean(lims{chind})+0.3, 4, ...
                {sprintf('Mean: %0.2f',nanmean(det_angles)),...
                sprintf('STD: %0.2f',nanstd(det_angles))})
            title({['Type-9 Histogram Pol ' pols{chind}],titles{pltind}})
            saveas(fig,fullfile(prefix,'figs',[fname '.png']))
        end
    end
end

%% Plot Pol A fits as a function of schedule number.

pols = {'A','B'};
fig = figure(1);
clf; hold on;
plts = [];
for chind = 1:2
lcolors = lines;
pltind = 1;
dirname = paths{pltind};
param_files = dir(fullfile(dirname,'param_09*'));
fname = ['type9_hist_pol_' pols{chind} '_' names{pltind} homenames{homeind}];
det_angles = [];

for fileind = 1:length(param_files)
    load(fullfile(dirname,param_files(fileind).name))
    
    det_angles(end+1) = atand(tand(i_param{chind}.phi_d));
    
end

plts(chind) = plot(1:length(det_angles),det_angles-nanmedian(det_angles),'Color',lcolors(chind,:));
plot(1:length(det_angles),det_angles-nanmedian(det_angles),'.','Color',lcolors(chind,:));
plts(3) = plot([10 14], det_angles([10,14])-nanmedian(det_angles),'k.','MarkerSize',10);
plts(4) = plot([12 16:18], det_angles([12 16:18])-nanmedian(det_angles),'r.','MarkerSize',10);
grid on
xlabel('Schedule Number')
ylabel('Pol Angle - median [Degrees]')
xlim([0,19])
end
pols{end+1} = 'Blown Fridge';
pols{end+1} = 'Poor Fits';
legend(plts,pols,'Location','southwest')
title('Type-9 Per Det Angles per schedule #')
saveas(fig,fullfile(prefix,'figs','angle_vs_schedule.png'))


%% Compare two modulation curves at 0 and 5.4545

rot = -180:30:180;
modfunc = @(x) (cosd(2*(rot-x))+1)/2;

load(fullfile(dirname,param_files(7).name))
mod0 = i_param{1}.bparam(6:18);
mod0 = mod0./max(mod0);
a0 = atand(tand(i_param{1}.phi_d+90));
load(fullfile(dirname,param_files(17).name))
a5 = atand(tand(i_param{1}.phi_d+90));
mod5 = i_param{1}.bparam(6:18);
mod5 = mod5./max(mod5);

close all
fig = figure(1);
set(fig,'Position',[700,100,1000,500])
clf;
subplot(3,1,2)
hold on;
plot(rot,modfunc(a0))
plot(rot,modfunc(a5))
ylim([0,1.2])
legend({'Sch04 Sim', 'Sch17 Sim'})
grid on
ylabel('Simmed Mod Curves')

subplot(3,1,1)
hold on;
plot(rot,mod0)
plot(rot,mod5)
ylim([0,1.2])
legend({'Sch04 Pol A', 'Sch17 Pol A'})
grid on
ylabel('Real Mod Curves')

subplot(3,1,3)
hold on;
plot(rot,mod0-mod5)
plot(rot,modfunc(a0)-modfunc(a5),'--')
ylim([-0.1, 0.2])
grid on
legend({'Sch04 - Sch14','Sim @ 0 - Sim @ -5.454545'})
xlabel('Grid Angle [Degrees]')
ylabel('Residuals')
saveas(fig,fullfile(prefix,'figs','mod_curves.png'))

%%
figure(2)
clf; hold on;
gof = [];
det_angle = [];
for fileind = 1:18

    load(fullfile(dirname,param_files(fileind).name))
    gof(fileind) = i_param{2}.agof;
    det_angle(fileind) = i_param{2}.phi_d;
end

plot(gof,det_angle)