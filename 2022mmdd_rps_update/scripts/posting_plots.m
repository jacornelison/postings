function posting_plots()


%% Make histograms per detector with (or without) tilt (or homing) correction
close all

prefix = fullfile('C:','Users','James','Documents','GitHub','postings','2022mmdd_rps_update','');
paths = {fullfile(prefix,'data','2021_pole_notilt_params','params',''),fullfile(prefix,'data','2021_pole_params','params','')};
schedtype = {'9','_','0'};
pols = {'A','B'};
titles = {'No Tilt Correction', 'Tilt Corrected'};
names = {'no_tilt','with_tilt'};
homenames = {'','_cut'};
resolution = 15;
lims = {[-2.5, -1.5],...
    [86.5 87.5]};

fridgeind = [1, 10, 14, 20];
poorind = [32];

cutind = [fridgeind poorind];
typeind = 2;
for homeind = 2%1:2
    
    for pltind = 2%1:length(paths)
        
        dirname = paths{pltind};
        
        param_files = dir(fullfile(dirname,['param_0' schedtype{typeind} '*']));
        for chind = 1:2
            
            fname = ['type9_hist_pol_' pols{chind} '_' names{pltind} homenames{homeind}];
            det_angles = [];
            
            for fileind = 1:length(param_files)
                load(fullfile(dirname,param_files(fileind).name))
                det_angles(end+1) = atand(tand(i_param{chind}.phi_d));
            end
            
            if homeind == 2
                ind = cutind;
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
            ylim([0,18])
            text(mean(lims{chind})+0.3, 4, ...
                {sprintf('Mean: %0.2f',nanmean(det_angles)),...
                sprintf('STD: %0.2f',nanstd(det_angles))})
            title({['Type-9 Histogram Pol ' pols{chind}],titles{pltind}})
            saveas(fig,fullfile(prefix,'figs2',[fname '.png']))
        end
    end
end

%% Plot Pol A fits as a function of schedule number.

pols = {'A','B'};
fig = figure(3);
clf; hold on;
plts = [];
for chind = 1:2
lcolors = lines;
pltind = 2;
dirname = paths{pltind};
typeind = 1;
param_files = dir(fullfile(dirname,['param_0' schedtype{typeind} '*']));
%param_files = dir(fullfile(dirname,'param_0*'));
fname = ['type9_hist_pol_' pols{chind} '_' names{pltind} homenames{homeind}];
det_angles = [];

for fileind = 1:length(param_files)
    load(fullfile(dirname,param_files(fileind).name))
    
    det_angles(end+1) = atand(tand(i_param{chind}.phi_d+0*mean(i_param{chind}.tilt)));
    
end

plts(chind) = plot(1:length(det_angles),det_angles-nanmedian(det_angles),'Color',lcolors(chind,:));
plot(1:length(det_angles),det_angles-nanmedian(det_angles),'.','Color',lcolors(chind,:));
plts(3) = plot(fridgeind, det_angles(fridgeind)-nanmedian(det_angles),'k.','MarkerSize',10);
plts(4) = plot(poorind, det_angles(poorind)-nanmedian(det_angles),'r.','MarkerSize',10);
plts(5) = plot([37 38], det_angles([37 38])-nanmedian(det_angles),'.','Color',lcolors(3,:),'MarkerSize',10);
plts(6) = plot([39 40], det_angles([39 40])-nanmedian(det_angles),'.','Color',lcolors(4,:),'MarkerSize',10);
plts(7) = plot([41 42], det_angles([41 42])-nanmedian(det_angles),'.','Color',lcolors(5,:),'MarkerSize',10);
plts(8) = plot([43], det_angles([43])-nanmedian(det_angles),'.','Color',lcolors(6,:),'MarkerSize',10);
plts(9) = plot([44 45], det_angles([44 45])-nanmedian(det_angles),'.','Color',lcolors(7,:),'MarkerSize',10);
plts(10) = plot([46], det_angles([46])-nanmedian(det_angles),'.','Color',lcolors(8,:),'MarkerSize',10);
grid on
xlabel('Schedule Number')
ylabel('Pol Angle - median [Degrees]')
xlim([0,length(param_files)+1])
end
pols{end+1} = 'Blown Fridge';
pols{end+1} = 'Poor Fits';
pols{end+1} = 'DK = 23';
pols{end+1} = 'DK = 45';
pols{end+1} = 'DK = 68';
pols{end+1} = 'DK = 90';
pols{end+1} = 'DK = 135';
pols{end+1} = 'DK = 174';
legend(plts,pols,'Location','southwest')
title('Type-9 Per Det Angles per schedule #')
saveas(fig,fullfile(prefix,'figs2','angle_vs_schedule.png'))


%% Compare two modulation curves at 0 and 5.4545

rot = -180:30:180;
modfunc = @(x) (cosd(2*(rot-x))+1)/2;
chind = 1;

load(fullfile(dirname,param_files(7).name))
mod0 = i_param{chind}.bparam(6:18);
mod0 = mod0./max(mod0);
a0 = atand(tand(i_param{chind}.phi_d+90));
modsim0 = i_param{chind}.amodel./max(i_param{chind}.amodel);

load(fullfile(dirname,param_files(17).name))
a5 = atand(tand(i_param{chind}.phi_d+90));
mod5 = i_param{chind}.bparam(6:18);
mod5 = mod5./max(mod5);
modsim5 = i_param{chind}.amodel./max(i_param{chind}.amodel);

close all
fig = figure(1);
set(fig,'Position',[700,100,1000,500])
clf;
subplot(3,1,2)
hold on;
plot(rot,modsim0)
plot(rot,modsim5)
%plot(rot,modfunc(a0))
%plot(rot,modfunc(a5))
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
plot(rot,mod0-modsim0)
plot(rot,mod5-modsim5,'--')
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