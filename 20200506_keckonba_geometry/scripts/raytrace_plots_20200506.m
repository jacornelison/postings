figdir = 'C:\Users\james\Documents\GitHub\postings\20200506_keckonba_geometry\figs\';

% Look at OG Keck

% Forebaffle dimensions from MSG's 20120727 posting
% Optical params from Lorenzo's instrument parameter spreadsheet
% Mount values from Colin's pointing model defaults
expt.fb_h = 0.7366+0.14; %plus dist to aperture
expt.fov = 9.6*2;
expt.win_d = 0.264; % aperture diameter
expt.dk_off = [-1.1, 1.5964];
expt.el_off = [0., 1.1750];
expt.az_off = [0., 0.];
expt.n_rx = 5;
expt.gs_dim = [7.62 5.1]; % DASI groundshield
expt.expt = 'Keck';

for i = 45:60

    expt.min_el = i;
    
    [parm, fig]  = s4_gs_study(expt,'PLOT',true,'LEGEND',true,'INTEXT',true,'fixwindist',0.54,'anim',true,'showalt',false);
    
    print(fig,[figdir 'keck_el_' num2str(i)],'-dpng')
    
end


%% Now Keck on the BA Mount

% FB dimensions are the same as above.
% Taken from SOLIDWorks "AsBuild" models in the repo.
% Optical params from Lorenzo's instrument parameter spreadsheet
expt.fb_h = 0.7366+0.14; %plus dist to aperture
expt.fov = 9.6*2;
expt.win_d = 0.264;
expt.dk_off = [0., 0.9271];
expt.el_off = [0., 2.3];
expt.az_off = [0., 0.];
expt.n_rx = 4;
expt.gs_dim = [7.62 5.1]; % DASI groundshield
expt.expt = 'Keck on BA';

for i = 45:60

    expt.min_el = i;
    
    [parm, fig]  = s4_gs_study(expt,'PLOT',true,'LEGEND',true,'INTEXT',true,'fixwindist',0.94,'anim',true,'showalt',false);
    
    print(fig,[figdir 'keckonba_el_' num2str(i)],'-dpng')
    
end



%% Now do just BA

% Mount dimensions are the same as above.
% FB dimensions from Fig. 3.7 of NWP's 20190207_BA_Forebaffle_Concepts posting in BA logbook.
expt.fb_h = 1.01;
expt.fov = 29.6;
expt.win_d = 0.69;
expt.dk_off = [0., 0.9271];
expt.el_off = [0., 2.3];
expt.az_off = [0., 0.];
expt.n_rx = 4;
expt.gs_dim = [7.62 5.1]; % DASI groundshield
expt.expt = 'BA';

for i = 45:60

    expt.min_el = i;
    
    [parm, fig]  = s4_gs_study(expt,'PLOT',true,'LEGEND',true,'INTEXT',true,'fixwindist',0.94,'anim',true,'showalt',false);
    
    print(fig,[figdir 'BA_el_' num2str(i)],'-dpng')
    
end

