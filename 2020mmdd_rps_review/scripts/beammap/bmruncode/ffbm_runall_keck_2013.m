function ffbm_runall_keck_2013(bmnum,demod,fullmap,winmap,winmaplarge,...
    html,plotmap)
% ffbm_runall_keck_2013(bmnum,demod,fullmap,winmap,winmaplarge,html,plotmap)
%
% Wrapper function to go from arcfiles to plotted maps for a single
% beammapping run.
% In your working directory, there should be a symlink "beammaps" pointing
% to /n/bicepfs2/keck/beammaps or /n/bicepfs3/bicep3/beammaps
%
% This function processes beammaps by the "run number" documented on our
% html logbook beammap tables.  It relies on a .csv file,
% beammaps/bmrunlist_201?.csv, which should be made and updated as beam maps
% come in and are processed.
%
% General usage notes: For demodulation and mapmaking, you can either run
% interactively (arg = 1) or farm it out (arg = 2).  Of course for things to
% work correctly, you need all the arcfiles demodded before mapmaking, and
% both full and windowed maps made before plotting (which is only done
% interactively for now).  Before plotting, you will need to change the
% 'filename' field in the .csv to correspond to the map/mapwin files.  So
% the least-complicated way (without farming) to run would be:
% >> ffbm_runall_keck_2013(1,1,1,1,0,0); % For bm run 1 in 2015
%    Then modify the .csv field 'filename' with the timestamp found in
%    beammaps/maps (can do this with make_expt_year_bmcsv.m)
% >> ffbm_runall_keck_2013(1,0,0,0,1,1); % And the plots should be made,
%    with the HTML already in there!
% >> for ii = 2:??; ffbm_runall_keck_2013(ii,2,0,0,0,0); end; % Farm demods

p = get_bm_info(2013);

% Select options.
bmopt.t1 = p.t1{bmnum};
bmopt.t2 = p.t2{bmnum};
bmopt.expt = 'keck';
bmopt.run = 'k2013';
bmopt.rxNum = 'all';
bmopt.notlastscan = 1; % just for diagnostic purposes
bmopt.demodtype = 'square';
bmopt.demoddir = 'beammaps/demod/';
%bmopt.demodshift = []; % Nominal shift in bm_demod

% Option to re-run demodulation on farmed jobs that failed
% ONLY USE IF YOU KNOW WHAT YOU'RE DOING
%bmopt.update = 1;
 
% Pointing model hack to avoid SINGULARITY
mjd = date2mjd(datenum('20130201','yyyymmdd'));
bmopt.pm = get_pointing_model(mjd);
bmopt.pm.az_tilt_ha = 0;
bmopt.pm.az_tilt_lat = 0;
bmopt.pm.el_tilt = 0;

% Turn off all warnings which swamp the output file
warning('off','all')

switch demod
  case 1
    bm_demod(bmopt);
  case 2
    s = gen_stamp();
    mkdir('farmfiles/demod');
    tmpfile = ['farmfiles/demod/k13demod_bm' sprintf('%02i',bmnum) '_' s];
    cmd = 'bm_demod(bmopt)';
    farmit(tmpfile,cmd,'var',{'bmopt'},'queue',...
	'general,serial_requeue,itc_cluster,regal','mem',80000,'maxtime',1440);
end


bmopt.el_res = 0.05;
bmopt.source_az = 0;
bmopt.source_el = 0;
bmopt.rot = [];
bmopt.component = 'cos';
bmopt.mapcoord = 'azel_ap';
bmopt.maptype = 'full';
bmopt.dofit = 0;
bmopt.save_map = 1;
bmopt.bmdir = 'beammaps/maps/';
bmopt.suffix = '';

% Option to re-run mapmaking on farmed jobs that failed
% ONLY USE IF YOU KNOW WHAT YOU'RE DOING
%bmopt.update = 1;

switch fullmap
  case 1
    bm_makemap(bmopt);
  case 2
    s = gen_stamp();
    mkdir('farmfiles/fullmap')
    tmpfile = ['farmfiles/fullmap/k13fullmap_bm' sprintf('%02i',bmnum) '_' s];
    cmd = 'bm_makemap(bmopt)';
    farmit(tmpfile,cmd,'var',{'bmopt'},'queue',...
	'general,serial_requeue,itc_cluster,regal','mem',100000,'maxtime',1440);
end

bmopt.source_az = 0.6;
bmopt.source_el = 3;
bmopt.maptype = 'window';
bmopt.dofit = 1;
bmopt.beamcen = 'obs';

switch winmap
  case 1
    bm_makemap(bmopt);
  case 2
    s = gen_stamp();
    mkdir('farmfiles/winmap')
    tmpfile = ['farmfiles/winmap/k13winmap_bm' sprintf('%02i',bmnum) '_' s];
    cmd = 'bm_makemap(bmopt)';
    farmit(tmpfile,cmd,'var',{'bmopt'},'queue',...
	'general,serial_requeue,itc_cluster,regal','mem',100000,'maxtime',1440);
end

if winmaplarge
  bmopt.winsize = 8;
  bmopt.suffix = '_8deg';
end

switch winmaplarge
  case 1
    bm_makemap(bmopt);
    bmopt.suffix = '';
  case 2
    s = gen_stamp();
    mkdir('farmfiles/winmap')
    tmpfile = ['farmfiles/winmap/k13winmapl_bm' sprintf('%02i',bmnum) '_' s];
    cmd = 'bm_makemap(bmopt)';
    farmit(tmpfile,cmd,'var',{'bmopt'},'queue',...
	'general,serial_requeue,itc_cluster,regal','mem',100000,'maxtime',1440);
    bmopt.suffix = '';
end

% Make sure the .csv file has the right filename, obtained from the
% timestamp on the maps generated
bmopt.filename = p.filename{bmnum};

% Make HTML before maps so we can look at them as they generate
bmopt.author = 'KSK';
bmopt.bmnum = p.number(bmnum);
bmopt.sched = p.sched{bmnum};
bmopt.pager = 1;

if html
  make_bmhtml(bmopt)
end

% Sometimes the large map coords have a very far away point, so plot just
% around the source
bmopt.fullaxis = [-20 20 -15 25];

switch plotmap
 case 1
  bm_plotter(bmopt)
 case 2
  s = gen_stamp();
  mkdir('farmfiles/plotter')
  tmpfile = ['farmfiles/plotter/k13plotter_bm' sprintf('%02i',bmnum) '_' s];
  cmd = 'bm_plotter(bmopt)';
  farmit(tmpfile,cmd,'var',{'bmopt'},'queue',...
         'general,serial_requeue,itc_cluster,regal','mem',40000,'maxtime',1440);
end



return
