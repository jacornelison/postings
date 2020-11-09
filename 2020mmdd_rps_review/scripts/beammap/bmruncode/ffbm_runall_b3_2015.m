function ffbm_runall_b3_2015(bmnum,demod,fullmap,winmap,plotmap,html)
% ffbm_runall_b3_2015(bmnum,demod,fullmap,winmap,plotmap,html)
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
% >> ffbm_runall_b3_2015(1,1,1,1,0,0); % For bm run 1 in 2015
%    Then modify the .csv field 'filename' with the timestamp found in
%    beammaps/maps (can do this with make_expt_year_bmcsv.m)
% >> ffbm_runall_b3_2015(1,0,0,0,1,1); % And the plots should be made, 
%    with the HTML already in there!

p = get_bm_info(2015);

% Select options.
bmopt.t1 = p.t1{bmnum};
bmopt.t2 = p.t2{bmnum};
bmopt.expt = 'b3';
bmopt.run = 'b3r6';
bmopt.rxNum = 0;
bmopt.notlastscan = 1;
bmopt.demodtype = 'square';
bmopt.demoddir = 'beammaps/demod/';
% bmopt.demodshift = 3; % Defaults based on dates in bm_demod

% Pointing model hack to avoid SINGULARITY
mjd = date2mjd(datenum('20150201','yyyymmdd'));
bmopt.pm = get_pointing_model(mjd);
bmopt.pm.az_tilt_ha = 0;
bmopt.pm.az_tilt_lat = 0;
bmopt.pm.el_tilt = 0;

switch demod
  case 1
    bm_demod(bmopt);
  case 2
    s = gen_stamp();
    tmpfile = ['scratch/demod/demod_bm' sprintf('%02i',bmnum) '_' s];
    cmd = 'bm_demod(bmopt)';
    farmit(tmpfile,cmd,'var',{'bmopt'},'queue',...
	'general,serial_requeue,itc_cluster','mem',40000,'maxtime',1440);
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

switch fullmap
  case 1
    bm_makemap(bmopt);
  case 2
    s = gen_stamp();
    tmpfile = ['scratch/fullmap/fullmap_bm' sprintf('%02i',bmnum) '_' s];
    cmd = 'bm_makemap(bmopt)';
    farmit(tmpfile,cmd,'var',{'bmopt'},'queue',...
	'general,serial_requeue,itc_cluster','mem',40000,'maxtime',1440);
end

bmopt.source_az = -177.5;
bmopt.source_el = 2.75;
bmopt.maptype = 'window';
bmopt.dofit = 1;

switch winmap
  case 1
    bm_makemap(bmopt);
  case 2
    s = gen_stamp();
    tmpfile = ['scratch/winmap/winmap_bm' sprintf('%02i',bmnum) '_' s];
    cmd = 'bm_makemap(bmopt)';
    farmit(tmpfile,cmd,'var',{'bmopt'},'queue',...
	'general,serial_requeue,itc_cluster','mem',40000,'maxtime',1440);
end

% Make sure the .csv file has the right filename, obtained from the
% timestamp on the maps generated
bmopt.filename = p.filename{bmnum};

if plotmap
  bm_plotter(bmopt);
end

bmopt.author = 'KSK';
bmopt.bmnum = p.number(bmnum);
bmopt.sched = p.sched{bmnum};

if html
  make_bmhtml(bmopt);
end

return
