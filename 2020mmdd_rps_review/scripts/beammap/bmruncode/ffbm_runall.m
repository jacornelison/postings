function ffbm_runall(expt,year,bmnums,realorsim,timestream,fullmap,buddymap,...
                     winmap,winbuddymap,winmapsizes,html,plotmap)
% ffbm_runall(expt,year,bmnums,timestream,demod,fullmap,buddymap,...
%             winmap,winbuddymap,winmapsizes,html,plotmap)
%
% Wrapper function to go from arcfiles to plotted maps:
% 1. Loads arcfiles, demodulates, and makes timestream masks (bm_timestream)
% 3. Bins demodded arcfiles into full maps (bm_makemap)
% 4. Takes full map and windows/fits around main beam (bm_winfitmap)
% 5. Repeats steps 2/3 for buddy beams (i.e. 180 deg around boresight)
% 6. Makes html page for standard posting
% 7. Makes plots for standard posting
%
% This function processes beam maps by the "run number" documented on our
% html logbook beam map tables.  It relies on a .csv file,
% beammaps/bmrunlist_201?.csv, which should be made and updated as beam
% maps come in and are processed.  
%
% In your working directory, there should be a symlink "beammaps"
% pointing to /n/bicepfs?/expt/beammaps
%
% expt: 'bicep2','keck','bicep3'
% year: 2012, etc
% bmnums: vector of beam maps to process (e.g. 1, 1:10, etc.)
% realorsim: 'real' use real data, 'sim' use simulated timestream data
% timestream:  1 interactive, 2 to farm per schedule
% fullmap: 1 interactive, 2 to farm per schedule 
% winmap:  1 interactive, 2 to farm per schedule
% winmapsizes: vector of size of windowed maps to make 
%              (rad from center, in degrees) -- for full set, would like
%              [1.2,2,4,6,8]
% buddymap: 1 interactive, 2 to farm (flips theta by 180 degrees)
% html, plotmap: 1 interactive

p = get_bm_info(year);

queues = 'general,serial_requeue,itc_cluster,kovac';
memtimestream = 50000;
memfullmap = 180000; % This can vary quite a bit
memwinmap = 120000; % have to load large fullmap...
memplot = 40000;
maxtime = 720;

% Any options that are common to all run numbers in a given year
% are set in this function
bmopt = get_bmopts(expt,year,realorsim,buddymap,winbuddymap);

% Timestream step for sim data requires more time/mem
if strcmp(realorsim,'sim') && timestream==2
  memtimestream = 80000;
  maxtime = 1200;
end

for ii = 1:length(bmnums)

  number = bmnums(ii);
  bmopt.number = number;
  bmopt.t1 = p.t1{number};
  bmopt.t2 = p.t2{number};
  bmopt.filename = p.filename{number};
  bmopt.bmnum = p.number(number);
  bmopt.sched = p.sched{number};
  bmopt.fullmaptime = p.filename{number};
  if strcmp(expt,'keck')
    bmopt.mirror = p.mirror{number};
  end
  [bmopt.source_az bmopt.source_el] = get_sourcepos(expt,year,number);
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Process raw arcfiles into demodulated timestream data
  switch timestream
    case 1
      disp(['Running bm_timestream on ' num2str(year) ' run ' ...
            num2str(number)])
      bm_timestream(bmopt);
    case 2
      disp(['Farming bm_timestream on ' num2str(year) ' run ' ...
            num2str(number)])
      s = gen_stamp();
      mkdir('farmfiles/timestream');
      tmpfile = ['farmfiles/timestream/' expt '_' num2str(year) ...
                 '_timestream_bm' sprintf('%02i',number) '_' s];
      cmd = 'bm_timestream(bmopt)';
      farmit(tmpfile,cmd,'var',{'bmopt'},'queue', ...
             queues,'mem',memtimestream,'maxtime',maxtime);
  end
  

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Make full map
  % Source az,el only need to be specified in this step if making
  % map in x,y coords.
  
  switch fullmap
    case 1
      disp(['Running bm_makemap on ' num2str(year) ' run ' num2str(number)])
      bm_makemap(bmopt);
    case 2
      disp(['Farming bm_makemap on ' num2str(year) ' run ' num2str(number)])
      s = gen_stamp();
      mkdir('farmfiles/fullmap');
      tmpfile = ['farmfiles/fullmap/' expt '_' num2str(year) '_fullmap_bm' ...
                 sprintf('%02i',number) '_' s];
      cmd = 'bm_makemap(bmopt)';
      farmit(tmpfile,cmd,'var',{'bmopt'},'queue', ...
             queues,'mem',memfullmap,'maxtime',maxtime);
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Make full BUDDY map

  switch buddymap
    case 1
      disp(['Running bm_makemap (for buddymaps) on ' num2str(year) ' run ' num2str(number)])
      bmopt.buddy = 1;
      bm_makemap(bmopt);
    case 2
      disp(['Farming bm_makemap (for buddymaps) on ' num2str(year) ' run ' num2str(number)])
      bmopt.buddy = 1;
      s = gen_stamp();
      mkdir('farmfiles/fullmap_buddy');
      tmpfile = ['farmfiles/fullmap/' expt '_' num2str(year) '_fullmap_bm' ...
                 sprintf('%02i',number) '_' s];
      cmd = 'bm_makemap(bmopt)';
      farmit(tmpfile,cmd,'var',{'bmopt'},'queue', ...
             queues,'mem',memfullmap,'maxtime',maxtime);
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Window full map and fit
  % For this to work, the full map must have been made and the
  % corresponding filename recorded in bmrunlist.csv
  
  for jj = 1:length(winmapsizes)
    
    bmopt.winmapsize = winmapsizes(jj);
    switch winmap
      case 1
        disp(['Running bm_winfitmap on ' num2str(year) ' run ' ...
              num2str(number) ', ' num2str(winmapsizes(jj)) ' deg'])
        bm_winfitmap(bmopt);
      case 2
        disp(['Farming bm_winfitmap on ' num2str(year) ' run ' ...
              num2str(number) ', ' num2str(winmapsizes(jj)) ' deg'])
        s = gen_stamp();
        mkdir('farmfiles/winmap');
        tmpfile = ['farmfiles/winmap/' expt '_' num2str(year) '_winmap_bm' ...
                   sprintf('%02i',number) '_' s];
        cmd = 'bm_winfitmap(bmopt)';
        farmit(tmpfile,cmd,'var',{'bmopt'},'queue', ...
               queues,'mem',memwinmap,'maxtime',maxtime);
    end
    
  end

  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Window full BUDDY map and fit
  % For this to work, the full map must have been made and the
  % corresponding filename recorded in bmrunlist.csv
  
  for jj = 1:length(winmapsizes)
    
    bmopt.winmapsize = winmapsizes(jj);
    switch winbuddymap
      case 1
        disp(['Running bm_winfitmap (for buddymaps) on ' num2str(year) ' run ' ...
              num2str(number) ', ' num2str(winmapsizes(jj)) ' deg'])
        bmopt.buddy = 1;
        bmopt.fullmaptime = p.filename{number};
        bm_winfitmap(bmopt);
      case 2
        disp(['Farming bm_winfitmap (for buddymaps) on ' num2str(year) ...
              'run ' num2str(number) ', ' num2str(winmapsizes(jj)) ' deg'])
        bmopt.fullmaptime = p.filename{number};
        bmopt.buddy = 1;
        s = gen_stamp();
        mkdir('farmfiles/winmap_buddy');
        tmpfile = ['farmfiles/winmap/' expt '_' num2str(year) '_winmap_bm' ...
                   sprintf('%02i',number) '_' s];
        cmd = 'bm_winfitmap(bmopt)';
        farmit(tmpfile,cmd,'var',{'bmopt'},'queue', ...
               queues,'mem',memwinmap,'maxtime',maxtime);
    end
    
  end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Make standard HTML page before the plots so you can look at them
  % while they generate
  
  if html
    make_bmhtml(bmopt);
  end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Plot maps
  
  % We chose the bmopt.fullaxis in get_bmopts
  % Standard to use 2deg maps for window
  switch plotmap
    case 1
      bm_plotter(bmopt)
    case 2
      s = gen_stamp();
      mkdir('farmfiles/plotter')
      tmpfile = ['farmfiles/plotter/' expt '_' num2str(year) '_plotter_bm' ...
                 sprintf('%02i',bmnum) '_' s];
      cmd = 'bm_plotter(bmopt)';
      farmit(tmpfile,cmd,'var',{'bmopt'},'queue',...
             queues,'mem',memplot,'maxtime',maxtime);
  end
  
  
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function bmopt = get_bmopts(expt,year,realorsim,buddymap,winbuddymap)

bmopt.year = year;

% Select experiment-specific options
switch expt
  case 'bicep2'
    bmopt.expt = 'bicep2';
  case 'keck'
    bmopt.expt = 'keck';
    switch year
      case 2015
        bmopt.fullaxis = [-20 20 -15 25];
    end
    bmopt.run = ['k' num2str(year)];
  case 'bicep3'
    bmopt.expt = 'bicep3';
    switch year
      case 2015
        bmopt.run = 'b3r6';
      case 2016
        bmopt.run = 'b3r8';
        bmopt.fullaxis = [-210 -150 -10 30]; % CHANGE MEEEEE
      case {2017,2018}
        bmopt.run = 'b3r9';
    end
end

% Elevation resolution determines pixelization (2*el_res) -- standard
% maps are made with 0.05 deg steps, but this could change...move to
% above block later if needed!
bmopt.el_res = 0.05;

% Timestream / demod options
bmopt.notlastscan = 1;
bmopt.demodtype = 'square';
bmopt.timestreamdir = 'beammaps/timestream/';
bmopt.pairdiffsum = 0;
bmopt.beamparammatdir = 'beammaps/beamparams/';
bmopt.beamparammatfile = ['beamparams_' num2str(year)];

% Fullmap options
bmopt.component = 'cos';
bmopt.mapcoord = 'xpyp';
bmopt.fullmapdir = 'beammaps/maps/';
bmopt.xpypmax = 40;
bmopt.suffix = '';

% Winmap options
bmopt.winmapdir = 'beammaps/maps_windowed/';
bmopt.shift = 'ab_centroid';

% If we're dealing with buddy maps, use these directories instead
if winbuddymap > 0 || buddymap > 0
  bmopt.fullmapdir = 'beammaps/buddymaps/';
  bmopt.winmapdir = 'beammaps/buddymaps_windowed/';
end

% If we're using sim data, specify options/directories
if strcmp(realorsim,'sim')
  bmopt.btsimopt.mapopt = 'B3GaussianNoChopperSignalOnly';
  bmopt.btsimopt.refopt.type = 'square';
  bmopt.btsimopt.refopt.freq = 16;
  bmopt.btsimopt.refopt.phase = 0;
  bmopt.cleanref = 0;
  bmopt.timestreamdir = 'beammaps/bmtodsim/timestream/';
  bmopt.fullmapdir = 'beammaps/bmtodsim/maps/';
  bmopt.winmapdir = 'beammaps/bmtodsim/maps_windowed/';
elseif strcmp(realorsim,'real')
  bmopt.btsimopt = 0;
else
  error('realorsim should be ''real'' or ''sim''')
end

% Make html options
bmopt.author = 'TSG';
bmopt.pager = 1;

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [source_az source_el] = get_sourcepos(expt,year,number)
% For windowing in azel_ap, provide the apparent az/el to center on.
% This changes slightly with each movement of the mirror, 
% so for Keck with multiple
% mirror positions we need to know the run number

switch expt
  case 'bicep2'
  case 'keck'
    switch year
      case 2012
        % Updated 2017-09-21 TSG
	source_az = 0.92;
        source_el = 2.94;
      case 2013
        % Updated 2017-09-21 TSG
	source_az = 0.62;
        source_el = 3.05;
      case 2014
	% Updated 2017-09-21 TSG
	source_az = 0.23;
	source_el = 3.69;
      case 2015
	% Updated 2017-09-21 TSG
        if number <= 29
          source_az = 0.17;
          source_el = 3.93;
        else
          source_az = -0.02;
          source_el = 5.29;
        end
      case 2016
	% Updated 2017-09-21 TSG
	if number <= 29
	  source_az = -0.11;
	  source_el = 4.13;
	else
	  source_az = -0.17;
	  source_el = 4.80;
	end
      case 2017
	% Updated 2017-09-21 TSG
	if number <= 21
	  source_az = -0.15;
	  source_el = 4.83;
	else
	  source_az = -0.26;
	  source_el = 4.33;
	end
      case 2018
	% Back position with BSNS source
 	% Two sets of blocks were swapped in mirror installation,
 	% so as of now there is a ~7 deg offset in el. 
 	% Rough guess as of now.
 	if number <= 10
 	  source_az = -0.26;
 	  source_el = -3.5;
 	% Forwardest with BSNS. Still with ~7 deg offset.
 	elseif number <=22 
 	  source_az = -0.3;
 	  source_el = -6.9;
	% Forwardest with chopper, after mirror legs adjustment.
        % Legs were adjusted to compensate for incorrect block installtion.
	elseif number <=42
	  source_az = -0.25;
	  source_el = 4.25;
	% Back with chopper, with legs adjustment
        else
	  source_az = -0.25;
	  source_el = 6.7;
 	end
    end
  case 'bicep3'
    switch year
      case 2015
      case 2016
        % Updated 2017-09-21 TSG
	source_az = -177.6;
        source_el = 2.87;
      case 2017
        % Updated 2017-09-21 TSG
        source_az = -176.8;
        source_el = 3.23;
      case 2018 
	% Works for BSNS data
	source_az = -178.3;
        source_el = 2.5;
    end
end

return