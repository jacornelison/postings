function ffbm_plotCompMaps(compopt)
% function ffbm_plotCompMaps(compopt)
%
% Plot all component maps going into a composite in one large panel, plus
% the thumbnail plot.  Can plot with/without the final composite map.  Can
% be done pre/post-rotation, and with/without cuts.  
%
% Thumbnail plot is currently only r<2 deg, even if the component/composite
% files are larger. 
%
% Different behavior depending on which stage in the procedure we're at
% 'component': component maps in raw apparent az/el
% 'masked':    masked/centered maps in detector-centered apparent az/el
% 'rotated':   rotated maps in 'dk0' or x'/y'
%   All the above plot thumbnails as the number of component maps for
%   that detector
% 'composite': composite + thumbnail, and rotated maps
%
% RECALL that imagesc(array) plots such that axis 1 of the array is the
% up/down direction - but if you send in axis limits, it wants
% imagesc(ax2,ax1,array), i.e. the up/down direction limit is the SECOND
% argument!!!
%
% INPUTS (should all be sent in with compopt)
%
%  expt:           'keck','b3'
%  year:           keck: 2014, 2015
%                  b3:   2015
%
% OPTIONAL INPUTS 
%
%   plottype:        'component','masked','rotated','composite'
%                    'component' default
%   makefull:        make full panel plots of components (different behavior
%                    depending on whether composite exists)
%   makesmall:       make thumbnail plot (different behavior depending on
%                    whether composite exists)
%   makediff:        Plot A/B diff beam (only if plotting composites)
%                    1 default
%   subtightscale:   spacing for 'subtightplot' subfunction, default 0.06
%   dblim:           clim for peak-normalized dB plots, default [-30 0]
%   plotfitxy:       Plot fit az/el coordinates plus amplitude in
%                    component map panel (default 1)
%   coord:           Coordinate system for rotated/composite maps:
%                    'dk0' for inputs to reduc_makesim
%                    'xpyp' for normal visualization (x'/y')
%   onlychflag:      Only make thumbnails for pairs which pass chflags
%                    1 default
%   onlycommon_ab:   0: use all schedules passing cuts (default)
%                    1: only keep schedules where A/B were both measured
%   plotdir:         directory in which to save plots
%   componentdir:    directory from which to load component struct
%                    (default beammaps/maps_component)
%   componentfile:   name of file in componentdir to load
%                    (default ffbm_year_allcomp)
%   compositedir:    directory from which to load composite
%                    (default beammaps/composites)
%   compositefile:   name of file in compositedir to load
%                    (default ffbm_year_all)
%   applycuts:       0: plot all 
%                    1: plot only maps which pass cuts (default)
%   cutdir:          directory from which to load the cut structure
%                    (default beammaps/cuts)
%   cutfile:         name of cut file to load
%                    (default ffbm_cuts_year)
%   cutlist:         criteria on which to cut 
%                    (default all in cutfile)

% Parse compopt
expt = compopt.expt;
year = compopt.year;
if ~isfield(compopt,'componentdir')
  compopt.componentdir = 'beammaps/maps_component';
  componentdir = compopt.componentdir;
else
  componentdir = compopt.componentdir;
end
if ~isfield(compopt,'componentfile')
  compopt.componentfile = ['ffbm_' num2str(year) '_allcomp'...
	compopt.suffix];
  componentfile = compopt.componentfile;
else
  componentfile = compopt.componentfile;
end
if ~isfield(compopt,'compositedir')
  compopt.compositedir = 'beammaps/maps_composite';
  compositedir = compopt.compositedir;
else
  compositedir = compopt.compositedir;
end
if ~isfield(compopt,'compositefile')
  compopt.compositefile = ['ffbm_' num2str(year) '_all'...
	compopt.suffix];
  compositefile = compopt.compositefile;
else
  compositefile = compopt.compositefile;
end
if ~isfield(compopt,'applycuts')
  compopt.applycuts = 1;
  applycuts = compopt.applycuts;
else
  applycuts = compopt.applycuts;
end
if ~isfield(compopt,'cutdir')
  compopt.cutdir = 'beammaps/cuts';
  cutdir = compopt.cutdir;
else
  cutdir = compopt.cutdir;
end
if ~isfield(compopt,'cutfile')
  compopt.cutfile = ['ffbm_cuts_' num2str(year)];
  cutfile = compopt.cutfile;
else
  cutfile = compopt.cutfile;
end
if ~isfield(compopt,'cutlist')
  compopt.cutlist = {'mirror','notlight','nanfit','peak',...
      'sigma','ellip','hand'};
  cutlist = compopt.cutlist;
else
  cutlist = compopt.cutlist;
end
% Plotting options
if ~isfield(compopt,'plottype')
  compopt.plottype = 'component';
  plottype = compopt.plottype;
else
  plottype = compopt.plottype;
end
switch plottype
  case {'component','masked'}
    coord = 'azel_ap';
    onlychflag = 0;
  case 'rotated'
    if ~isfield(compopt,'coord')
      compopt.coord = 'xpyp';
      coord = compopt.coord;
    else
      coord = compopt.coord;
    end
    onlychflag = 0;
  case 'composite'
    if ~isfield(compopt,'makediff')
      compopt.makediff = 1;
      makediff = compopt.makediff;
    else
      makediff = compopt.makediff;
    end
    if ~isfield(compopt,'onlychflag')
      compopt.onlychflag = 1;
      onlychflag = compopt.onlychflag;
    else
      onlychflag = compopt.onlychflag;
    end
    if ~isfield(compopt,'coord')
      compopt.coord = 'xpyp';
      coord = compopt.coord;
    else
      coord = compopt.coord;
    end
end
if ~isfield(compopt,'subtightscale')
  compopt.subtightscale = 0.07;
  subtightscale = compopt.subtightscale;
else
  subtightscale = compopt.subtightscale;
end
if ~isfield(compopt,'dblim')
  compopt.dblim = [-30 0];
  dblim = compopt.dblim;
else
  dblim = compopt.dblim;
end
if ~isfield(compopt,'plotfitxy')
  compopt.plotfitxy = 1;
  plotfitxy = compopt.plotfitxy;
else
  plotfitxy = compopt.plotfitxy;
end
if ~isfield(compopt,'onlycommon_ab')
  compopt.onlycommon_ab = 0;
  onlycommon_ab = compopt.onlycommon_ab;
else
  onlycommon_ab = compopt.onlycommon_ab;
end
if ~isfield(compopt,'makefull')
  makefull = 1;
else
  makefull = compopt.makefull;
end
if ~isfield(compopt,'makesmall')
  makesmall = 1;
else
  makesmall = compopt.makesmall;
end
if ~isfield(compopt,'plotdir')
  compopt.plotdir = 'plots/';
  plotdir = compopt.plotdir;
else
  plotdir = compopt.plotdir;
end

% Get chflags if required
chflags = [];
if onlychflag
  chflags = get_default_chflags(expt,year);
end

[p ind] = get_array_info([num2str(year) '0201'],...
                         [],[],[],[],chflags);
[p0 ind0] = get_array_info([num2str(year) '0201']);

% Make directories
switch expt
  case 'b3'
    mkdir(plotdir);
    mkdir([plotdir '/small']);
    mkdir([plotdir '/full']);
  case 'keck'
    for ii = 0:4
      mkdir([plotdir '_rx' num2str(ii)]);
      mkdir([plotdir '_rx' num2str(ii) '/small']);
      mkdir([plotdir '_rx' num2str(ii) '/full']);
    end
end

% Load cuts if needed
if applycuts
  load([cutdir '/' cutfile]);
  cut = cuts.cuts;
end

% Load up component map file
load([componentdir '/' componentfile]);
disp(['Plotting components ' componentdir '/' componentfile]);

% Load up composites if we have them
switch plottype
  case 'composite'
    load([compositedir '/' compositefile]);
    disp(['Plotting composites ' compositedir '/' compositefile]);
end

for ii = 1:length(ind.la)
  
  if mod(ii,10) == 0
    disp(['Plotting pair ' int2str(ii) ' / ' int2str(length(ind.la))]);
  end
  
  % Make a total cut mask
  totmask_a = ones(length(comp.bm.number),1);
  totmask_b = ones(length(comp.bm.number),1);
  
  if applycuts
    for jj = 1:length(comp.bm.number);
      a = structcut(cut(jj),ind.la(ii)); % all cuts for this detector
      b = structcut(cut(jj),ind.lb(ii)); 
      totcut_a = 0;
      totcut_b = 0;
      for kk = 1:length(cutlist)
	totcut_a = totcut_a + getfield(a,cutlist{kk});
	totcut_b = totcut_b + getfield(b,cutlist{kk});
      end
      % Handcut == 1 means the fit may suck, but it's probably okay for
      % compositing
      if totcut_a == 0 | a.hand == 1 
	totmask_a(jj) = 0;
      end
      if totcut_b == 0 | b.hand == 1
	totmask_b(jj) = 0;
      end
    end
  end
  
  % Indices of good maps
  plotinds_a = find(totmask_a == 0);
  plotinds_b = find(totmask_b == 0);
  
  % If we want to keep only schedules where A/B were both measured
  if onlycommon_ab
    plotinds_a = intersect(plotinds_a,plotinds_b);
    plotinds_b = plotinds_a;
  end

  ro = int2str(p.det_row(ind.la(ii)));
  co = int2str(p.det_col(ind.la(ii)));
  ti = int2str(p.tile(ind.la(ii)));
  rx = int2str(p.rx(ind.la(ii)));
  mce = int2str(p.mce(ind.la(ii)));
  
  % Plotting directories - keep by rx for Keck
  switch expt
    case {'b2','b3'}
      subplotdir = plotdir;
    case 'keck'
      subplotdir = [plotdir '_rx' num2str(p.rx(ind.la(ii)))];
  end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % COMPOSITE MAPS + COMPONENTS
  
  switch plottype
    case 'composite'
      if makefull
	
	% This works for up to 23 component maps
	sp1 = 4;
	sp2 = 9;
	plotpos = [4,5,6,7,8,9,...
	      13,14,15,16,17,18,...
	      22,23,24,25,26,27,...
	      31,32,33,34,35,36];
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% A detector
	fh = figure(1);
	set(fh,'Position',[0 0 1600 900],'Visible','off')
	clf
	colormap('jet')

	gcp = int2str(p.gcp(ind.la(ii)));
	
	% Peak normalize full map, and normalize std by same value
	map = composite.map_med{ind.la(ii)};
	norm = composite.map_fit{ind.la(ii)}(1); %nanmax(map(:)); 
	map = map./norm;
	%map(map < 0) = 0; % for log plots
	std = composite.map_std{ind.la(ii)};
	std = std./norm;
	
        ax1 = comp.ad.t_val_deg{1}; % See documentation in
        ax2 = comp.ad.t_val_deg{2}; % ffbm_rotatemaps
        
	% Main panel composite
	subplot(2,3,1)
	imagescnan(ax2,ax1,10*log10(abs(map)));
	title(['GCP ' gcp ' composite'])
	set(gca,'YDir','normal')
	caxis(dblim)
	colorbar
	axis square
        switch coord
          case 'dk0'
            xlabel('Azimuth')
            ylabel('Elevation')
          case 'xpyp'
            set(gca,'XDir','reverse')
            xlabel('y prime')
            ylabel('x prime') 
        end
	
	% Main panel rms
	subplot(2,3,4)
	imagescnan(ax2,ax1,std)
	title(['GCP ' gcp ' RMS'])
	set(gca,'YDir','normal')
	caxis([0 0.05])
	axis square
	colorbar
        switch coord
          case 'dk0'
            xlabel('Azimuth')
            ylabel('Elevation')
          case 'xpyp'
            set(gca,'XDir','reverse')
            xlabel('y prime')
            ylabel('x prime') 
        end

	% Up to 24 component plots
	for jj = 1:min([length(plotinds_a),23])
          maptoplot = comp.map{ind.la(ii)}.component{plotinds_a(jj)};
	  %maptoplot(maptoplot < 0) = 0;
	  dk = comp.bm.dk(plotinds_a(jj));
	  
	  clim_sum = get_clim_sum(expt,p,ind.la,ii,...
                                       comp.bm.filename{plotinds_a(jj)});
          mA = comp.map{ind.la(ii)}.fit{plotinds_a(jj)}(1);
	  
	  subplot(sp1,sp2,plotpos(jj))
          if isnan(mA) 
            imagescnan(ax2,ax1,10*log10(abs(maptoplot)));
            prettify_large([10*log10(clim_sum(2)) - range(dblim),...
                            10*log10(clim_sum(2))],coord);
          else
            imagescnan(ax2,ax1,10*log10(abs(maptoplot./mA)));
            prettify_large(dblim,coord);
          end
	  title(['Map ' num2str(plotinds_a(jj)) ' dk ' num2str(dk)]);

        end
	
        % Add a legend in the last panel
        switch expt
          case {'b2','keck'}
            legendstr1 = ['Rx ' rx ' GCP ' gcp ];
          case 'b3'
            legendstr1 = ['MCE' mce ' GCP ' gcp ];
        end
	legendstr2 = ['Tile ' ti ' Row ' ro ' Col ' co ];
	subplot(sp1,sp2,plotpos(end))
	text(0,0.8,legendstr1,'FontSize',16);
	text(0,0.5,legendstr2,'FontSize',16);
	axis off
	
	% Save image
	img_str = [subplotdir '/full/det_row_' ro '_col_' co '_tile_' ti ...
	      '_' p.pol{ind.la(ii)}];
	mkpng(img_str);

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% B detector
	fh = figure(1);
	set(fh,'Position',[0 0 1600 900],'Visible','off')
	clf
	colormap('jet')

	gcp = int2str(p.gcp(ind.lb(ii)));
	
	% Peak normalize full map, and normalize std by same value
	map = composite.map_med{ind.lb(ii)};
	norm = composite.map_fit{ind.lb(ii)}(1); %nanmax(map(:)); 
	map = map./norm;
	%map(map < 0) = 0; % for log plots
	std = composite.map_std{ind.lb(ii)};
	std = std./norm;

        ax1 = comp.ad.t_val_deg{1}; % See documentation in
        ax2 = comp.ad.t_val_deg{2}; % ffbm_rotatemaps

        % Main panel composite
	subplot(2,3,1)
	imagescnan(ax2,ax1,10*log10(abs(map)));
	title(['GCP ' gcp ' composite'])
	set(gca,'YDir','normal')
	caxis(dblim)
	colorbar
	axis square
        switch coord
          case 'dk0'
            xlabel('Azimuth')
            ylabel('Elevation')
          case 'xpyp'
            set(gca,'XDir','reverse')
            xlabel('y prime')
            ylabel('x prime') 
        end

        % Main panel rms
	subplot(2,3,4)
	imagescnan(ax2,ax1,std)
	title(['GCP ' gcp ' RMS'])
	set(gca,'YDir','normal')
	caxis([0 0.05])
	axis square
	colorbar
        switch coord
          case 'dk0'
            xlabel('Azimuth')
            ylabel('Elevation')
          case 'xpyp'
            set(gca,'XDir','reverse')
            xlabel('y prime')
            ylabel('x prime') 
        end
	
	% Up to 24 component plots
	for jj = 1:min([length(plotinds_b),23])
          maptoplot = comp.map{ind.lb(ii)}.component{plotinds_b(jj)};
	  %maptoplot(maptoplot < 0) = 0;
	  dk = comp.bm.dk(plotinds_b(jj));
	  
	  clim_sum = get_clim_sum(expt,p,ind.lb,ii,...
                                       comp.bm.filename{plotinds_b(jj)});
          mB = comp.map{ind.lb(ii)}.fit{plotinds_b(jj)}(1);
	  
	  subplot(sp1,sp2,plotpos(jj))
          if isnan(mB) 
            imagescnan(ax2,ax1,10*log10(abs(maptoplot)));
            prettify_large([10*log10(clim_sum(2)) - range(dblim),...
                            10*log10(clim_sum(2))],coord);
          else
            imagescnan(ax2,ax1,10*log10(abs(maptoplot./mB)));
            prettify_large(dblim,coord);
          end
	  title(['Map ' num2str(plotinds_b(jj)) ' dk ' num2str(dk)]);

        end
	% Add a legend in the last panel
        switch expt
          case {'b2','keck'}
            legendstr1 = ['Rx ' rx ' GCP ' gcp ];
          case 'b3'
            legendstr1 = ['MCE' mce ' GCP ' gcp ];
        end
	legendstr2 = ['Tile ' ti ' Row ' ro ' Col ' co ];
	subplot(sp1,sp2,plotpos(end))
	text(0,0.8,legendstr1,'FontSize',16);
	text(0,0.5,legendstr2,'FontSize',16);
	axis off
	
	% Save image
	img_str = [subplotdir '/full/det_row_' ro '_col_' co '_tile_' ti ...
	      '_' p.pol{ind.lb(ii)}];
	mkpng(img_str);
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Diff
	if makediff
	  fh = figure(1);
          set(fh,'Position',[0 0 1600 900],'Visible','off')
	  clf
	  colormap('jet')

	  % Normalized diff between A/B
	  map = composite.map_med{ind.la(ii)}./...
                composite.map_fit{ind.la(ii)}(1) ...
	      - composite.map_med{ind.lb(ii)}./...
                composite.map_fit{ind.lb(ii)}(1);
	  
          % Get the fit sigma
          sig = nanmedian([composite.map_fit{ind.la(ii)}(4),...
                           composite.map_fit{ind.la(ii)}(5),...
                           composite.map_fit{ind.lb(ii)}(4),...
                           composite.map_fit{ind.lb(ii)}(5)]);
          
          % Now deproject the diff map
          [c T] = ffbm_deprojmap(inpaint_nans(map),ad,sig);
          dp1100 = map - c(1).*T{1} - c(2).*T{2} - c(3).*T{3};
          dp1101 = dp1100 - c(5).*T{5} - c(6).*T{6};
          dp1111 = dp1101 - c(4).*T{4};
          
          % Standard difference
          subplot(2,3,1)
	  imagescnan(ax2,ax1,map);
	  title('0000')
	  set(gca,'YDir','normal')
	  caxis([-0.1 0.1])
	  colorbar
	  axis square
          switch coord
            case 'dk0'
              xlabel('Azimuth')
              ylabel('Elevation')
            case 'xpyp'
              set(gca,'XDir','reverse')
              xlabel('y prime')
              ylabel('x prime') 
          end

          % Deprojecting relgain/diff pointing
          subplot(2,3,2)
          imagescnan(ax2,ax1,dp1100);
	  title('1100: Relgain/diff pointing')
	  set(gca,'YDir','normal')
	  caxis([-0.025 0.025])
	  colorbar
	  axis square
          switch coord
            case 'dk0'
              xlabel('Azimuth')
              ylabel('Elevation')
            case 'xpyp'
              set(gca,'XDir','reverse')
              xlabel('y prime')
              ylabel('x prime') 
          end
          
          % Deprojecting relgain/diff pointing/diffellip
          subplot(2,3,4)
          imagescnan(ax2,ax1,dp1100);
	  title('1101: adding diff ellip')
	  set(gca,'YDir','normal')
	  caxis([-0.025 0.025])
	  colorbar
	  axis square
          switch coord
            case 'dk0'
              xlabel('Azimuth')
              ylabel('Elevation')
            case 'xpyp'
              set(gca,'XDir','reverse')
              xlabel('y prime')
              ylabel('x prime') 
          end
          
          % Deprojecting all
          subplot(2,3,5)
          imagescnan(ax2,ax1,dp1111);
	  title('1111: adding diff bw')
	  set(gca,'YDir','normal')
	  caxis([-0.025 0.025])
	  colorbar
	  axis square
          switch coord
            case 'dk0'
              xlabel('Azimuth')
              ylabel('Elevation')
            case 'xpyp'
              set(gca,'XDir','reverse')
              xlabel('y prime')
              ylabel('x prime') 
          end
          
          % Add a legend in the last panel
          switch expt
            case {'b2','keck'}
              legendstr1 = ['Rx ' rx ' Diff'];
            case 'b3'
              legendstr1 = ['MCE' mce ' Diff'];
          end
          legendstr2 = ['Tile ' ti ' Row ' ro ' Col ' co ];
          subplot(sp1,sp2,plotpos(end))
          text(0,0.8,legendstr1,'FontSize',16);
          text(0,0.5,legendstr2,'FontSize',16);
          axis off
          
	  img_str = [subplotdir '/full/det_row_' ro '_col_' co '_tile_' ti ...
		'_dif'];
	  mkpng(img_str);
	  
	end
	
      end % makelarge

      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % Thumbnail plot of composite
      if makesmall

	% Plot Xs through maps not going into CMB analysis
	% ind0 is before channel cuts, ind is after
	if onlychflag
	  
	  notrgl = 0;
	  cutbych = 0;
	  mapsnogood = 0;
	  
	  % Case: not on RGL list before or after ch cuts
	  if ~ismember(ind.la(ii),ind0.rgl) & ~ismember(ind.la(ii),ind.rgl)
	    notrgl = 1;
	  end
	  % Case: on RGL list before ch cuts but not after
	  if ismember(ind.la(ii),ind0.rgl) & ~ismember(ind.la(ii),ind.rgl)
	    cutbych = 1;
	  end
	  % Case: on RGL list after ch cuts, but one or both of A/B have not
	  % maps
	  if ismember(ind.la(ii),ind.rgl) 
	    if isempty(plotinds_a) | isempty(plotinds_b)
	      mapsnogood = 1;
	    end
	  end
          
        else
            
          notrgl = 0;
          cutbych = 0;
          mapsnogood = 0;
          
	end
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% A detector
	fh = figure(1);
	set(fh,'Position',[10 10 670 500],'Visible','off')
	clf
	colormap('jet')

	map = composite.map_med{ind.la(ii)};
	norm = composite.map_fit{ind.la(ii)}(1); %nanmax(map(:)); 
	map = map./norm;
	%map(map < 0) = 0; % for log plots
	
	% Single main image of composite
	imagescnan(ax2,ax1,10*log10(abs(map)))
	caxis(dblim)
	if notrgl
	  hold on
	  plot([-2 2],[-2 2],'r','linewidth',20)
	  plot([-2 2],[2 -2],'r','linewidth',20)
	end
	if cutbych
	  hold on
	  plot([-2 2],[-2 2],'y','linewidth',20)
	  plot([-2 2],[2 -2],'y','linewidth',20)
	end
	if mapsnogood
	  hold on
	  plot([-2 2],[-2 2],'m','linewidth',20)
	  plot([-2 2],[2 -2],'m','linewidth',20)
	end
	axis tight
	axis equal
        axis([-2 2 -2 2])
	set(gca,'YDir','normal')
        if strcmp(coord,'xpyp')
          set(gca,'XDir','reverse')
        end
	set(gca,'xtick',[],'ytick',[])
	
	% Save image
	img_str = [subplotdir '/small/det_row_' ro '_col_' co '_tile_' ti ...
	      '_' p.pol{ind.la(ii)} '_sm'];
	mkpng(img_str);
	% 400x400 window starting at (145,40)
	system(['mogrify -crop ''400x400+145+40'' ' subplotdir ...
                '/small/det_row_' ro '_col_' co ...
                '_tile_' ti '_' p.pol{ind.la(ii)} '*']);
	% Resize to 75 pixels (400 * 0.18) 
	system(['mogrify -resize ''18%x18%'' ' subplotdir ...
                '/small/det_row_' ro '_col_' co ...
                '_tile_' ti '_' p.pol{ind.la(ii)} '*']);
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% B detector

	fh = figure(1);
	set(fh,'Position',[10 10 670 500],'Visible','off')
	clf
	colormap('jet')
	
	map = composite.map_med{ind.lb(ii)};
	norm = composite.map_fit{ind.lb(ii)}(1); %nanmax(map(:)); 
	map = map/norm;
	%map(map < 0) = 0; % for log plots
	
	% Single main image of composite
	imagescnan(ax2,ax1,10*log10(abs(map)))
	caxis(dblim)
	if notrgl
	  hold on
	  plot([-2 2],[-2 2],'r','linewidth',20)
	  plot([-2 2],[2 -2],'r','linewidth',20)
	end
	if cutbych
	  hold on
	  plot([-2 2],[-2 2],'y','linewidth',20)
	  plot([-2 2],[2 -2],'y','linewidth',20)
	end
	if mapsnogood
	  hold on
	  plot([-2 2],[-2 2],'m','linewidth',20)
	  plot([-2 2],[2 -2],'m','linewidth',20)
	end

	axis tight
	axis equal
        axis([-2 2 -2 2])
	set(gca,'YDir','normal')
        if strcmp(coord,'xpyp')
          set(gca,'XDir','reverse')
        end
	set(gca,'xtick',[],'ytick',[])
	
	% Save image
	img_str = [subplotdir '/small/det_row_' ro '_col_' co '_tile_' ti ...
	      '_' p.pol{ind.lb(ii)} '_sm'];
	mkpng(img_str);
	% 400x400 window starting at (145,40)
	system(['mogrify -crop ''400x400+145+40'' ' subplotdir ...
                '/small/det_row_' ro '_col_' co ...
                '_tile_' ti '_' p.pol{ind.lb(ii)} '*']);
	% Resize to 75 pixels (400 * 0.18) 
	system(['mogrify -resize ''18%x18%'' ' subplotdir ...
                '/small/det_row_' ro '_col_' co ...
                '_tile_' ti '_' p.pol{ind.lb(ii)} '*']);
	
	if makediff

	  fh = figure(1);
	  set(fh,'Position',[10 10 670 500],'Visible','off')
	  clf
	  colormap('jet')

	  % Normalized diff between A/B
	  map = composite.map_med{ind.la(ii)}./...
                composite.map_fit{ind.la(ii)}(1) ...
	      - composite.map_med{ind.lb(ii)}./...
                composite.map_fit{ind.lb(ii)}(1);
	  
	  % Single main image of composite
	  imagescnan(ax2,ax1,map);
	  caxis([-0.1 0.1])
	  if notrgl
	    hold on
	    plot([-2 2],[-2 2],'r','linewidth',20)
	    plot([-2 2],[2 -2],'r','linewidth',20)
	  end
	  if cutbych
	    hold on
	    plot([-2 2],[-2 2],'y','linewidth',20)
	    plot([-2 2],[2 -2],'y','linewidth',20)
	  end
	  if mapsnogood
	    hold on
	    plot([-2 2],[-2 2],'m','linewidth',20)
	    plot([-2 2],[2 -2],'m','linewidth',20)
	  end
	  axis tight
	  axis equal
          axis([-2 2 -2 2])
	  set(gca,'YDir','normal')
          if strcmp(coord,'xpyp')
            set(gca,'XDir','reverse')
          end
	  set(gca,'xtick',[],'ytick',[])

	  img_str = [subplotdir '/small/det_row_' ro '_col_' co '_tile_' ti ...
		'_dif_sm'];
	  mkpng(img_str);
	  % 400x400 window starting at (145,40)
	  system(['mogrify -crop ''400x400+145+40'' ' subplotdir ...
                  '/small/det_row_' ro '_col_' co '_tile_' ti '_dif*']);
	  % Resize to 75 pixels (400 * 0.18) 
	  system(['mogrify -resize ''18%x18%'' ' subplotdir ...
                  '/small/det_row_' ro '_col_' co '_tile_' ti '_dif*']);

	end
	
      end % makesmall
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
    case {'component','masked','rotated'} % plottype
    
      if makefull
	
	% This works for up to 23 component maps
	sp1 = 4;
	sp2 = 6;
	plotpos = 1:(sp1*sp2);

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% A Detector
      	fh = figure(1);
	set(fh,'Position',[0 0 1600 900],'Visible','off')
	clf
	colormap('jet')

	gcp = int2str(p.gcp(ind.la(ii)));
	
	% Up to 24 component plots
        for jj = 1:min([length(plotinds_a),24])

          maptoplot = comp.map{ind.la(ii)}.component{plotinds_a(jj)};
	  dk = comp.bm.dk(plotinds_a(jj));
	  clim_sum = get_clim_sum(expt,p,ind.la,ii,...
                                       comp.bm.filename{plotinds_a(jj)});
	  fit(jj,:) = comp.map{ind.la(ii)}.fit{plotinds_a(jj)};
          
          switch plottype
            case 'component'
              ax1 = comp.el_ap{plotinds_a(jj)};
              ax2 = comp.az_ap{plotinds_a(jj)};
            case 'masked'
              ax1 = comp.ad.t_val_deg{2}; % See documentation in
              ax2 = comp.ad.t_val_deg{1}; % ffbm_maskground
            case 'rotated'
              ax1 = comp.ad.t_val_deg{1}; % See documentation in
              ax2 = comp.ad.t_val_deg{2}; % ffbm_rotatemaps
          end

	  subtightplot(sp1,sp2,plotpos(jj),subtightscale)
          if isnan(fit(jj,1)) 
            imagescnan(ax2,ax1,10*log10(abs(maptoplot)));
            prettify_large([10*log10(clim_sum(2)) - range(dblim),...
                            10*log10(clim_sum(2))],coord);
          else
            imagescnan(ax2,ax1,10*log10(abs(maptoplot./fit(jj,1))));
            prettify_large(dblim,coord);
          end
	  title(['Map ' num2str(plotinds_a(jj)) ' dk ' num2str(dk)]);
          
          if jj == 1
            switch coord
              case 'azel_ap'
                xlabel('Az apparent')
                ylabel('El apparent')
              case 'dk0'
                xlabel('Azimuth')
                ylabel('Elevation')
              case 'xpyp'
                xlabel('y prime')
                ylabel('x prime')
            end
          end
	end % plotinds
    
	% Add a legend in the last panel
	legendstr1 = ['Rx ' rx ' GCP ' gcp];
	legendstr2 = ['Tile ' ti ' Row ' ro ' Col ' co ];
	subtightplot(sp1,sp2,plotpos(end),subtightscale)
	text(0,0.8,legendstr1,'FontSize',16);
	text(0,0.5,legendstr2,'FontSize',16);
	axis off
	
        if strcmp(plottype,'component') & exist('fit','var')
          subtightplot(sp1,sp2,plotpos(end-1),subtightscale)
          scatter(fit(:,2),fit(:,3),[],fit(:,1),'filled');
          axis equal;
          cb = colorbar();
          xlabel(cb,'Amp');
          xlabel('Fit az');
          ylabel('Fit el');
        end
        clear fit
        
	% Save image
	img_str = [subplotdir '/full/det_row_' ro '_col_' co '_tile_' ti ...
	      '_' p.pol{ind.la(ii)}];
	mkpng(img_str);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% B Detector
      	fh = figure(1);
	set(fh,'Position',[0 0 1600 900],'Visible','off')
	clf
	colormap('jet')

	gcp = int2str(p.gcp(ind.lb(ii)));
	
	% Up to 24 component plots
        for jj = 1:min([length(plotinds_b),24])

          maptoplot = comp.map{ind.lb(ii)}.component{plotinds_b(jj)};
	  dk = comp.bm.dk(plotinds_b(jj));
	  clim_sum = get_clim_sum(expt,p,ind.lb,ii,...
                                       comp.bm.filename{plotinds_b(jj)});
	  fit(jj,:) = comp.map{ind.lb(ii)}.fit{plotinds_b(jj)};
          
          switch plottype
            case 'component'
              ax1 = comp.el_ap{plotinds_b(jj)};
              ax2 = comp.az_ap{plotinds_b(jj)};
            case 'masked'
              ax1 = comp.ad.t_val_deg{2}; % See documentation in
              ax2 = comp.ad.t_val_deg{1}; % ffbm_maskground
            case 'rotated'
              ax1 = comp.ad.t_val_deg{1}; % See documentation in
              ax2 = comp.ad.t_val_deg{2}; % ffbm_rotatemaps
          end

	  subtightplot(sp1,sp2,plotpos(jj),subtightscale)
          if isnan(fit(jj,1)) 
            imagescnan(ax2,ax1,10*log10(abs(maptoplot)));
            prettify_large([10*log10(clim_sum(2)) - range(dblim),...
                            10*log10(clim_sum(2))],coord);
          else
            imagescnan(ax2,ax1,10*log10(abs(maptoplot./fit(jj,1))));
            prettify_large(dblim,coord);
          end
	  title(['Map ' num2str(plotinds_b(jj)) ' dk ' num2str(dk)]);
          
          if jj == 1
            switch coord
              case 'azel_ap'
                xlabel('Az apparent')
                ylabel('El apparent')
              case 'dk0'
                xlabel('Azimuth')
                ylabel('Elevation')
              case 'xpyp'
                xlabel('y prime')
                ylabel('x prime')
            end
          end
	end % plotinds
    
	% Add a legend in the last panel
	legendstr1 = ['Rx ' rx ' GCP ' gcp];
	legendstr2 = ['Tile ' ti ' Row ' ro ' Col ' co ];
	subtightplot(sp1,sp2,plotpos(end),subtightscale)
	text(0,0.8,legendstr1,'FontSize',16);
	text(0,0.5,legendstr2,'FontSize',16);
	axis off
	
        if strcmp(plottype,'component') & exist('fit','var')
          subtightplot(sp1,sp2,plotpos(end-1),subtightscale)
          scatter(fit(:,2),fit(:,3),[],fit(:,1),'filled');
          axis equal;
          cb = colorbar();
          xlabel(cb,'Amp');
          xlabel('Fit az');
          ylabel('Fit el');
        end
        clear fit
        
	% Save image
	img_str = [subplotdir '/full/det_row_' ro '_col_' co '_tile_' ti ...
	      '_' p.pol{ind.lb(ii)}];
	mkpng(img_str);

      end % makefull
       
      if makesmall

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% A Detector
	fh = figure(1);
	set(fh,'Position',[10 10 670 500],'Visible','off')
	clf
	colormap('jet')
    
	% Dumb block of solid color
        if length(plotinds_a) > 0
          imagescnan(ones(size(maptoplot))*length(plotinds_a));
        else
          imagescnan(NaN(size(maptoplot)));
        end
	prettify_small_keck([0 24]);
	
	% Save image
	img_str = [subplotdir '/small/det_row_' ro '_col_' co '_tile_' ti ...
	      '_' p.pol{ind.la(ii)} '_sm'];
	%print('-dpng',img_str);
	mkpng(img_str)
	% Crop - this is good enough for a solid color!
	system(['mogrify -crop ''50x50+250+200'' ' subplotdir ...
                '/small/det_row_' ro '_col_' co ...
                '_tile_' ti '_' p.pol{ind.la(ii)} '*']);
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% B Detector
	fh = figure(1);
	set(fh,'Position',[10 10 670 500],'Visible','off')
	clf
	colormap('jet')
    
	% Dumb block of solid color
        if length(plotinds_b) > 0
          imagescnan(ones(size(maptoplot))*length(plotinds_b));
        else
          imagescnan(NaN(size(maptoplot)));
        end
	prettify_small_keck([0 24]);
	
	% Save image
	img_str = [subplotdir '/small/det_row_' ro '_col_' co '_tile_' ti ...
	      '_' p.pol{ind.lb(ii)} '_sm'];
	%print('-dpng',img_str);
	mkpng(img_str)
	% Crop - this is good enough for a solid color!
	system(['mogrify -crop ''50x50+250+200'' ' subplotdir ...
                '/small/det_row_' ro '_col_' co ...
                '_tile_' ti '_' p.pol{ind.lb(ii)} '*']);
      end  % makesmall
  
  end % switch plottype
    
end % ind.la
  
return % Main function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function clim_sum = get_clim_sum(expt,p,goodind,ii,filename)
% Color axes are unfortunately experiment and frequency-specific, so choose
% them here

bmtime = date2mjd(str2num(filename(1:4)),str2num(filename(5:6)),...
    str2num(filename(7:8)),str2num(filename(10:11)),...
    str2num(filename(12:13)),str2num(filename(14:15)));

switch expt
  case 'b2'
    clim_sum = [0 300];
  case 'keck'
    if bmtime < date2mjd(2016)
      switch p.band(goodind(ii))
        case {100,220}
          clim_sum = [0 200];
        case 150
          clim_sum = [0 300];
      end
    else
      clim_sum = [0 800];
    end
  case 'b3'
    % Before 2015-03-17 we had the 24" and after the 18"
    % For 2016 24"
    if bmtime < date2mjd(2015,03,17,00,00,00)
      clim_sum = [0 200];
    else
      clim_sum = [0 100];
    end
end

return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function prettify_large(cax,coord)

set(gca,'YDir','normal')
switch coord
  case 'xpyp'
    set(gca,'XDir','reverse')
end
caxis(cax)
axis equal
axis tight

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function prettify_small_keck(cax)

caxis(cax)
axis equal
set(gca,'xtick',[],'ytick',[])
set(gca,'YDir','normal')

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% OLD CODE - useful for reference

%{
% Set up plotting - sp1 = number of plots down, sp2 across
if length(qq) <= 12
  sp1 = 3;
  sp2 = 6;
  plotpos = [3,4,5,6,...
	     9,10,11,12,...
	     15,16,17,18];
elseif length(qq) > 12 && length(qq) <= 16
  sp1 = 4;
  sp2 = 6;
  plotpos = [3,4,5,6,...
	     9,10,11,12,...
	     15,16,17,18,...
	     21,22,23,24];
elseif length(qq) > 16 && length(qq) <= 20
  sp1 = 4;
  sp2 = 7;
  plotpos = [3,4,5,6,7,...
	     10,11,12,13,14,...
	     17,18,19,20,21,...
	     24,25,26,27,28];
elseif length(qq) > 20 
  sp1 = 4;
  sp2 = 8;
  plotpos = [3,4,5,6,7,8,...
	     11,12,13,14,15,16,...
	     19,20,21,22,23,24,...
	     27,28,29,30,31,32];
end

clim_sum = [0 1];
clim_sum2 = [-4 0];
clim_std = [0 0.05];
clim_dif = [-0.1 0.1];

markch = 1;


if ~skip_large
  % Plot all component maps
  figure('Position',[0 0 1600 900],'Visible','off')
  colormap('jet')
  
  for ii = 1:length(ind.a)
    marking = 0;
    if markch
      if isempty(intersect(ind.a(ii),ind.rgla))
	marking = 1;
      end
    end
    
    ii
    ro = int2str(pp.det_row(ind.a(ii)));
    co = int2str(pp.det_col(ind.a(ii)));
    ti = int2str(pp.tile(ind.a(ii)));
    gcpA = int2str(pp.gcp(ind.a(ii)));
    gcpB = int2str(pp.gcp(ind.b(ii)));
    
    A = map(:,:,ind.a(ii));
    B = map(:,:,ind.b(ii));
    Astd = map_std(:,:,ind.a(ii));
    Bstd = map_std(:,:,ind.b(ii));
    
    %plot A
    clf
    %  subplot(3,6,[7 14])
    %  imagesc(x_bin,y_bin,log10(abs(A)))
    %  rectangle('Position',[-4 -4 8 8],'LineWidth',2);
    %  title(['Det row ' ro ', col ' co ', tile ' ti ',  GCP ' gcpA ...
    %	', pol A, log10 scale'])
    %  set(gca,'YDir','normal')
    %  caxis(clim_sum2)
    %  colorbar
    %  axis image
    %  xlabel('x'', Deg')
    %  ylabel('y'', Deg')
    
    subplot(2,3,1)
    imagesc(x_bin,y_bin,log10(abs(A)))
    axis([-4 4 -4 4])
    title(['gcp' gcpA])
    set(gca,'YDir','normal')
    caxis(clim_sum2)
    colorbar
    axis square
    xlabel('Deg')
    ylabel('Deg')
    if marking
      hold on
      plot([-5 5],[-5 5],'k','linewidth',7)
      plot([-5 5],[5 -5],'k','linewidth',7)
    end
    
    subplot(2,3,4)
    imagesc(x_bin,y_bin,Astd)
    axis([-4 4 -4 4])
    title(['gcp' gcpA ' rms'])
    set(gca,'YDir','normal')
    caxis(clim_std)
    axis square
    colorbar
    xlabel('Deg')
    ylabel('Deg')    
    
    for jj = 1:length(qq)
      if jj < 24
	
	% Figure out if we cut this component
	handcut = 0;
	if cutmaps
	  if mapstocut(jj,ind.a(ii)) == 1
	    handcut = 1;
	  end
	end
	poscut = 0;
	bmcut = 0;
	if strcmp(experiment,'keck')
	  if isnan(qq{jj}.poscut(ind.a(ii)))
	    poscut = 1;
	  elseif isnan(qq{jj}.badmapcut(ind.a(ii)))
	    bmcut = 1;
	  end
	end
	subplot(sp1,sp2,plotpos(jj))
	amap = qq{jj}.map(:,:,ind.a(ii));
	imagesc(x_bin,y_bin,log10(abs(amap)))
	axis([-4 4 -4 4])
	title(['gcp' gcpA ', dk' num2str(qq{jj}.dk) ', sc' num2str(qq{jj}.sc)])
	set(gca,'YDir','normal')
	caxis(clim_sum2)
	axis square
	xlabel('Deg')
	ylabel('Deg'); 
	if handcut == 1
	  hold on
	  plot([-5 5],[-5 5],'g','linewidth',7)
	  plot([-5 5],[5 -5],'g','linewidth',7)
	end
	if poscut == 1
	  hold on
	  plot([-5 5],[-5 5],'r','linewidth',7)
	  plot([-5 5],[5 -5],'r','linewidth',7)
	end
	if bmcut == 1
	  hold on
	  plot([-5 5],[-5 5],'y','linewidth',7)
	  plot([-5 5],[5 -5],'y','linewidth',7)
	end
      end
    end
    
    img_str=[subdir '/full/det_row_' ro '_col_' co '_tile_' ti '_A'];
    print('-dpng', img_str);
    
    %B    
    clf
    %  subplot(3,6,[7 14])
    %  imagesc(x_bin,y_bin,log10(abs(B)))
    %  rectangle('Position',[-4 -4 8 8],'LineWidth',2);
    %  title(['Det row ' ro ', col ' co ', tile ' ti ',  GCP ' gcpB ...
    %	', pol A, log10 scale'])
    %  set(gca,'YDir','normal')
    %  caxis(clim_sum2)
    %  colorbar
    %  axis image
    %  xlabel('x'', Deg')
    %  ylabel('y'', Deg')
    
    subplot(2,3,1)
    imagesc(x_bin,y_bin,log10(abs(B)))
    axis([-4 4 -4 4])
    title(['gcp' gcpB])
    set(gca,'YDir','normal')
    caxis(clim_sum2)
    colorbar
    axis square
    xlabel('Deg')
    ylabel('Deg')
    if marking
      hold on
      plot([-5 5],[-5 5],'k','linewidth',7)
      plot([-5 5],[5 -5],'k','linewidth',7)
    end
    
    subplot(2,3,4)
    imagesc(x_bin,y_bin,Bstd)
    axis([-4 4 -4 4])
    title(['gcp' gcpB ' rms'])
    set(gca,'YDir','normal')
    caxis(clim_std)
    colorbar
    axis square
    xlabel('Deg')
    ylabel('Deg')
    
    for jj=1:length(qq)
      if jj < 24
	% Figure out if we cut this component
	handcut = 0;
	if cutmaps
	  if mapstocut(jj,ind.a(ii)) == 1
	    handcut = 1;
	  end
	end
	poscut = 0;
	bmcut = 0;
	if strcmp(experiment,'keck')
	  if isnan(qq{jj}.poscut(ind.a(ii)))
	    poscut = 1;
	  elseif isnan(qq{jj}.badmapcut(ind.a(ii)))
	    bmcut = 1;
	  end
	end
	subplot(sp1,sp2,plotpos(jj))
	amap = qq{jj}.map(:,:,ind.b(ii));
	imagesc(x_bin,y_bin,log10(abs(amap)))
	axis([-4 4 -4 4])
	title(['gcp' gcpB ', dk' num2str(qq{jj}.dk) ', sc' num2str(qq{jj}.sc)])
	set(gca,'YDir','normal')
	caxis(clim_sum2)
	axis square
	xlabel('Deg')
	ylabel('Deg')
	if handcut == 1
	  hold on
	  plot([-5 5],[-5 5],'g','linewidth',7)
	  plot([-5 5],[5 -5],'g','linewidth',7)
	end
	if poscut == 1
	  hold on
	  plot([-5 5],[-5 5],'r','linewidth',7)
	  plot([-5 5],[5 -5],'r','linewidth',7)
	end
	if bmcut == 1
	  hold on
	  plot([-5 5],[-5 5],'y','linewidth',7)
	  plot([-5 5],[5 -5],'y','linewidth',7)
	end
      end
    end
    
    img_str=[subdir '/full/det_row_' ro '_col_' co '_tile_' ti '_B'];
    print('-dpng', img_str);
    
    %{
    %dif    
    clf
    %  subplot(3,6,[7 14])
    %  imagesc(x_bin,y_bin,log10(abs(B)))
    %  rectangle('Position',[-4 -4 8 8],'LineWidth',2);
    %  title(['Det row ' ro ', col ' co ', tile ' ti ',  GCP ' gcpB ...
    %	', pol A, log10 scale'])
    %  set(gca,'YDir','normal')
    %  caxis(clim_sum2)
    %  colorbar
    %  axis image
    %  xlabel('x'', Deg')
    %  ylabel('y'', Deg')
    
    subplot(2,3,1)
    imagesc(x_bin,y_bin,A-B)
    axis([-4 4 -4 4])
    title(['gcp' gcpA '-gcp' gcpB])
    set(gca,'YDir','normal')
    caxis(clim_dif)
    colorbar
    axis square
    xlabel('x'', Deg')
    ylabel('y'', Deg')
    
    subplot(2,3,4)
    imagesc(x_bin,y_bin,sqrt(Astd^2+Bstd^2))
    axis([-4 4 -4 4])
    title(['gcp' gcpB '-gcp' gcpB ' rms'])
    set(gca,'YDir','normal')
    caxis([0 0.05])
    colorbar
    axis square
    xlabel('x'', Deg')
    ylabel('y'', Deg')
    
    plotpos = [3,4,5,6,9,10,11,12,15,16,17,18];
    for jj = 1:length(qq)
      if jj < 24
	handcut = 0;
	if cutmaps
	  if mapstocut(jj,ind.l(ii)) == 1
	    handcut = 1;
	  end
	end
	subplot(3,6,plotpos(jj))
	amap=qq{jj}.map(:,:,ind.b(ii));
	imagesc(x_bin,y_bin,log10(abs(amap)))
	axis([-4 4 -4 4])
	title(['gcp' gcpA '- gcp' gcpB ', dk' num2str(qq{jj}.dk) ', sc' num2str(qq{jj}.sc)])
	set(gca,'YDir','normal')
	caxis(clim_dif)
	axis square
	xlabel('x'', Deg')
	ylabel('y'', Deg')
	if handcut == 1
	  hold on
	  plot([-5 5],[-5 5],'r','linewidth',7)
	  plot([-5 5],[5 -5],'r','linewidth',7)
	end
	if poscut == 1
	  hold on
	  plot([-5 5],[-5 5],'r','linewidth',7)
	  plot([-5 5],[5 -5],'r','linewidth',7)
	end
      end
    end
    
    img_str=[subdir '/full/det_row_' ro '_col_' co '_tile_' ti '_dif'];
    print('-dpng', img_str);
    %}
    
    
    
  end
end
%}

