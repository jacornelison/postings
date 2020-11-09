function ffbm_plot_compdetail(compopt)
% function ffbm_plot_compdetail(compopt)
%
% Plot composite beam maps such that we can analyze residuals with
% varying levels of deprojection across the FPU.  Designed for a new
% pager which allows the thumbnails to dynamicall change too.
%
% Make different thumbnail plots for each of the options, plus a large
% plot with tons of information and beam functions
%
% Input composite file should have 'ad' and 'map' struct, where 'map' is
% a struct of length n_det and has (at least) a field 'T' 
% 
% INPUTS (should all be sent in with compopt)
%
%   expt:           'bicep2','keck'
%   year:           b2: unncecessary
%                   keck: 2015
% 
% OPTIONAL INPUTS
%
%   makesmall:      make thumbnail plot for all dp options
%   makefull:       make full plot 
%   compositedir:   directory from which to load composite
%                   (default beammaps/composites)
%   compositefile:  name of file in compositedir to load
%                   (default ffbm_year_all)
%   coord:          Coordinate system for rotated/composite maps:
%                   'dk0' for inputs to reduc_makesim
%                   'xpyp' for normal visualization (x'/y')
%   mapaxis:        axis limits (defaults to ad struct)
%   maxscale:       peak raw value by which to scale colors
%                   default 0.03 (integral normalizing usually gets
%                   around 0.029 ish)
%   dblim:          dB down from maxscale to plot, default -40
%   subtightscale:  spacing for 'subtightplot' subfunction, default 0.06
%   onlychflag:     load channel flags and mark those removed
%   bmflag:         load beam map flags and mark those removed
%   bmflagtype:     regular flags 'peryear' or x-year 'xyear'
%   plotdir:        directory in which to save plots

% Parse compopt
expt = compopt.expt;
year = compopt.year;
if ~isfield(compopt,'makesmall')
  compopt.makesmall = 1;
  makesmall = compopt.makesmall;
else
  makesmall = compopt.makesmall;
end
if ~isfield(compopt,'makefull')
  compopt.makefull = 1;
  makefull = compopt.makefull;
else
  makefull = compopt.makefull;
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
if ~isfield(compopt,'coord')
  compopt.coord = 'xpyp';
  coord = compopt.coord;
else
  coord = compopt.coord;
end
if ~isfield(compopt,'mapaxis')
  compopt.mapaxis = [-2 2 -2 2];
  mapaxis = compopt.mapaxis;
else
  mapaxis = compopt.mapaxis;
end
if ~isfield(compopt,'maxscale')
  compopt.maxscale = 0.03;
  maxscale = compopt.maxscale;
else
  maxscale = compopt.maxscale;
end
if ~isfield(compopt,'dblim')
  compopt.dblim = -40;
  dblim = compopt.dblim;
else
  dblim = compopt.dblim;
end
if ~isfield(compopt,'subtightscale')
  compopt.subtightscale = 0.06;
  subtightscale = compopt.subtightscale;
else
  subtightscale = compopt.subtightscale;
end
if ~isfield(compopt,'onlychflag')
  compopt.onlychflag = 1;
  onlychflag = compopt.onlychflag;
else
  onlychflag = compopt.onlychflag;
end
if ~isfield(compopt,'bmflag')
  compopt.bmflag = 0;
  bmflag = compopt.bmflag;
else
  bmflag = compopt.bmflag;
end
% Choose bmflagtype
if bmflag
  if ~isfield(compopt,'bmflagtype')
    compopt.bmflagtype = 'peryear';
    bmflagtype = 'peryear';
  else
    bmflagtype = compopt.bmflagtype;
  end
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
if bmflag % Do these independently of chflags
  switch expt
    case 'bicep2'
      chflags_bm.filebase{1} = 'fp_data/fp_data_tocut';
      chflags_bm.par_name{1} = 'pair41';
      chflags_bm.low(1) = -1;
      chflags_bm.high(1) = 0.5;
    case 'keck'
      switch bmflagtype
        case 'peryear'
          chflags_bm.filebase{1} = 'fp_data/fp_data_tocut';
          chflags_bm.par_name{1} = 'nobeams';
          chflags_bm.low(1) = -1;
          chflags_bm.high(1) = 0.5;
          chflags_bm.filebase{2} = 'fp_data/fp_data_tocut';
          chflags_bm.par_name{2} = 'addcuts';
          chflags_bm.low(2) = -1;
          chflags_bm.high(2) = 0.5;
        case 'xyear'
          chflags_bm.filebase{1} = 'fp_data/fp_data_tocut';
          chflags_bm.par_name{1} = 'nobeams_xyear';
          chflags_bm.low(1) = -1;
          chflags_bm.high(1) = 0.5;
          chflags_bm.filebase{2} = 'fp_data/fp_data_tocut';
          chflags_bm.par_name{2} = 'addcuts_xyear';
          chflags_bm.low(2) = -1;
          chflags_bm.high(2) = 0.5;
      end
  end
else
  chflags_bm = [];
end

[p ind] = get_array_info([num2str(year) '0201'],...
                         [],[],[],[],chflags);
[p_bm ind_bm] = get_array_info([num2str(year) '0201'],...
                               [],[],[],[],chflags_bm);
[p0 ind0] = get_array_info([num2str(year) '0201']);

% Make directories
switch expt
  case {'bicep2','bicep3'}
    mkdir(plotdir);
    mkdir([plotdir '/dpA/small']);
    mkdir([plotdir '/dpA/full']);
    mkdir([plotdir '/dpB/small']);
    mkdir([plotdir '/dpB/full']);
    mkdir([plotdir '/dp0000/small']);
    mkdir([plotdir '/dp0000/full']);
    mkdir([plotdir '/dp1100/small']);
    mkdir([plotdir '/dp1100/full']);
    mkdir([plotdir '/dp1101/small']);
    mkdir([plotdir '/dp1101/full']);
    mkdir([plotdir '/dp1111/small']);
    mkdir([plotdir '/dp1111/full']);
  case 'keck'
    for ii = 0:4
      mkdir([plotdir '/rx' num2str(ii)]);
      mkdir([plotdir '/rx' num2str(ii) '/dpA/small']);
      mkdir([plotdir '/rx' num2str(ii) '/dpA/full']);
      mkdir([plotdir '/rx' num2str(ii) '/dpB/small']);
      mkdir([plotdir '/rx' num2str(ii) '/dpB/full']);
      mkdir([plotdir '/rx' num2str(ii) '/dp0000/small']);
      mkdir([plotdir '/rx' num2str(ii) '/dp0000/full']);
      mkdir([plotdir '/rx' num2str(ii) '/dp1100/small']);
      mkdir([plotdir '/rx' num2str(ii) '/dp1100/full']);
      mkdir([plotdir '/rx' num2str(ii) '/dp1101/small']);
      mkdir([plotdir '/rx' num2str(ii) '/dp1101/full']);
      mkdir([plotdir '/rx' num2str(ii) '/dp1111/small']);
      mkdir([plotdir '/rx' num2str(ii) '/dp1111/full']);
    end
end

% Get yo composite (I guess this also works with regular maps)
load([compositedir '/' compositefile]);

for ii = 1:length(ind.la)
 
  if mod(ii,10) == 0
    disp(['Plotting pair ' int2str(ii) ' / ' int2str(length(ind.la))]);
  end
  
  ro = int2str(p.det_row(ind.la(ii)));
  co = int2str(p.det_col(ind.la(ii)));
  ti = int2str(p.tile(ind.la(ii)));
  rx = int2str(p.rx(ind.la(ii)));
  mce = int2str(p.mce(ind.la(ii)));
  gcp_a = int2str(p.gcp(ind.la(ii)));
  gcp_b = int2str(p.gcp(ind.lb(ii)));

  map_a = map(ind.la(ii)).T;
  map_b = map(ind.lb(ii)).T;
  % Diff beam
  dp0000 = map_a - map_b;
  % Fit sigma
  sig = nanmedian([p.fwhm_maj(ind.la(ii)),...
                   p.fwhm_min(ind.la(ii)),...
                   p.fwhm_maj(ind.lb(ii)),...
                   p.fwhm_min(ind.lb(ii))])...
        / (2*sqrt(2*log(2)));
  % Deproject diff map
  [c T] = ffbm_deprojmap(inpaint_nans(dp0000),ad,sig);
  dp1100 = dp0000 - c(1).*T{1} - c(2).*T{2} - c(3).*T{3};
  dp1101 = dp1100 - c(5).*T{5} - c(6).*T{6};
  dp1111 = dp1101 - c(4).*T{4};
  
  % Plotting directories - keep by rx for Keck
  switch expt
    case {'bicep2','bicep3'}
      subplotdir = plotdir;
    case 'keck'
      subplotdir = [plotdir '/rx' num2str(p.rx(ind.la(ii)))];
  end

  ax1 = ad.t_val_deg{1}; % See documentation in
  ax2 = ad.t_val_deg{2}; % ffbm_rotatemaps
  
  if makefull
     
    fh = figure(1);
    set(fh,'Position',[0 0 1800 600],'Visible','off')
    clf
    colormap('jet')
        
    % A detector  
    subtightplot(2,6,1,subtightscale)
    imagescnan(ax2,ax1,10*log10(abs(map_a)));
    axis(mapaxis)
    title(['A: GCP ' gcp_a])
    set(gca,'YDir','normal')
    colorbar
    caxis([10*log10(maxscale)+dblim 10*log10(maxscale)])
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
    
    % A beam function
    [l,B_l] = b2bl(ad,map_a,1);
    subtightplot(2,6,7,subtightscale)
    plot(l,B_l);
    xlim([0 1000]);
    ylim([0 3.5e-6]);
    text(200,3.25e-6,...
         ['Map sum: ' num2str(nansum(map_a(:)),'%.3s')]);
    
    % B detector  
    subtightplot(2,6,2,subtightscale)
    imagescnan(ax2,ax1,10*log10(abs(map_b)));
    axis(mapaxis)
    title(['B: GCP ' gcp_b])
    set(gca,'YDir','normal')
    colorbar
    caxis([10*log10(maxscale)+dblim 10*log10(maxscale)])
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
    
    % B beam function
    [l,B_l] = b2bl(ad,map_b,1);
    subtightplot(2,6,8,subtightscale)
    plot(l,B_l);
    xlim([0 1000]);
    ylim([0 3.5e-6]);
    text(200,3.25e-6,...
         ['Map sum: ' num2str(nansum(map_b(:)),'%.3s')]);
    
    % Standard difference
    subtightplot(2,6,3,subtightscale)
    imagescnan(ax2,ax1,dp0000);
    axis(mapaxis)
    title('0000: Diff beam')
    set(gca,'YDir','normal')
    caxis([-maxscale*0.1 maxscale*0.1])
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
    
    % Beam function
    [l,B_l] = b2bl(ad,dp0000,1);
    subtightplot(2,6,9,subtightscale)
    plot(l,B_l);
    xlim([0 1000]);
    ylim([0 1e-7]);
    text(200,0.9e-7,...
         ['Map sum: ' num2str(nansum(dp0000(:)),'%.3s')]);
    text(200,0.8e-7,...
         ['Map std: ' num2str(nanstd(dp0000(:)),'%.3s')]);
    
    % Deprojecting relgain/diff pointing
    subtightplot(2,6,4,subtightscale)
    imagescnan(ax2,ax1,dp1100);
    axis(mapaxis)
    title('1100: Relgain/diff pointing')
    set(gca,'YDir','normal')
    caxis([-maxscale*0.025 maxscale*0.025])
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
    
    % Beam function
    [l,B_l] = b2bl(ad,dp1100,1);
    subtightplot(2,6,10,subtightscale)
    plot(l,B_l);
    xlim([0 1000]);
    ylim([0 0.5e-7]);
    text(200,0.45e-7,...
         ['Map sum: ' num2str(nansum(dp1100(:)),'%.3s')]);
    text(200,0.4e-7,...
         ['Map std: ' num2str(nanstd(dp1100(:)),'%.3s')]);
    
    % Deprojecting relgain/diff pointing/diffellip
    subtightplot(2,6,5,subtightscale)
    imagescnan(ax2,ax1,dp1101);
    axis(mapaxis)
    title('1101: + diff ellip')
    set(gca,'YDir','normal')
    caxis([-maxscale*0.025 maxscale*0.025])
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
          
    % Beam function
    [l,B_l] = b2bl(ad,dp1101,1);
    subtightplot(2,6,11,subtightscale)
    plot(l,B_l);
    xlim([0 1000]);
    ylim([0 0.5e-7]);
    text(200,0.45e-7,...
         ['Map sum: ' num2str(nansum(dp1101(:)),'%.3s')]);
    text(200,0.4e-7,...
         ['Map std: ' num2str(nanstd(dp1101(:)),'%.3s')]);
    
    % Deprojecting all
    subtightplot(2,6,6,subtightscale)
    imagescnan(ax2,ax1,dp1111);
    axis(mapaxis)
    title('1111: + diff bw')
    set(gca,'YDir','normal')
    caxis([-maxscale*0.025 maxscale*0.025])
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
    
    % Beam function
    [l,B_l] = b2bl(ad,dp1111,1);
    subtightplot(2,6,12,subtightscale)
    plot(l,B_l);
    xlim([0 1000]);
    ylim([0 0.5e-7]);
    text(200,0.45e-7,...
         ['Map sum: ' num2str(nansum(dp1111(:)),'%.3s')]);
    text(200,0.4e-7,...
         ['Map std: ' num2str(nanstd(dp1111(:)),'%.3s')]);
    
    gtitle(['Tile ' ti ' Row ' ro ' Col ' co ]);
    
    % Save plot
    img_str = [subplotdir '/dpA/full/det_row_' ro '_col_' co '_tile_' ti];
    mkpng(img_str);
    img_str = [subplotdir '/dpB/full/det_row_' ro '_col_' co '_tile_' ti];
    mkpng(img_str);
    img_str = [subplotdir '/dp0000/full/det_row_' ro '_col_' co '_tile_' ti];
    mkpng(img_str);
    img_str = [subplotdir '/dp1100/full/det_row_' ro '_col_' co '_tile_' ti];
    mkpng(img_str);
    img_str = [subplotdir '/dp1101/full/det_row_' ro '_col_' co '_tile_' ti];
    mkpng(img_str);
    img_str = [subplotdir '/dp1111/full/det_row_' ro '_col_' co '_tile_' ti];
    mkpng(img_str);
                   
  end % makefull
    
  if makesmall
    
    % Plot Xs through maps not going into CMB analysis
    % ind0 is before channel cuts, ind is after
    if onlychflag
        
      notrgl = 0;
      cutbych = 0;
      
      % Case: not on RGL list before or after ch cuts
      if ~ismember(ind.la(ii),ind0.rgl) & ~ismember(ind.la(ii),ind.rgl)
        notrgl = 1;
      end
      % Case: on RGL list before ch cuts but not after
      if ismember(ind.la(ii),ind0.rgl) & ~ismember(ind.la(ii),ind.rgl)
        cutbych = 1;
      end
      
    else
            
      notrgl = 0;
      cutbych = 0;

    end % chflag shit
      
    % Plot +s through maps cut by beam map channel flags
    % Too annoying to find bad beam maps which are NOT rgls...oh well
    if bmflag
       
      badbm = 0;
      
      % Case: on RGL list before bm flags but not after
      if ismember(ind_bm.la(ii),ind0.rgl) & ...
              ~ismember(ind_bm.la(ii),ind_bm.rgl)
        badbm = 1;
      end
      
    else
        
      badbm = 0;
      
    end
    
    % A detector
    fh = figure(1);
    set(fh,'Position',[10 10 670 500],'Visible','off')
    clf
    colormap('jet')
    imagescnan(ax2,ax1,10*log10(abs(map_a)));
    caxis([10*log10(maxscale)+dblim 10*log10(maxscale)])

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
    if badbm
      hold on
      plot([0 0],[-2 2],'k','linewidth',20)
      plot([-2 2],[0 0],'k','linewidth',20)
    end

    axis tight
    axis equal
    axis(mapaxis)
    set(gca,'YDir','normal')
    if strcmp(coord,'xpyp')
        set(gca,'XDir','reverse')
    end
    set(gca,'xtick',[],'ytick',[])

    % Save image
    img_str = [subplotdir '/dpA/small/det_row_' ro '_col_' co '_tile_' ti];
    mkpng(img_str);
    % 400x400 window starting at (145,40)
    system(['mogrify -crop ''400x400+145+40'' ' subplotdir ...
            '/dpA/small/det_row_' ro '_col_' co ...
            '_tile_' ti '*']);
    % Resize to 75 pixels (400 * 0.18) 
    system(['mogrify -resize ''18%x18%'' ' subplotdir ...
            '/dpA/small/det_row_' ro '_col_' co ...
            '_tile_' ti '*']);
    
    % B detector
    fh = figure(1);
    set(fh,'Position',[10 10 670 500],'Visible','off')
    clf
    colormap('jet')
    imagescnan(ax2,ax1,10*log10(abs(map_b)));
    caxis([10*log10(maxscale)+dblim 10*log10(maxscale)])

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
    if badbm
      hold on
      plot([0 0],[-2 2],'k','linewidth',20)
      plot([-2 2],[0 0],'k','linewidth',20)
    end
    
    axis tight
    axis equal
    axis(mapaxis)
    set(gca,'YDir','normal')
    if strcmp(coord,'xpyp')
        set(gca,'XDir','reverse')
    end
    set(gca,'xtick',[],'ytick',[])

    % Save image
    img_str = [subplotdir '/dpB/small/det_row_' ro '_col_' co '_tile_' ti];
    mkpng(img_str);
    % 400x400 window starting at (145,40)
    system(['mogrify -crop ''400x400+145+40'' ' subplotdir ...
            '/dpB/small/det_row_' ro '_col_' co ...
            '_tile_' ti '*']);
    % Resize to 75 pixels (400 * 0.18) 
    system(['mogrify -resize ''18%x18%'' ' subplotdir ...
            '/dpB/small/det_row_' ro '_col_' co ...
            '_tile_' ti '*']);
    
    % dp0000
    fh = figure(1);
    set(fh,'Position',[10 10 670 500],'Visible','off')
    clf
    colormap('jet')
    imagescnan(ax2,ax1,dp0000);
    caxis([-maxscale*0.1 maxscale*0.1])

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
    if badbm
      hold on
      plot([0 0],[-2 2],'k','linewidth',20)
      plot([-2 2],[0 0],'k','linewidth',20)
    end
    
    axis tight
    axis equal
    axis(mapaxis)
    set(gca,'YDir','normal')
    if strcmp(coord,'xpyp')
        set(gca,'XDir','reverse')
    end
    set(gca,'xtick',[],'ytick',[])

    % Save image
    img_str = [subplotdir '/dp0000/small/det_row_' ro '_col_' co '_tile_' ti];
    mkpng(img_str);
    % 400x400 window starting at (145,40)
    system(['mogrify -crop ''400x400+145+40'' ' subplotdir ...
            '/dp0000/small/det_row_' ro '_col_' co ...
            '_tile_' ti '*']);
    % Resize to 75 pixels (400 * 0.18) 
    system(['mogrify -resize ''18%x18%'' ' subplotdir ...
            '/dp0000/small/det_row_' ro '_col_' co ...
            '_tile_' ti '*']);
    
    % dp1100
    fh = figure(1);
    set(fh,'Position',[10 10 670 500],'Visible','off')
    clf
    colormap('jet')
    imagescnan(ax2,ax1,dp1100);
    caxis([-maxscale*0.025 maxscale*0.025])

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
    if badbm
      hold on
      plot([0 0],[-2 2],'k','linewidth',20)
      plot([-2 2],[0 0],'k','linewidth',20)
    end
    
    axis tight
    axis equal
    axis(mapaxis)
    set(gca,'YDir','normal')
    if strcmp(coord,'xpyp')
        set(gca,'XDir','reverse')
    end
    set(gca,'xtick',[],'ytick',[])

    % Save image
    img_str = [subplotdir '/dp1100/small/det_row_' ro '_col_' co '_tile_' ti];
    mkpng(img_str);
    % 400x400 window starting at (145,40)
    system(['mogrify -crop ''400x400+145+40'' ' subplotdir ...
            '/dp1100/small/det_row_' ro '_col_' co ...
            '_tile_' ti '*']);
    % Resize to 75 pixels (400 * 0.18) 
    system(['mogrify -resize ''18%x18%'' ' subplotdir ...
            '/dp1100/small/det_row_' ro '_col_' co ...
            '_tile_' ti '*']);
    
    % dp1101
    fh = figure(1);
    set(fh,'Position',[10 10 670 500],'Visible','off')
    clf
    colormap('jet')
    imagescnan(ax2,ax1,dp1101);
    caxis([-maxscale*0.025 maxscale*0.025])

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
    if badbm
      hold on
      plot([0 0],[-2 2],'k','linewidth',20)
      plot([-2 2],[0 0],'k','linewidth',20)
    end
    
    axis tight
    axis equal
    axis(mapaxis)
    set(gca,'YDir','normal')
    if strcmp(coord,'xpyp')
        set(gca,'XDir','reverse')
    end
    set(gca,'xtick',[],'ytick',[])

    % Save image
    img_str = [subplotdir '/dp1101/small/det_row_' ro '_col_' co '_tile_' ti];
    mkpng(img_str);
    % 400x400 window starting at (145,40)
    system(['mogrify -crop ''400x400+145+40'' ' subplotdir ...
            '/dp1101/small/det_row_' ro '_col_' co ...
            '_tile_' ti '*']);
    % Resize to 75 pixels (400 * 0.18) 
    system(['mogrify -resize ''18%x18%'' ' subplotdir ...
            '/dp1101/small/det_row_' ro '_col_' co ...
            '_tile_' ti '*']);
    
    % dp1111
    fh = figure(1);
    set(fh,'Position',[10 10 670 500],'Visible','off')
    clf
    colormap('jet')
    imagescnan(ax2,ax1,dp1111);
    caxis([-maxscale*0.025 maxscale*0.025])

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
    if badbm
      hold on
      plot([0 0],[-2 2],'k','linewidth',20)
      plot([-2 2],[0 0],'k','linewidth',20)
    end
    
    axis tight
    axis equal
    axis(mapaxis)
    set(gca,'YDir','normal')
    if strcmp(coord,'xpyp')
        set(gca,'XDir','reverse')
    end
    set(gca,'xtick',[],'ytick',[])

    % Save image
    img_str = [subplotdir '/dp1111/small/det_row_' ro '_col_' co '_tile_' ti];
    mkpng(img_str);
    % 400x400 window starting at (145,40)
    system(['mogrify -crop ''400x400+145+40'' ' subplotdir ...
            '/dp1111/small/det_row_' ro '_col_' co ...
            '_tile_' ti '*']);
    % Resize to 75 pixels (400 * 0.18) 
    system(['mogrify -resize ''18%x18%'' ' subplotdir ...
            '/dp1111/small/det_row_' ro '_col_' co ...
            '_tile_' ti '*']);
    

  end % makesmall

end

return