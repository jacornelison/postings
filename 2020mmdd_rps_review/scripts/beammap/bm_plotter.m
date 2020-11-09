function bm_plotter(bmopt)
% bm_plotter(bmopt)
%
% Function to plot standard far field beam maps.
% Makes small thumbnails (from windowed map) and large plots with fits (from
% full and windowed maps)
%
% INPUTS (should all be passed in with bmopt)
%
%   expt:      'b2','keck','bicep3'
%   rxNum:     For Keck: rx number, array of rx, or 'all' (default)
%   mapcoord:  'azel_ap' (default),'raw' (used to get map filename)
%   demodtype:   'square' (default), 'choplet'
%   fullmapdir:directory from which to load full maps, beammaps/maps
%   winmapdir: directory from which to load win maps, 
%              beammaps/maps_windowed default
%   suffix:    optional string to add to map name, '' default
%   filename:  string with timestamp of map file (i.e. '20150101_000000' or
%              p = get_bm_info(year); filename = p.filename{bmnum};)
%   plotdir:   directory in which to save images.  Will create if doesn't
%              exist, default plots/filename
%
% OPTIONAL
%
%   clim_sum:  color limits for normal maps, default is expt-dependent
%   clim_dif:  color limits for difference maps, default [-0.1 0.1]
%   dblim:     color limits for peak-normalized dB plots, default [-30 0]
%   subtightscale: spacing for 'subtightplot' subfunction, default 0.04
%   dosmall:   1 (default) to make small maps, 0 to skip
%   dofull:    1 (default) to make full maps, 0 to skip
%   fullaxis:  limits for full map, i.e. [-20 20 -13 23].  Defaults is
%              max/min of the fullmap, which may extend annoyingly far. 

% Parse bmopt and set defaults
expt = bmopt.expt;
if ~isfield(bmopt,'rxNum')
  switch expt
    case 'keck'
      rxNum = 0:4;
  end
else
  if strcmp(bmopt.rxNum,'all')
    rxNum = 0:4;
  else
    rxNum = bmopt.rxNum;
  end
end
if ~isfield(bmopt,'mapcoord')
  mapcoord = 'azel_ap';
else
  mapcoord = bmopt.mapcoord;
end
if ~isfield(bmopt,'demodtype')
  demodtype = 'square';
else
  demodtype = bmopt.demodtype;
end
if ~isfield(bmopt,'fullmapdir')
  fullmapdir = 'beammaps/maps';
else
  fullmapdir = bmopt.fullmapdir;
end
if ~isfield(bmopt,'winmapdir')
  winmapdir = 'beammaps/maps_windowed';
else
  winmapdir = bmopt.winmapdir;
end
if ~isfield(bmopt,'suffix')
  suffix = '';
else
  suffix = bmopt.suffix;
end
filename = bmopt.filename;
if ~isfield(bmopt,'plotdir')
  plotdir = ['plots/' filename];
else
  plotdir = bmopt.plotdir;
end
if ~isfield(bmopt,'dblim')
  dblim = [-30 0];
else
  dblim = bmopt.dblim;
end
if ~isfield(bmopt,'subtightscale')
  subtightscale = 0.04;
else
  subtightscale = bmopt.subtightscale;
end
if ~isfield(bmopt,'dosmall')
  dosmall = 1;
else
  dosmall = bmopt.dosmall;
end
if ~isfield(bmopt,'dolarge')
  dolarge = 1;
else
  dolarge = bmopt.dolarge;
end
% If empty, call subfunction in plot loop to figure out clim_sum
% We save to a different name here because it can change in the loop
if ~isfield(bmopt,'clim_sum')
  clim_sum_in = [];
else
  clim_sum_in = bmopt.clim_sum;
end
if ~isfield(bmopt,'clim_dif')
  clim_dif = [-0.1 0.1];
else
  clim_dif = bmopt.clim_dif;
end
if ~isfield(bmopt,'fullaxis')
  fullaxis = [];
else
  fullaxis = bmopt.fullaxis;
end

% Load up saved maps
fmapname = [fullmapdir,'/map_',filename,'_',demodtype,...
      '_',mapcoord,suffix];
f = load(fmapname);
f = f.bm;
wmapname = [winmapdir,'/mapwin_',filename,'_',demodtype,...
      '_',mapcoord,suffix,'_2deg']; % hard code 2deg win for now
w = load(wmapname);
w = w.bm;

% Make directories
switch expt
  case {'b2','bicep3'}
    mkdir(plotdir);
    mkdir([plotdir '/small']);
    mkdir([plotdir '/full']);
  case 'keck'
    for ii = 1:length(rxNum)
      mkdir([plotdir '_rx' num2str(rxNum(ii))]);
      mkdir([plotdir '_rx' num2str(rxNum(ii)) '/small']);
      mkdir([plotdir '_rx' num2str(rxNum(ii)) '/full']);
    end
end

[pp ind] = get_array_info(filename(1:8));

% Find the right rotation angle to put thumbnail plots in x'/y'
% See BKCMB posting 20150714_bmaxes for details
switch expt
  case {'b2'}
    %f.dk = f.dk + 90 + 180;
    %f.dk = double(f.dk);
    dk = f.dk + pp.drumangle;
  case {'keck','bicep3'}
    cutind = ismember(pp.rx,rxNum);
    pp = rmfield(pp,'expt'); % for structcut to work
    pp = structcut(pp,cutind);
    ind = make_ind(pp);
    % Cut down more stuff for each rx
    f.map = f.map(:,:,cutind);
    w.map = w.map(:,:,cutind);
    w.A = w.A(:,cutind);
    %f.dk = 90 + (f.dk + pp.drumangle(1));  % CHANGE for individual dets
    %dk = 90 + (f.dk + pp.drumangle);
    rot = -(90 + f.dk + pp.drumangle);
end 

% Gaussian fit from bm_makemap
F = w.A;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Small plots (thumbnails)
if dosmall

  fh = figure(1);
  set(fh,'Position',[10 10 670 500],'Visible','off') % For cropping to work 
  clf
  colormap('jet')
  
  for ii = 1:length(ind.la)
    if mod(ii,10) == 0
      disp([filename sprintf(': Small map, rx%i, %i/%i',pp.rx(ind.la(ii)),...
	  ii,length(ind.la))]);
    end
    
    ro = int2str(pp.det_row(ind.la(ii)));
    co = int2str(pp.det_col(ind.la(ii)));
    ti = int2str(pp.tile(ind.la(ii)));
    rx = int2str(pp.rx(ind.la(ii)));

    %A = imrotate(w.map(:,:,ind.la(ii)),dk(ind.la(ii)),'bilinear','crop');
    %B = imrotate(w.map(:,:,ind.lb(ii)),dk(ind.la(ii)),'bilinear','crop');
    
    WA = w.map(:,:,ind.la(ii));
    WB = w.map(:,:,ind.lb(ii));
    
    % Flip and rotate to x'/y' (see BM axes posting)
    WA_xpyp = flipud(WA);
    WB_xpyp = flipud(WB);
    WA_xpyp = imrotate(WA_xpyp,rot(ind.la(ii)),'bilinear','crop');
    WB_xpyp = imrotate(WB_xpyp,rot(ind.la(ii)),'bilinear','crop');
    % These maps have dimension 1 = increasing x'
    %                 dimension 2 = decreasing y'

    % Now change axes to match above map
    yp = fliplr(w.x_bin); % Decreasing y prime
    xp = w.y_bin; % increasing x prime
    
    switch expt
      case {'b2','bicep3'}
	subplotdir = plotdir;
      case {'keck'}
	subplotdir = [plotdir '_rx' num2str(pp.rx(ind.la(ii)))];
    end
    
    % Get clim_sum
    if isempty(clim_sum_in)
      clim_sum = get_clim_sum(expt,pp,ind,ii,filename);
    else
      clim_sum = clim_sum_in;
    end
    
    % Plot and save A pol
    imagescnan(yp,xp,WA_xpyp);
    prettify_small(clim_sum,'xpyp');
    img_str = [subplotdir '/small/det_row_' ro '_col_' co ...
	  '_tile_' ti '_A_sm'];
    print('-dpng', img_str);
    
    % Plot and save B pol
    imagescnan(yp,xp,WB_xpyp);
    prettify_small(clim_sum,'xpyp');
    img_str = [subplotdir '/small/det_row_' ro '_col_' co ...
	  '_tile_' ti '_B_sm'];
    print('-dpng', img_str);
    
    % Plot and save A/|A|-B/|B|
    imagescnan(yp,xp,...
               WA_xpyp./w.A(1,ind.la(ii)) - WB_xpyp./w.A(1,ind.lb(ii)));
    prettify_small(clim_dif,'xpyp');
    img_str = [subplotdir '/small/det_row_' ro '_col_' co ...
	  '_tile_' ti '_dif_sm'];
    print('-dpng', img_str);
    
    % Pause for a split second to let the files appear
    pause(0.5)
    
    % Crop
    system(['mogrify -crop ''666x666+244+69'' ' subplotdir ...
            '/small/det_row_' ro '_col_' co '_tile_' ti '_*']);
    
    % Shrink
    system(['mogrify -resize ''7.5%x7.5%'' ' subplotdir ...
            '/small/det_row_' ro '_col_' co '_tile_' ti '_*']);
    
  end % over detectors
end % dosmall

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if dolarge
  fh = figure(1);
  clf

  % Get coordinate axes
  az_ap_full = f.x_bin;
  el_ap_full = f.y_bin;
  
  az_ap_win = w.x_bin - mean(w.x_bin);
  el_ap_win = w.y_bin - mean(w.y_bin);

  yp = fliplr(w.x_bin - mean(w.x_bin)); % Decreasing y prime
  xp = w.y_bin - mean(w.y_bin); % increasing x prime

  % Create full plots
  set(fh,'Position',[0 0 1600 800],'visible','off')

  for ii = 1:length(ind.la)

    clf
    
    det_a = ind.la(ii);
    det_b = ind.lb(ii);

    if mod(ii,10) == 0
      disp([filename sprintf(': Full map, rx%i, %i/%i',pp.rx(det_a),...
	  ii,length(ind.la))]);
    end

    ro = int2str(pp.det_row(det_a));
    co = int2str(pp.det_col(det_a));
    ti = int2str(pp.tile(det_a));
    rx = int2str(pp.rx(det_a));
    gcpA = int2str(pp.gcp(det_a));
    gcpB = int2str(pp.gcp(det_b));
    mce = int2str(pp.mce(det_a));
    
    A = f.map(:,:,det_a);
    B = f.map(:,:,det_b);
    ZA = gauss2d(F(:,det_a),w.x_bin,w.y_bin,'std');
    ZB = gauss2d(F(:,det_b),w.x_bin,w.y_bin,'std');
    WA = w.map(:,:,det_a);
    WB = w.map(:,:,det_b);
    
    % Flip and rotate to x'/y' (see BM axes posting)
    WA_xpyp = flipud(WA);
    WB_xpyp = flipud(WB);
    WA_xpyp = imrotate(WA_xpyp,rot(det_a),'bilinear','crop');
    WB_xpyp = imrotate(WB_xpyp,rot(det_b),'bilinear','crop');
    % These maps have dimension 1 = increasing x'
    %                 dimension 2 = decreasing y'
    
    % Fit parameters
    mA = w.A(1,det_a);
    mB = w.A(1,det_b);
    sig_a = sqrt(((w.A(4,det_a))^2 + (w.A(5,det_a))^2)/2);
    sig_b = sqrt(((w.A(4,det_b))^2 + (w.A(5,det_b))^2)/2);
    ell_a = ((w.A(4,det_a))^2 - (w.A(5,det_a))^2) / ...
            ((w.A(4,det_a))^2 + (w.A(5,det_a))^2);
    ell_b = ((w.A(4,det_b))^2 - (w.A(5,det_b))^2) / ...
            ((w.A(4,det_b))^2 + (w.A(5,det_b))^2);
    
    switch expt
      case {'b2','bicep3'}
	subplotdir = plotdir;
      case {'keck'}
	subplotdir = [plotdir '_rx' num2str(pp.rx(det_a))];
    end  

    % Get clim_sum
    if isempty(clim_sum_in)
      clim_sum = get_clim_sum(expt,pp,ind,ii,filename);
    else
      clim_sum = clim_sum_in;
    end
      
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Plot and save A pol
    % full map
    subtightplot(3,5,[1 13],subtightscale)
    if isnan(mA)
      imagescnan(az_ap_full,el_ap_full,10*log10(abs(A)));
      prettify_huge(...
          [10*log10(clim_sum(2)) - range(dblim),...
           10*log10(clim_sum(2))],fullaxis);
    else
      imagescnan(az_ap_full,el_ap_full,10*log10(abs(A./mA)));
      prettify_huge(dblim,fullaxis);
    end
    xlabel('Az apparent, Degrees')
    ylabel('El apparent, Degrees')
    switch expt
      case 'b2'
        title(['Tile ',ti,' row ',ro,' col ',co,' GCP ',gcpA ...
               ', pol A, Peak-normalized dB scale'])
      case 'keck'
        title(['Rx ',rx,' tile ',ti,' row ',ro,' col ',co,' GCP ',gcpA ...
               ', pol A, Peak-normalized dB scale'])
      case 'bicep3'
        title(['MCE ',mce,' tile ',ti,' row ',ro,' col ',co,' GCP ',gcpA ...
               ', pol A, Peak-normalized dB scale'])
    end


    % Windowed map 
    subtightplot(3,5,4,subtightscale)
    imagescnan(az_ap_win,el_ap_win,WA)
    title('Map, Linear, az/el ap')
    prettify_large(clim_sum,'azel_ap') 
    % Beam fit
    subtightplot(3,5,9,subtightscale)
    imagescnan(az_ap_win,el_ap_win,ZA)
    title('Fit, az/el ap')
    prettify_large(clim_sum,'azel_ap') 
    % Residual
    subtightplot(3,5,14,subtightscale)
    imagescnan(az_ap_win,el_ap_win,(WA - ZA)./mA)  
    title('Data - Fit, az/el ap')
    xlabel('Az ap')
    ylabel('El ap')
    prettify_large(clim_dif,'azel_ap') 

    % Now in x'/y'
    subtightplot(3,5,5,subtightscale)
    imagescnan(yp,xp,WA_xpyp);
    title('Map, Linear, xp/yp')
    prettify_large(clim_sum,'xpyp')
    
    subtightplot(3,5,10,subtightscale)
    imagescnan(yp,xp,10*log10(abs(WA_xpyp./mA)));
    title('Map, dBpeak, xp/yp')
    prettify_large(dblim,'xpyp')
    xlabel('Y prime')
    ylabel('X prime')
    
    % And basic statistics
    subtightplot(3,5,15,subtightscale)
    axis([-2 2 -2 2])
    text(-1.8,1.1,['Amplitude: ' sprintf('%3.0f',mA)]);
    text(-1.8,-1.1,['Beamwidth: ' sprintf('%0.3f',sig_a) ' deg']);
    %text(-1.8,-1.2,['Ellipticity: ' sprintf('%0.2f',ell_a)]);
    axis off
    axis tight
    
    img_str = [subplotdir '/full/det_row_' ro '_col_' co '_tile_' ti '_A'];
    mkpng(img_str);
 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Plot and save B pol
    % full map
    clf;
    subtightplot(3,5,[1 13],subtightscale)
    if isnan(mB)
      imagescnan(az_ap_full,el_ap_full,10*log10(abs(B)));
      prettify_huge(...
          [10*log10(clim_sum(2)) - range(dblim),...
           10*log10(clim_sum(2))],fullaxis);
    else
      imagescnan(az_ap_full,el_ap_full,10*log10(abs(B./mB)));
      prettify_huge(dblim,fullaxis);
    end
    xlabel('Az apparent, Degrees')
    ylabel('El apparent, Degrees')
    switch expt
      case 'b2'
        title(['Tile ',ti,' row ',ro,' col ',co,' GCP ',gcpB ...
               ', pol B, Peak-normalized dB scale'])
      case 'keck'
        title(['Rx ',rx,' tile ',ti,' row ',ro,' col ',co,' GCP ',gcpB ...
               ', pol B, Peak-normalized dB scale'])
      case 'bicep3'
        title(['MCE ',mce,' tile ',ti,' row ',ro,' col ',co,' GCP ',gcpB ...
               ', pol B, Peak-normalized dB scale'])
    end

    % Windowed map 
    subtightplot(3,5,4,subtightscale)
    imagescnan(az_ap_win,el_ap_win,WB)
    title('Map, Linear, az/el ap')
    prettify_large(clim_sum,'azel_ap') 
    % Beam fit
    subtightplot(3,5,9,subtightscale)
    imagescnan(az_ap_win,el_ap_win,ZB)
    title('Fit, az/el ap')
    prettify_large(clim_sum,'azel_ap') 
    % Residual
    subtightplot(3,5,14,subtightscale)
    imagescnan(az_ap_win,el_ap_win,(WB - ZB)./mB)  
    title('Data - Fit, az/el ap')
    xlabel('Az ap')
    ylabel('El ap')
    prettify_large(clim_dif,'azel_ap') 

    % Now in x'/y'
    subtightplot(3,5,5,subtightscale)
    imagescnan(yp,xp,WB_xpyp);
    title('Map, Linear, xp/yp')
    prettify_large(clim_sum,'xpyp')
    
    subtightplot(3,5,10,subtightscale)
    imagescnan(yp,xp,10*log10(abs(WB_xpyp./mB)));
    title('Map, dBpeak, xp/yp')
    prettify_large(dblim,'xpyp')
    xlabel('Y prime')
    ylabel('X prime')
    
    % And basic statistics
    subtightplot(3,5,15,subtightscale)
    axis([-2 2 -2 2])
    text(-1.8,1.1,['Amplitude: ' sprintf('%3.0f',mB)]);
    text(-1.8,-1.1,['Beamwidth: ' sprintf('%0.3f',sig_b) ' deg']);
    %text(-1.8,-1.2,['Ellipticity: ' sprintf('%0.2f',ell_b)]);
    axis off
    axis tight
    
    img_str = [subplotdir '/full/det_row_' ro '_col_' co '_tile_' ti '_B'];
    mkpng(img_str);


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Plot and save A/|A|-B/|B|  
    % Full map
    clf;
    subtightplot(3,5,[1 13],subtightscale)
    imagescnan(az_ap_full,el_ap_full,A./mA - B./mB)
    switch expt
      case 'b2'
        title(['Tile ',ti,' row ',ro,' col ',co,' GCP ',...
               gcpA,'-',gcpB,' A-B']);
      case 'keck'
        title(['Rx ',rx,' tile ',ti,' row ',ro,' col ',co,' GCP ',...
               gcpA,'-',gcpB,' A-B']);
      case 'bicep3'
        title(['MCE ',mce,' tile ',ti,' row ',ro,' col ',co,' GCP ',...
               gcpA,'-',gcpB,' A-B']);
    end
    prettify_huge(clim_dif,fullaxis);
    
    % Windowed map 
    subtightplot(3,5,4,subtightscale)
    imagescnan(az_ap_win,el_ap_win,WA./mA - WB./mB)
    title('Map, Linear, az/el ap')
    prettify_large(clim_dif,'azel_ap')
    % Beam fit
    subtightplot(3,5,9,subtightscale)
    imagescnan(az_ap_win,el_ap_win,ZA./mA - ZB./mB)
    title('Fit, az/el ap')
    prettify_large(clim_dif,'azel_ap')
    % Residual
    subtightplot(3,5,14,subtightscale)
    imagescnan(az_ap_win,el_ap_win,(WA - ZA)./mA - (WB - ZB)./mB)
    title('Data - Fit, az/el ap')
    xlabel('Az ap')
    ylabel('El ap')
    prettify_large(clim_dif./5,'azel_ap')
  
    % Now in x'/y'
    subtightplot(3,5,5,subtightscale)
    imagescnan(yp,xp,WA_xpyp./mA - WB_xpyp./mB)
    title('Map, Linear, xp/yp')
    prettify_large(clim_dif,'xpyp');
    xlabel('Y prime')
    ylabel('X prime')
    
    img_str = [subplotdir '/full/det_row_' ro '_col_' co '_tile_' ti '_dif'];
    mkpng(img_str);
    
  end
end
close all;

return % Main function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function clim_sum = get_clim_sum(expt,pp,ind,ii,filename)
% Color axes are unfortunately experiment and frequency-specific, so choose
% them here

bmtime = date2mjd(str2num(filename(1:4)),str2num(filename(5:6)),...
    str2num(filename(7:8)),str2num(filename(10:11)),...
    str2num(filename(12:13)),str2num(filename(14:15)));

switch expt
  case 'b2'
    clim_sum = [0 300];
  case 'keck'
    % Before 2016 we used the 18" chopper
    if bmtime < date2mjd(2016)
      switch pp.band(ind.la(ii))
	case {100,220}
	  clim_sum = [0 200];
	case 150
	  clim_sum = [0 300];
      end
    else
      clim_sum = [0 800];
    end
  case 'bicep3'
    % Before 2015-03-17 we had the 24" and after the 18"
    % For 2016 24" 
    if bmtime < date2mjd(2015,03,17,00,00,00)
      clim_sum = [0 200];
    else
      clim_sum = [0 150];
    end
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function prettify_small(cax,coord)

set(gca,'YDir','normal');
switch coord
  case 'xpyp'
    set(gca,'XDir','reverse');
end
caxis(cax)
axis equal
set(gca,'xtick',[],'ytick',[])

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function prettify_huge(cax,fullaxis)

set(gca,'YDir','normal')
caxis(cax)
colorbar('location','eastoutside')
axis equal
if ~isempty(fullaxis)
  axis(fullaxis)
else
  axis tight
end


return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function prettify_large(cax,coord)

set(gca,'YDir','normal');
switch coord
  case 'xpyp'
    set(gca,'XDir','reverse');
end
caxis(cax)
colorbar('location','eastoutside')
axis equal
axis tight

return
