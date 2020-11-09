function ffbm_plot_beamparam_tile(compopt)
% ffbm_plot_beamparam_tile(compopt)
%
% Load up beam parameters and plot them in FPU format.  This is closer to
% a publication-quality plot. 
%
% Note that this uses just the final saved csv file, whereas
% ffbm_plot_beamparam uses the full beamparam file
%
% INPUTS (should all be sent in with compopt)
%
%  expt:               'bicep2','keck','bicep3'
%  year:               bicep2: unnecessary
%                      keck: 2012-2016
%                      bicep3: 2015, 2016
%
% OPTIONAL INPUTS
%
%  rxNum:              rx to plot - not yet equipped to do 'all'
%  source:             data source for diffpoint: 'cmb','ffbm'
%                      (default 'ffbm')
%  chflag:             0 = none
%                      1 = apply defaults (greyed out in plots)
%  goodplot:           field in ind struct to use for plotting
%                      'rgl','gl','l' (default 'rgl')
%  plotdiffpoint:      1 to plot diff pointing
%                      0 (default)
%  plotbw:             1 to plot beamwidth (uses plot_tiles but nicer)
%                      0 (default)
%  plotdiffbw:         1 to plot diff beamwidth
%                      0 (default)
%  plotdiffellip:      1 to plot diff ellipticity
%                      0 (default)
%  plotellip:          1 to plot total ellipticity (per-det)
%                      0 (default)
%  plotellipdiffellip: 1 to plot diff ellip as ellipse (UNCLEAR IF I WORK)
%                      0 (default)
%  makepng:            1 to make nice .pngs
%                      0 (default) to just make .eps
% 
% OUTPUTS
%  
%  A bunch of .eps/png files are saved in the local directory

expt = compopt.expt;
year = compopt.year;
if ~isfield(compopt,'rxNum')
  compopt.rxNum = 'all';
  rxNum = compopt.rxNum;
else
  rxNum = compopt.rxNum;
end
if ~isfield(compopt,'source')
  compopt.source = 'ffbm';
  source = compopt.source;
else
  source = compopt.source;
end
if ~isfield(compopt,'chflag')
  compopt.chflag = 1;
  chflag = compopt.chflag;
else
  chflag = compopt.chflag;
end
if ~isfield(compopt,'goodplot')
  compopt.goodplot = 'rgl';
  goodplot = compopt.goodplot;
else
  goodplot = compopt.goodplot;
end
if ~isfield(compopt,'plotdiffpoint')
  compopt.plotdiffpoint = 0;
  plotdiffpoint = compopt.plotdiffpoint;
else
  plotdiffpoint = compopt.plotdiffpoint;
end
if ~isfield(compopt,'plotbw')
  compopt.plotbw = 0;
  plotbw = compopt.plotbw;
else
  plotbw = compopt.plotbw;
end
if ~isfield(compopt,'plotdiffbw')
  compopt.plotdiffbw = 0;
  plotdiffbw = compopt.plotdiffbw;
else
  plotdiffbw = compopt.plotdiffbw;
end
if ~isfield(compopt,'plotellip')
  compopt.plotellip = 0;
  plotellip = compopt.plotellip;
else
  plotellip = compopt.plotellip;
end
if ~isfield(compopt,'plotdiffellip')
  compopt.plotdiffellip = 0;
  plotdiffellip = compopt.plotdiffellip;
else
  plotdiffellip = compopt.plotdiffellip;
end
if ~isfield(compopt,'plotellipdiffellip')
  compopt.plotellipdiffellip = 0;
  plotellipdiffellip = compopt.plotellipdiffellip;
else
  plotellipdiffellip = compopt.plotellipdiffellip;
end
if ~isfield(compopt,'makepng')
  compopt.makepng = 0;
  makepng = compopt.makepng;
else
  makepng = compopt.makepng;
end
if ~isfield(compopt,'subdir')
  compopt.subdir = './';
  subdir = compopt.subdir;
else
  subdir = compopt.subdir;
end

% First get pid/indid ('ideal') and then p/ind with chflags if requested
% For Keck/BICEP3, we call [p ind] with
% 'beamcen' = 'ideal', 'chi' = 'ideal', 'beamwid' = 'obs', 
% 'diffpoint' = 'ideal'
% So for beamwidth, it  calls up the 'beamwid' file meaning our FFBM
% measurements are in p.  Since 'beamwid' contains sigma/p/c,
% p.fwhm_maj,p.fwhm_min,p.alpha are also updated based on our beam
% measurements, which is why we don't have to read directly from
% 'beamwid' file.  Totally obfuscating.
switch expt
  case 'bicep2'
    [pid,indid] = get_array_info('20120303');
    if chflag
      flags = get_default_chflags(expt,2012);
    else
      flags = [];
    end
    [p,ind] = get_array_info('20120303','obs','obs','obs','obs',flags);
    pp = ParameterRead('aux_data/beams/beams_bicep2_obs_rwa_20130607.csv');
    p.p = pp.p;
    p.c = pp.c;
    p.sigma = pp.sigma;
  case 'keck'
    date = [num2str(year) '0201'];
    [pid,indid] = get_array_info(date);
    pid = rmfield(pid,'expt');
    if chflag
      flags = get_default_chflags(expt,year);
    else
      flags = [];
    end
    [p,ind] = get_array_info(date,'ideal','ideal','obs','ideal',flags);
    pp = ParameterRead(['aux_data/beams/beamwid_' num2str(year) ...
                        '0101.csv']);
    % The next 3 lines were commented out?
    %p.p = pp.p;
    %p.c = pp.c;
    %p.sig = pp.sigma;
    %
    p.err_sigma = pp.err_sigma;
    p.err_p = pp.err_p;
    p.err_c = pp.err_c;
    p.err_dx = pp.aboffset_err_x;
    p.err_dy = pp.aboffset_err_y;
    switch rxNum
      case 'all'
      otherwise
	cutind = (p.rx == rxNum);
        p = rmfield(p,'expt'); % for structcut to work
	p = structcut(p,cutind);
	ind = make_ind(p);
	pid = structcut(pid,cutind);
	indid = make_ind(pid);
	pp = structcut(pp,cutind);
    end
  case 'bicep3'
    date = [num2str(year) '0201'];
    [pid,indid] = get_array_info(date);
    pid = rmfield(pid,'expt');
    if chflag
      flags = get_default_chflags(expt,year);
    else
      flags = [];
    end
    [p,ind] = get_array_info(date,'ideal','ideal','obs','ideal',flags);
    pp = ParameterRead(['aux_data/beams/beamwid_' num2str(year) ...
                        '0101.csv']);
    p.err_sigma = pp.err_sigma;
    p.err_p = pp.err_p;
    p.err_c = pp.err_c;
    p.err_dx = pp.aboffset_err_x;
    p.err_dy = pp.aboffset_err_y;
end


% Typical CLW stream of consciousness:
% dkangle increases clockwise.
% To get the usual tile layout, we have to rotate dkangle = -90.
% But theta increases CCW.
% Drumangle is in the same direction as dkangle.
% What you want to do is remove drumangle from theta,
% so you add +drumangle to remove the drumangle
% and -dkangle since you want to rotate by dkangle

% Get ideal pixel centers (pid = "ideal") for plotting locations
% Rotate by 90 deg since +x' is UP in our standard layout
rotangle = -90;
switch expt
  case 'bicep2'
    pid.x = pid.r.*cosd(pid.theta - rotangle);
    pid.y = pid.r.*sind(pid.theta - rotangle);
  case {'keck','bicep3'}
    pid.x = pid.r.*cosd(pid.theta + pid.drumangle - rotangle);
    pid.y = pid.r.*sind(pid.theta + pid.drumangle - rotangle);
end

% Get measured beam locations for diff pointing
% If from ffbm, get from beams file (encoded in pp)
% If from CMB, get_array_info (encoded in p) automatically loads from
% diffpoint file
switch source
  case 'ffbm'
    % Values from beams file are in x'/y'
    p.x = pp.r.*cosd(pp.theta - rotangle);
    p.y = pp.r.*sind(pp.theta - rotangle);
  case 'cmb'
    % Values from get_array_info have drum angle baked in
    p.x = p.r.*cosd(p.theta + p.drumangle - rotangle);
    p.y = p.r.*sind(p.theta + p.drumangle - rotangle);
end

[p.sig, p.c, p.p] = egauss2_mmt2scp(p.fwhm_maj,p.fwhm_min,...
    p.alpha + p.theta + p.drumangle - rotangle);
p.ellip = sqrt(p.c.^2 + p.p.^2);
p.ellipanglefromx = p.alpha + p.theta + p.drumangle - rotangle; %CCW

p.dsig = NaN(size(p.gcp));
p.dx = NaN(size(p.gcp));
p.dy = NaN(size(p.gcp));
p.dp = NaN(size(p.gcp));
p.dc = NaN(size(p.gcp));


p.dsig(ind.a) = p.sig(ind.a) - p.sig(ind.b);
p.dx(ind.a) = p.x(ind.a) - p.x(ind.b);
p.dy(ind.a) = p.y(ind.a) - p.y(ind.b);
p.dp(ind.a) = p.p(ind.a) - p.p(ind.b);
p.dc(ind.a) = p.c(ind.a) - p.c(ind.b);

% KLUGE to not plot BICEP3 dets which have large diff parameters and
% which are not (yet) cut by the rgl list
switch expt
  case 'bicep3'
    p.dx(abs(p.dx) > 0.1) = NaN;
    p.dy(abs(p.dx) > 0.1) = NaN;
    p.dx(abs(p.dy) > 0.1) = NaN;
    p.dy(abs(p.dy) > 0.1) = NaN;
    p.dsig(abs(p.dsig) > 0.01) = NaN;
end

% Choose which indices to plot in a more general way
% Bad pixels (plotted but greyed out) are those in 'ideal' but not in
% index with cuts applied
switch goodplot
  case 'rgl'
    goodpix_a = ind.rgla;
    goodpix_b = ind.rglb;
    badpix_a = setdiff(indid.rgla,ind.rgla);
    badpix_b = setdiff(indid.rglb,ind.rglb);
  case 'gl'
    goodpix_a = ind.gla;
    goodpix_b = ind.glb;
    badpix_a = setdiff(indid.gla,ind.gla);
    badpix_b = setdiff(indid.glb,ind.glb);
  case 'l'
    goodpix_a = ind.la;
    goodpix_b = ind.lb;
    badpix_a = setdiff(indid.la,ind.la);
    badpix_b = setdiff(indid.lb,ind.lb);
end

figure('Visible','off') 

if plotdiffpoint
  % Good detectors in black
  clf; 
  switch expt
    case {'bicep2','keck'}
      setwinsize(gcf,400,400);
      scalefac = 20;
    case 'bicep3'
      setwinsize(gcf,800,800);
      scalefac = 40;
  end
  toplot = goodpix_a;
  
  q = quiver(pid.x(toplot),pid.y(toplot),...
      p.dx(toplot)*scalefac,p.dy(toplot)*scalefac,0,'o','filled');
  set(q,'MarkerSize',1);
  set(q,'LineWidth',0.5);
  set(q,'ShowArrowHead','on');
  set(q,'MaxHeadSize',1);
  set(q,'Color',[0,0,0]);
  
  % Bad detectors in grey
  toplot = badpix_a;
  hold on
  q = quiver(pid.x(toplot),pid.y(toplot),...
      p.dx(toplot)*scalefac,p.dy(toplot)*scalefac,0,'o','filled');
  set(q,'MarkerSize',1);
  set(q,'LineWidth',0.5);
  set(q,'ShowArrowHead','on');
  set(q,'MaxHeadSize',1);
  set(q,'Color',[0.6,0.6,0.6]);
  
  % Legend
  switch expt
    case {'bicep2','keck'}
      legx = sqrt(2^2/2)/60; % 5 arcmin legend
      q = quiver([-7.8],[7.3],[legx*scalefac],[legx*scalefac],0,'o','filled',...
                 'MarkerSize',1,'LineWidth',0.5,'ShowArrowHead','on',...
                 'MaxHeadSize',1,'Color','k');
      text(-7.3,7.5,'$\sqrt{\delta x^2+~\delta y^2}=~2~\mathrm{arcm~in}$',...
           'Interpreter','latex','Fontsize',10);
      q = quiver([7.1],[7.3],[-0.4],[0],0,...
                 'MarkerSize',1,'LineWidth',0.5,'ShowArrowHead','on',...
                 'MaxHeadSize',1,'Color','k');
      q = quiver([7.1],[7.3],[0],[0.4],0,...
                 'MarkerSize',1,'LineWidth',0.5,'ShowArrowHead','on',...
                 'MaxHeadSize',1,'Color','k');
      text(7.1,7.5,' +x'' ','Fontsize',7);
      text(6.5,7,' +y'' ','Fontsize',7);
    case 'bicep3'
      legx = sqrt(2^2/2)/60; % 5 arcmin legend
      q = quiver([-12.5],[-7.95],[legx*scalefac],[legx*scalefac],...
                 0,'o','filled',...
                 'MarkerSize',1,'LineWidth',0.5,'ShowArrowHead','on',...
                 'MaxHeadSize',1,'Color','k');
      text(-12,-7.9,'$\sqrt{\delta x^2+~\delta y^2}=~2~\mathrm{arcm~in}$',...
           'Interpreter','latex','Fontsize',10);
      q = quiver([-10],[-5],[-1],[0],0,...
                 'MarkerSize',1,'LineWidth',0.5,'ShowArrowHead','on',...
                 'MaxHeadSize',1,'Color','k');
      q = quiver([-10],[-5],[0],[1],0,...
                 'MarkerSize',1,'LineWidth',0.5,'ShowArrowHead','on',...
                 'MaxHeadSize',1,'Color','k');
      text(-10,-4.5,' +x'' ','Fontsize',9);
      text(-10.8,-5.3,' +y'' ','Fontsize',9);
  end
  
  % Plot size
  switch expt
    case {'bicep2','keck'}
      axis([-8 8 -8 8])
      axis square
      set(gca,'YTick',[-5,0,5]);
      set(gca,'XTick',[-5,0,5]);
      set(gca,'YTickLabel',char({'-5','0','5'}));
      set(gca,'XTickLabel',char({'-5','0','5'}));      
    case 'bicep3'
      axis([-13 13 -13 13])
      axis square
      set(gca,'YTick',[-10 -5 0 5 10]);
      set(gca,'XTick',[-10 -5 0 5 10]);
      set(gca,'YTickLabel',char({'-10','-5','0','5','10'}));
      set(gca,'XTickLabel',char({'-10','-5','0','5','10'}));
  end
  xlabel('Degrees');
  ylabel('Degrees');
  set(gca,'Box','on');
  grid on
  
  switch expt
    case 'bicep2'
      title(['BICEP2 Differential Pointing']);
      filename = [subdir '/diffpoint_' source '_' expt];
    case 'keck'
      title(['Keck ' num2str(year) ' rx' num2str(rxNum) ...
	    ' Differential Pointing'])
      filename = [subdir '/diffpoint_' source '_' ...
	    expt '_' num2str(year) '_rx' num2str(rxNum)];
    case 'bicep3'
      title(['BICEP3 ' num2str(year) ' Differential Pointing']);
      filename = [subdir '/diffpoint_' source '_' expt '_' num2str(year)];
  end

  print('-depsc2',filename);
  
  if makepng
    system(['rm ' filename '.png']);
    system(['eps2png -B ' filename '.eps']);
  end
  
end % Diffpoint

if plotbw

  sigma = pp.sigma;
  % Remove the automatically-inserted median
  sigma(sigma == nanmedian(sigma)) = NaN;

  clf; 
  % plot_tiles has problems with the default renderer
  set(0,'DefaultFigureRenderer','zbuffer')
  switch expt
    case {'bicep2','keck'}
      setwinsize(gcf,400,400);
    case 'bicep3'
      plot_tiles(sigma,p,...
                 'clim',[0.152 0.182],...
                 'clab','Sigma (degrees)',...
                 'title',['BICEP3 ' num2str(year) ' Beamwidth'])
      setwinsize(gcf,1000,800);
      filename = [subdir '/sigma_' expt '_' num2str(year)];
      print(gcf,filename,'-depsc2','-painters');
      %printfig(gcf,filename,'png');
  end

  if makepng
    system(['rm ' filename '.png']);
    system(['eps2png -B ' filename '.eps']);
  end
  
  set(0,'DefaultFigureRenderer','painters')
end % Beamwidth

if plotdiffbw
  
  clf; 
  switch expt
    case {'bicep2','keck'}
      setwinsize(gcf,400,400);
    case 'bicep3'
      setwinsize(gcf,800,800);
  end
  
  scale = 800;
  clear toplot;
  toplot{1} = goodpix_a;
  colr1{1} = [0 0 1]; % Blue
  colr2{1} = [1 0 0]; % Red
  toplot{2} = badpix_a;
  colr1{2} = [0.7 0.7 1]; % Light blue
  colr2{2} = [1 0.7 0.7]; % Light red
  
  for kk = 1:length(toplot)
    subplot{1} = find(p.dsig(toplot{kk}) > 0); % Blue: A > B
    subplot{2} = find(p.dsig(toplot{kk}) < 0);
    hold on
    for ii = 1:length(subplot{1})
      q = plot(pid.x(toplot{kk}(subplot{1}(ii))),...
	  pid.y(toplot{kk}(subplot{1}(ii))),...
	  'o','MarkerSize',abs(p.dsig(toplot{kk}(subplot{1}(ii))))*scale);
      set(q,'MarkerFaceColor',colr1{kk});
      set(q,'MarkerEdgeColor',colr1{kk});
    end
    for ii = 1:length(subplot{2})
      q = plot(pid.x(toplot{kk}(subplot{2}(ii))),...
	  pid.y(toplot{kk}(subplot{2}(ii))),...
	  'o','MarkerSize',abs(p.dsig(toplot{kk}(subplot{2}(ii))))*scale);
      set(q,'MarkerFaceColor',colr2{kk});
      set(q,'MarkerEdgeColor',colr2{kk});
    end
  end

  % Legend
  switch expt
    case {'bicep2','keck'}
      plot(-7.5,7.3,'ok','MarkerFaceColor','k','MarkerSize',0.005*scale);
      text(-7.2,7.2,'d\sigma = 0.005 degrees')
      q = quiver([7.1],[7.3],[-0.4],[0],0,...
                 'MarkerSize',1,'LineWidth',0.5,'ShowArrowHead','on',...
                 'MaxHeadSize',1,'Color','k');
      q = quiver([7.1],[7.3],[0],[0.4],0,...
                 'MarkerSize',1,'LineWidth',0.5,'ShowArrowHead','on',...
                 'MaxHeadSize',1,'Color','k');
      text(7.1,7.5,' +x'' ','Fontsize',7);
      text(6.5,7,' +y'' ','Fontsize',7);
    case 'bicep3'
      plot(-12,-7.9,'ok','MarkerFaceColor','k','MarkerSize',0.005*scale);
      text(-11.5,-7.9,'d\sigma = 0.005 degrees')
      q = quiver([-10],[-5],[-1],[0],0,...
                 'MarkerSize',1,'LineWidth',0.5,'ShowArrowHead','on',...
                 'MaxHeadSize',1,'Color','k');
      q = quiver([-10],[-5],[0],[1],0,...
                 'MarkerSize',1,'LineWidth',0.5,'ShowArrowHead','on',...
                 'MaxHeadSize',1,'Color','k');
      text(-10,-4.5,' +x'' ','Fontsize',9);
      text(-10.8,-5.3,' +y'' ','Fontsize',9);
  end

  % Plot size
  switch expt
    case {'bicep2','keck'}
      axis([-8 8 -8 8])
      axis square
      set(gca,'YTick',[-5,0,5]);
      set(gca,'XTick',[-5,0,5]);
      set(gca,'YTickLabel',char({'-5','0','5'}));
      set(gca,'XTickLabel',char({'-5','0','5'}));      
    case 'bicep3'
      axis([-13 13 -13 13])
      axis square
      set(gca,'YTick',[-10 -5 0 5 10]);
      set(gca,'XTick',[-10 -5 0 5 10]);
      set(gca,'YTickLabel',char({'-10','-5','0','5','10'}));
      set(gca,'XTickLabel',char({'-10','-5','0','5','10'}));
  end
  xlabel('Degrees');
  ylabel('Degrees');
  set(gca,'Box','on');
  grid on

  switch expt
    case 'bicep2'
      title(['BICEP2 Differential Beamwidth']);
      filename = [subdir '/dsig_' expt];
    case 'keck'
      title(['Keck ' num2str(year) ' rx' num2str(rxNum) ...
	    ' Differential Beamwidth'])
      filename = [subdir '/dsig_' ...
	    expt '_' num2str(year) '_rx' num2str(rxNum)];
    case 'bicep3'
      title(['BICEP3 ' num2str(year) ' Differential Beamwidth']);
      filename = [subdir '/dsig_' expt '_' num2str(year)];
  end

  print('-depsc2',filename)
  
  if makepng
    system(['rm ' filename '.png']);
    system(['eps2png -B ' filename '.eps']);
  end
 
end

if plotdiffellip

  clf; 
  switch expt
    case {'bicep2','keck'}
      setwinsize(gcf,400,400);
      scale = 0.2;
    case 'bicep3'
      setwinsize(gcf,800,800);
      scale = 0.15;
  end

  % Set up ellipses
  scale2 = 75;
  [a,b,alpha] = egauss2_scp2mmt(0.05,p.dc,p.dp);
  de = sqrt(p.dc.^2 + p.dp.^2);
  a = sqrt(1 + scale2*de);
  b = 1./a;
  [x,y] = ellipse(pid.x,pid.y,a*scale,b*scale,alpha*pi/180,100,1);
  
  % For the magnitude legend
  % Need different locations for B2/Keck/B3
  sizescale = 0.02;
  sizescalex = -12;
  sizescaley = -7.9;
  sizescaleangle = 45*ones(size(sizescale));
  a = sqrt(1 + scale2*sizescale);
  b = 1./a;
  [xx,yy] = ellipse(sizescalex,sizescaley,a*scale,b*scale,...
      sizescaleangle*pi/180,[],1);

  % For the direction legend - one +p, -p, +c, -c
  sizescale = 0.05;
  sizescalex = -12;
  sizescaley = -9.5;
  sizescaleangle = 90*ones(size(sizescale));
  a = sqrt(1 + scale2*sizescale);
  b = 1./a;
  [xx1,yy1] = ellipse(sizescalex,sizescaley,a*scale,b*scale,...
                      sizescaleangle*pi/180,[],1);
  
  sizescale = 0.05;
  sizescalex = -12;
  sizescaley = -10.5;
  sizescaleangle = zeros(size(sizescale));
  a = sqrt(1 + scale2*sizescale);
  b = 1./a;
  [xx2,yy2] = ellipse(sizescalex,sizescaley,a*scale,b*scale,...
                      sizescaleangle*pi/180,[],1);
  
  sizescale = 0.05;
  sizescalex = -12;
  sizescaley = -11.5;
  sizescaleangle = 135*ones(size(sizescale));
  a = sqrt(1 + scale2*sizescale);
  b = 1./a;
  [xx3,yy3] = ellipse(sizescalex,sizescaley,a*scale,b*scale,...
                      sizescaleangle*pi/180,[],1);
  
  sizescale = 0.05;
  sizescalex = -12;
  sizescaley = -12.5;
  sizescaleangle = 45*ones(size(sizescale));
  a = sqrt(1 + scale2*sizescale);
  b = 1./a;
  [xx4,yy4] = ellipse(sizescalex,sizescaley,a*scale,b*scale,...
                      sizescaleangle*pi/180,[],1);
  
  clear toplot;
  toplot{1} = goodpix_a;
  colr1{1} = [0 0 1]; % Blue
  colr2{1} = [1 0 0]; % Red
  toplot{2} = badpix_a;
  colr1{2} = [0.6 0.6 1];
  colr2{2} = [1 0.6 0.6];
  
  hold on
  for ii = 1:length(toplot)
    q = patch(x(:,toplot{ii}),y(:,toplot{ii}),'b');
    set(q,'LineWidth',0.5);
    set(q,'FaceColor',colr1{ii});
    set(q,'EdgeColor',colr1{ii});
    %q=plot(xx(:,toplot{ii}),yy(:,toplot{ii}),'r','linewidth',1);
  end
  
  % Legend
  switch expt
    case {'bicep2','keck'}
      q = patch(xx,yy,'k');
      set(q,'LineWidth',1);
      set(q,'FaceColor',[0,0,1]);
      set(q,'EdgeColor',[0,0,1]);
      text(-7.3,7.5,'$\sqrt{\delta p^2+\delta c^2}=~2\%$','Interpreter','latex');
      q = quiver([7.1],[7.3],[-0.4],[0],0,...
                 'MarkerSize',1,'LineWidth',0.5,'ShowArrowHead','on',...
                 'MaxHeadSize',1,'Color','k');
      q = quiver([7.1],[7.3],[0],[0.4],0,...
                 'MarkerSize',1,'LineWidth',0.5,'ShowArrowHead','on',...
                 'MaxHeadSize',1,'Color','k');
      text(7.1,7.5,' +x'' ','Fontsize',7);
      text(6.5,7,' +y'' ','Fontsize',7);
    case 'bicep3'
      q = patch(xx1,yy1,'k');
      set(q,'LineWidth',1);
      set(q,'FaceColor',[0,0,1]);
      set(q,'EdgeColor',[0,0,1]);
      text(-11.5,-9.5,'+p')
      q = patch(xx2,yy2,'k');
      set(q,'LineWidth',1);
      set(q,'FaceColor',[0,0,1]);
      set(q,'EdgeColor',[0,0,1]);
      text(-11.5,-10.5,'-p')
      q = patch(xx3,yy3,'k');
      set(q,'LineWidth',1);
      set(q,'FaceColor',[0,0,1]);
      set(q,'EdgeColor',[0,0,1]);
      text(-11.5,-11.5,'+c')
      q = patch(xx4,yy4,'k');
      set(q,'LineWidth',1);
      set(q,'FaceColor',[0,0,1]);
      set(q,'EdgeColor',[0,0,1]);
      text(-11.5,-12.5,'-c')
      q = patch(xx,yy,'k');
      set(q,'LineWidth',1);
      set(q,'FaceColor',[0,0,1]);
      set(q,'EdgeColor',[0,0,1]);
      text(-11.5,-7.9,'$\sqrt{\delta p^2+\delta c^2}=~2\%$','Interpreter','latex');
      q = quiver([-10],[-5],[-1],[0],0,...
                 'MarkerSize',1,'LineWidth',0.5,'ShowArrowHead','on',...
                 'MaxHeadSize',1,'Color','k');
      q = quiver([-10],[-5],[0],[1],0,...
                 'MarkerSize',1,'LineWidth',0.5,'ShowArrowHead','on',...
                 'MaxHeadSize',1,'Color','k');
      text(-10,-4.5,' +x'' ','Fontsize',9);
      text(-10.8,-5.3,' +y'' ','Fontsize',9);
  end

  % Plot size
  switch expt
    case {'bicep2','keck'}
      axis([-8 8 -8 8])
      axis square
      set(gca,'YTick',[-5,0,5]);
      set(gca,'XTick',[-5,0,5]);
      set(gca,'YTickLabel',char({'-5','0','5'}));
      set(gca,'XTickLabel',char({'-5','0','5'}));      
    case 'bicep3'
      axis([-13 13 -13 13])
      axis square
      set(gca,'YTick',[-10 -5 0 5 10]);
      set(gca,'XTick',[-10 -5 0 5 10]);
      set(gca,'YTickLabel',char({'-10','-5','0','5','10'}));
      set(gca,'XTickLabel',char({'-10','-5','0','5','10'}));
  end
  xlabel('Degrees');
  ylabel('Degrees');
  set(gca,'Box','on');
  grid on
  
  switch expt
    case 'bicep2'
      title(['BICEP2 Differential Ellipticity']);
      filename = [subdir '/dellip_' expt];
    case 'keck'
      title(['Keck ' num2str(year) ' rx' num2str(rxNum) ...
	    ' Differential Ellipticity'])
      filename = [subdir '/dellip_' ...
	    expt '_' num2str(year) '_rx' num2str(rxNum)];
    case 'bicep3'
      title(['BICEP3 ' num2str(year) ' Differential Ellipticity']);
      filename = [subdir '/dellip_' expt '_' num2str(year)];
  end

  print('-depsc2',filename)

  if makepng
    system(['rm ' filename '.png']);
    system(['eps2png -B ' filename '.eps']);
  end
  
end

if plotellip

  clf;
  switch expt
    case {'bicep2','keck'}
      setwinsize(gcf,400,400);
      scale = 0.2;
    case 'bicep3'
      setwinsize(gcf,800,800);
      scale = 0.15;
  end
  
  scale2 = 20;
  a = sqrt(1 + scale2*p.ellip);
  b = 1./a;
  %b=ones(size(p.gcp))*0.07;
  [x,y] = ellipse(pid.x,pid.y,a.*scale,b.*scale,...
      p.ellipanglefromx*pi/180,[],1);
  
  clf; 
  clear toplot;
  toplot{1} = goodpix_a;
  colr{1} = [1 0 0]; % Red
  toplot{2} = goodpix_b;
  colr{2} = [0 0 1]; % Blue
  toplot{3} = badpix_a;
  colr{3} = [1 0.6 0.6];
  toplot{4} = badpix_b;
  colr{4} = [0.6 0.6 1];
  
  hold on
  for ii = 1:length(toplot)
    q = plot(x(:,toplot{ii}),y(:,toplot{ii}));
    set(q,'LineWidth',0.5);
    set(q,'Color',colr{ii});    
  end
  
  % Legend
  switch expt
    case {'bicep2','keck'}
      a = sqrt(1 + scale2*0.05);
      b = 1./a;
      [x,y] = ellipse(-7.6,7.5,a.*scale,b.*scale,45*pi/180,[],1);
      q = plot(x,y);
      set(q,'MarkerSize',1);
      set(q,'LineWidth',1);
      set(q,'Color',[0,0,0]);
      %text(-7,7.5,'sqrt(p^2+c^2)= 5%');
      text(-7.3,7.5,'$\sqrt{p^2+c^2}=~5\%$','Interpreter','latex');
      q = quiver([7.1],[7.3],[-0.4],[0],0,...
                 'MarkerSize',1,'LineWidth',0.5,'ShowArrowHead','on',...
                 'MaxHeadSize',1,'Color','k');
      q = quiver([7.1],[7.3],[0],[0.4],0,...
                 'MarkerSize',1,'LineWidth',0.5,'ShowArrowHead','on',...
                 'MaxHeadSize',1,'Color','k');
      text(7.1,7.5,' +x'' ','Fontsize',7);
      text(6.5,7,' +y'' ','Fontsize',7);
    case 'bicep3'
      a = sqrt(1 + scale2*0.05);
      b = 1./a;
      [x,y] = ellipse(-12,-7.9,a.*scale,b.*scale,45*pi/180,[],1);
      q = plot(x,y);
      set(q,'MarkerSize',1);
      set(q,'LineWidth',1);
      set(q,'Color',[0,0,0]);
      text(-11.5,-7.9,'$\sqrt{p^2+c^2}=~5\%$','Interpreter','latex');
      q = quiver([-10],[-5],[-1],[0],0,...
                 'MarkerSize',1,'LineWidth',0.5,'ShowArrowHead','on',...
                 'MaxHeadSize',1,'Color','k');
      q = quiver([-10],[-5],[0],[1],0,...
                 'MarkerSize',1,'LineWidth',0.5,'ShowArrowHead','on',...
                 'MaxHeadSize',1,'Color','k');
      text(-10,-4.5,' +x'' ','Fontsize',9);
      text(-10.8,-5.3,' +y'' ','Fontsize',9);
  end
  
  % Plot size
  switch expt
    case {'bicep2','keck'}
      axis([-8 8 -8 8])
      axis square
      set(gca,'YTick',[-5,0,5]);
      set(gca,'XTick',[-5,0,5]);
      set(gca,'YTickLabel',char({'-5','0','5'}));
      set(gca,'XTickLabel',char({'-5','0','5'}));      
    case 'bicep3'
      axis([-13 13 -13 13])
      axis square
      set(gca,'YTick',[-10 -5 0 5 10]);
      set(gca,'XTick',[-10 -5 0 5 10]);
      set(gca,'YTickLabel',char({'-10','-5','0','5','10'}));
      set(gca,'XTickLabel',char({'-10','-5','0','5','10'}));
  end
  xlabel('Degrees');
  ylabel('Degrees');
  set(gca,'Box','on');
  grid on

  switch expt
    case 'bicep2'
      title(['BICEP2 Total Ellipticity']);
      filename = [subdir '/ellip_' expt];
    case 'keck'
      title(['Keck ' num2str(year) ' rx' num2str(rxNum) ...
	    ' Total Ellipticity'])
      filename = [subdir '/ellip_' ...
	    expt '_' num2str(year) '_rx' num2str(rxNum)];
    case 'bicep3'
      title(['BICEP3 ' num2str(year) ' Total Ellipticity']);
      filename = [subdir '/ellip_' expt '_' num2str(year)];
  end

  print('-depsc2',filename)

  if makepng
    system(['rm ' filename '.png']);
    system(['eps2png -B ' filename '.eps']);
  end
  
end


if plotellipdiffellip
  
  clf; setwinsize(gcf,400,700)
  subplot_grid(2,1,1);
  scale = 0.2;
  scale2 = 20;
  a = sqrt(1 + scale2*p.ellip);
  b = 1./a;
  %b=ones(size(p.gcp))*0.07;
  [x,y] = ellipse(pid.x,pid.y,a.*scale,b.*scale,...
      p.ellipanglefromx*pi/180,[],1);
  
  clear toplot;
  toplot{1} = goodpix_a;
  colr{1} = [0 0 1]; % Blue
  toplot{2} = goodpix_b;
  colr{2} = [1 0 0]; % Red
  toplot{3} = badpix_a;
  colr{3} = [0.6 0.6 1];
  toplot{4} = badpix_b;
  colr{4} = [1 0.6 0.6];
  
  hold on
  for ii = 1:length(toplot)
    q = plot(x(:,toplot{ii}),y(:,toplot{ii}));
    set(q,'LineWidth',0.5);
    set(q,'Color',colr{ii});    
  end
  xlabel('Degrees');
  ylabel('Degrees');
  axis([-8 8 -8 8]);
  axis square
  set(gca,'YTick',[-5,0,5]);
  set(gca,'XTick',[-5,0,5]);
  set(gca,'YTickLabel',char({'-5','0','5'}));
  set(gca,'XTickLabel',char({'-5','0','5'}));
  grid on
  
  % legend
  a = sqrt(1 + scale2*0.05);
  b = 1./a;
  [x,y] = ellipse(-7.6,7.5,a.*scale,b.*scale,45*pi/180,[],1);
  q=plot(x,y);
  set(q,'MarkerSize',1);
  set(q,'LineWidth',1);
  set(q,'Color',[0,0,0]);
  %text(-7,7.5,'sqrt(p^2+c^2)= 5%');
  text(-7.3,7.5,'$\sqrt{p^2+c^2}=~5\%$','Interpreter','latex');
  q = quiver([7.1],[7.3],[-0.4],[0],0,...
      'MarkerSize',1,'LineWidth',0.5,'ShowArrowHead','on',...
      'MaxHeadSize',1,'Color','k');
  q = quiver([7.1],[7.3],[0],[0.4],0,...
      'MarkerSize',1,'LineWidth',0.5,'ShowArrowHead','on',...
      'MaxHeadSize',1,'Color','k');
  text(7.1,7.5,' +x'' ','Fontsize',7);
  text(6.5,7,' +y'' ','Fontsize',7);

  subplot_grid2(2,1,1);
  subplot_grid(2,1,2);

  scale = 0.2;
  scale2 = 75;
  [a,b,alpha] = egauss2_scp2mmt(0.05,p.dc,p.dp);
  de = sqrt(p.dc.^2 + p.dp.^2);
  a = sqrt(1 + scale2*de);
  b = 1./a;
  [x,y] = ellipse(pid.x,pid.y,a*scale,b*scale,alpha*pi/180,100,1);
  
  %sizescale=[0.02 0.05 0.10];
  %sizescalex=[-7.6 -7 -6.4];
  %sizescaley=[7.3 7.3 7.3];
  sizescale = 0.02;
  sizescalex = -7.6;
  sizescaley = 7.5;
  sizescaleangle = 45*ones(size(sizescale));
  a = sqrt(1 + scale2*sizescale);
  b = 1./a;
  [xx,yy] = ellipse(sizescalex,sizescaley,a*scale,b*scale,...
      sizescaleangle*pi/180,[],1);
  
  clear toplot;
  toplot{1} = goodpix_a;
  colr1{1} = [0 0 1]; % Blue
  colr2{1} = [1 0 0]; % Red
  toplot{2} = badpix_a;
  colr1{2} = [0.6 0.6 1];
  colr2{2} = [1 0.6 0.6];
    
  hold on
  for ii = 1:length(toplot)
    q = patch(x(:,toplot{ii}),y(:,toplot{ii}),'b');
    set(q,'LineWidth',0.5);
    set(q,'FaceColor',colr1{ii});
    set(q,'EdgeColor',colr1{ii});
    %q=plot(xx(:,toplot{ii}),yy(:,toplot{ii}),'r','linewidth',1);
  end
  xlabel('Degrees');
  %ylabel('Degrees');
  axis([-8 8 -8 8]);
  axis square
  set(gca,'YTick',[-5,0,5]);
  set(gca,'XTick',[-5,0,5]);
  set(gca,'YTickLabel',char({'-5','0','5'}));
  set(gca,'XTickLabel',char({'-5','0','5'}));
  grid on
  
  % Legend
  q=patch(xx,yy,'k');
  set(q,'LineWidth',1);
  set(q,'FaceColor',[0,0,1]);
  set(q,'EdgeColor',[0,0,1]);
  text(-7.3,7.5,'$\sqrt{dp^2+dc^2}=~2\%$','Interpreter','latex');
  
  subplot_grid2(2,1,2);

  switch expt
    case 'bicep2'
      title(['BICEP2 Ellipticity']);
      filename = [subdir '/ellipellip_' expt];
    case 'keck'
      title(['Keck ' num2str(year) ' rx' num2str(rxNum) ...
	    ' Ellipticity'])
      filename = [subdir '/ellipellip_' ...
	    expt '_' num2str(year) '_rx' num2str(rxNum)];
  end

  print('-depsc2',filename)

  if makepng
    system(['rm ' filename '.png']);
    system(['eps2png -B ' filename '.eps']);
  end

end


return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dim = get_paneldims()
% Parameters for six-panel layout

% Overall width, in inches.
dim.W = 7.3; 

% Small gap, in inches.
%dim.thin = 0.05;
% Medium gap, in inches.
%dim.med = 0.24; 
% Wide gap, in inches.
%dim.wide = 0.4; 
% Width of color bar, in inches.
%dim.cbar = 0.15; 
dim.med = 0.25;
dim.wide = 0.5;
% Map width, in inches.
dim.mapw = (dim.W - 2 * dim.wide - 1 * dim.med) / 2;
% Map height, in inches.
dim.maph = dim.mapw; 
% Overall height, in inches.
dim.H = 3 * dim.maph + 2*dim.wide + 2 * dim.wide;
% Left edge of column 1.
dim.x1 = dim.wide / dim.W; 
% Left edge of column 2.
dim.x2 = (dim.wide + dim.mapw + dim.med) / dim.W; 
% Left edge of column 3 (colorbar).
%dim.x3 = (dim.wide + 2 * dim.mapw + 2 * dim.thin) / dim.W; 
dim.x3=(dim.x1+dim.x2)/2;
% Bottom edge of row 1.
dim.y1 = (dim.wide + 2 * dim.maph + 2 * dim.wide) / dim.H; 
% Bottom edge of row 2.
dim.y2 = (dim.wide + 1 * dim.maph + 1 * dim.wide) / dim.H;
% Bottom edge of row 3.
dim.y3 = dim.wide / dim.H; 
% Map width as fraction of image width.
dim.xmap = dim.mapw / dim.W; 
% Map height as fraction of image height.
dim.ymap = dim.maph / dim.H; 

dim.x = [dim.x1 dim.x2 dim.x1 dim.x2 dim.x3];
dim.y = [dim.y1 dim.y1 dim.y2 dim.y2 dim.y3];

return