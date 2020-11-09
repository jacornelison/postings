function plot_map(m,map,convention,cropwindow,nancolor,proj)
% plot_map(m,map,convention,cropwindow,nancolor,proj)
%
% plot a map with correct aspect ratio and "sphere viewed from the inside"
% If proj is set, plot with different projections.
% The following projections (except for 'mollweid') are centered
%  at the center of the observing field defined by m structure.
%   'radec'    : (default) radec projection
%   Below are plotted with pcolorm(). To save the plot, use recent version of  matlab. 2009a version doesn't work.
%   'ortho'    : Orthographic projection. Originally "sphere viewed from the infinite distance outside"
%                 But reversed RA direction so that RA increases to the left just as we observe the sphere from the inside.
%   'gnomonic' : Gnomonic projection. "sphere viewed from the inside at the center"
%   'mollweid' : Mollweide projection.
%   'lambert'  : Lambert projection.

if(~exist('convention','var') || isempty(convention))
  convention='iau';
end

if(~exist('nancolor','var') || isempty(nancolor))
  nancolor=[1,1,1];
end

if(~exist('cropwindow','var') || isempty(cropwindow))
  cropwindow=false;
elseif(length(cropwindow)~=4)
  cropwindow=[-30,30,-66,-49];
end

if(~exist('proj','var') || isempty(proj))
  proj='radec';
end


if any(cropwindow)
  % crop down to the specified ra/dec
  xl=find(m.x_tic<cropwindow(1),1,'last');
  xh=find(m.x_tic>cropwindow(2),1,'first');
  yl=find(m.y_tic<cropwindow(3),1,'last');
  yh=find(m.y_tic>cropwindow(4),1,'first');
  xr=range(m.x_tic);
  yr=range(m.y_tic);
  m.x_tic=m.x_tic(xl:xh);
  m.y_tic=m.y_tic(yl:yh);
  m.xdos=m.xdos*range(m.x_tic)/xr;
  m.ydos=m.ydos*range(m.y_tic)/yr;
  map_tmp=map(yl:yh,xl:xh);
  clear map
  % taper at the very edge
  taper=0.1*ones(length(m.y_tic),length(m.x_tic));
  taper(3:(length(m.y_tic)-2),3:(length(m.x_tic)-2))=1;
  taper=conv2(taper,ones(5,5),'same');
  taper=taper+flipud(taper); taper=taper+fliplr(taper);
  taper=taper/nanmedian(taper(:));
  map=map_tmp.*taper;
end

cla % this needed to avoid weird behavior..
switch proj
 case 'radec'
  if verLessThan('matlab', '8.4')
    imagescnan(m.x_tic,m.y_tic,full(map),'nancolor', nancolor);
  else
    % NaN handling much easier/better after graphics rewrite in R2014b
    im = imagesc(m.x_tic, m.y_tic, full(map));
    set(im, 'AlphaData', ~isnan(map));
  end
  if(isfield(m,'xdos'))
    xdos=m.xdos; ydos=m.ydos;
  else
    xdos=(m.x_tic(end)-m.x_tic(1))*cosd(mean(m.y_tic));
    ydos=(m.y_tic(end)-m.y_tic(1));
  end
  set(gca,'PlotBoxAspectRatio',[xdos,ydos,1])
  axis xy;
  switch convention
   case 'iau'
    % north is up
    set(gca,'XDir','reverse');
   case 'spole'
    % sky as seen by someone standing at South Pole
    set(gca,'YDir','reverse');
  end
  xlabel('RA (deg)'); ylabel('Dec (deg)');

 case 'ortho'
  [lon,lat]=meshgrid(m.x_tic,m.y_tic); % in deg
  axlim = [m.ly m.hy  m.lx m.hx];
  meanx=mean([m.lx,m.hx]);
  meany=mean([m.ly,m.hy]);
  axesm('MapProjection','ortho','Origin',[meany,meanx])
  axis off
  pcolorm(lat,lon,map)
  % Make sure it is the view from the inside the sphere (RA increases to the left of the figure)
  %  (assumes RA range is from -180 to +180 and field doesn't overlap RA = +/- 180)
  if m.x_tic(1)<m.x_tic(end)
    set(gca,'XDir','reverse');
  end
  %setwinsize(gcf(), 700, 700*(m.ydos/m.xdos));
  % Add grid lines that correcpond to RA/Dec lines
  ragrids=[-180:30:180];ragrids=ragrids(ragrids>=m.lx-15 & ragrids<=m.hx+15);
  decgrids=[-90:10:90];decgrids=decgrids(decgrids>=m.ly-5 & decgrids<=m.hy+5);
  for raidx=1:length(ragrids)
    % draw constant ra grid line
    x=ones(1,200)*ragrids(raidx);
    y=linspace(decgrids(1),decgrids(end),200);
    plotm(y,x,':')
    textm(y(1)-0.05*(y(end)-y(1)),x(1),[num2str(ragrids(raidx)) char(176)])
  end
  for decidx=1:length(decgrids)
    % draw constant dec grid line
    x=linspace(ragrids(1),ragrids(end),200);
    y=ones(1,200)*decgrids(decidx);
    plotm(y,x,':')
    textm(y(1),x(end)+0.05*(x(end)-x(1)),[num2str(decgrids(decidx)) char(176)])
  end

 case 'gnomonic'
  [lon,lat]=meshgrid(m.x_tic,m.y_tic); % in deg
  axlim = [m.ly m.hy  m.lx m.hx];
  meanx=mean([m.lx,m.hx]);
  meany=mean([m.ly,m.hy]);
  axesm('MapProjection','gnomonic','Origin',[meany,meanx])
  axis off
  pcolorm(lat,lon,map)
  % Make sure it is the view from the inside the sphere (RA increases to the left of the figure)
  %  (assumes RA range is from -180 to +180 and field doesn't overlap RA = +/- 180)
  if m.x_tic(1)<m.x_tic(end)
    set(gca,'XDir','reverse');
  end
  %setwinsize(gcf(), 700, 700*(m.ydos/m.xdos));
  % Add grid lines that correcpond to RA/Dec lines
  ragrids=[-180:30:180];ragrids=ragrids(ragrids>=m.lx-15 & ragrids<=m.hx+15);
  decgrids=[-90:10:90];decgrids=decgrids(decgrids>=m.ly-5 & decgrids<=m.hy+5);
  for raidx=1:length(ragrids)
    % draw constant ra grid line
    x=ones(1,200)*ragrids(raidx);
    y=linspace(decgrids(1),decgrids(end),200);
    plotm(y,x,':')
    textm(y(1)-0.05*(y(end)-y(1)),x(1),[num2str(ragrids(raidx)) char(176)])
  end
  for decidx=1:length(decgrids)
    % draw constant dec grid line
    x=linspace(ragrids(1),ragrids(end),200);
    y=ones(1,200)*decgrids(decidx);
    plotm(y,x,':')
    textm(y(1),x(end)+0.05*(x(end)-x(1)),[num2str(decgrids(decidx)) char(176)])
  end

 case 'mollweid'
  % Note: for now, this works for traditional map definition structure. Not applicable for healpix map definition, yet.
  % generate matrix grid
  [lon, lat]=meshgrid(m.x_tic,m.y_tic); % in deg
  xwidth=m.hx-m.lx; ywidth=m.hy-m.ly;
  axlim = [m.ly-0.06*ywidth m.hy+0.06*ywidth  m.lx-0.05*xwidth m.hx+0.05*xwidth];
  axesm('MapProjection','mollweid','Grid','on','GLineStyle',':','GLineWidth',0.1,...
    'MapLatLimit',[axlim(1),axlim(2)],'MapLonLimit',[axlim(3),axlim(4)],...
    'MLabelParallel',m.hy+0.05*ywidth,'PLabelMeridian',axlim(4),...
    'ParallelLabel','on','MeridianLabel','on','LabelFormat','none');
  axis off
  pcolorm(lat, lon, map)
  % Make sure it is the view from the inside the sphere (RA increases to the left of the figure)
  %  (assumes RA range is from -180 to +180 and field doesn't overlap RA = +/- 180)
  if m.x_tic(1)<m.x_tic(end)
    set(gca,'XDir','reverse');
  end
  %setwinsize(gcf(), 700, 700*(m.ydos/m.xdos));

 case 'lambert'
  % Note: for now, this works for traditional map definition structure. Not applicable for healpix map definition, yet.
  % generate matrix grid
  [lon, lat]=meshgrid(m.x_tic,m.y_tic); % in deg
  % set axis limit a bit wider than the map
  xwidth=m.hx-m.lx; ywidth=m.hy-m.ly;
  axlim = [m.ly-0.06*ywidth m.hy+0.06*ywidth  m.lx-0.05*xwidth m.hx+0.05*xwidth];
  axesm('MapProjection','lambert','MapParallels',[],'Grid','on','GLineStyle',':','GLineWidth',0.1,...
    'MapLatLimit',[axlim(1),axlim(2)],'MapLonLimit',[axlim(3),axlim(4)],...
    'MLabelParallel',m.hy+0.05*ywidth,'PLabelMeridian',axlim(4),...
    'ParallelLabel','on','MeridianLabel','on','LabelFormat','none')
  axis off
  pcolorm(lat,lon,map);
  % Make sure it is the view from the inside the sphere (RA increases to the left of the figure)
  %  (assumes RA range is from -180 to +180 and field doesn't overlap RA = +/- 180)
  if m.x_tic(1)<m.x_tic(end)
    set(gca,'XDir','reverse');
  end
  %setwinsize(gcf(), 700, 700*(m.ydos/m.xdos));
end

return
