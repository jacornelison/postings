function lines = plot_pol(m,Q,U,Qvar,Uvar,scalefac,noicut,downsamp,showlegend,cropwindow,noplot)
% plot_pol(m,Q,U,Qvar,Uvar,scalefac,noicut,downsamp,noplot)
%
% Plot polarization pseudo vectors
%
% scalefac scales size of vectors - select to task
% noicut selects region of map over which vectors will be plotted - up
%   to noicut times the median noise
% downsamp selects the downsampling of input map - vector plots are
%   way too cluttered if we have one vector for every map pixel - need
%   to downsample to way fewer - however this implies that the input maps
%   need to have been smoothed or Fourier filtered before input or we
%   will be aliasing
% noplot: don't plot but hand back the vectors as lines
%
% e.g.
% [m,realn]=reduc_plotebmap('1450/real_a_filtp3_weight3_gs_dp1100_jack0',[],'normal');
% plot_pol(m,realn.Q,realn.U,realn.Qvar,realn.Uvar,0.2);

if(~exist('scalefac','var'))
  scalefac=[];
end
if(~exist('noicut','var'))
  noicut=[];
end
if(~exist('downsamp','var'))
  downsamp=[];
end
if(~exist('showlegend','var'))
  showlegend=true;
end
if(~exist('cropwindow','var'))
  cropwindow=false;
end
if(~exist('noplot','var'))
  noplot=false;
end


if(isempty(scalefac))
  scalefac=1;
end
if(isempty(noicut))
  noicut=2;
end
if(isempty(downsamp))
  downsamp=5;
end
if(any(cropwindow)&&length(cropwindow)~=4)
  cropwindow=[-30,30,-66,-49];
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
  m.lx=m.x_tic(1)-mean(diff(m.x_tic))/2;
  m.hx=m.x_tic(end)+mean(diff(m.x_tic))/2;
  m.ly=m.y_tic(1)-mean(diff(m.y_tic))/2;
  m.hy=m.y_tic(end)+mean(diff(m.y_tic))/2;
  Q_tmp=Q(yl:yh,xl:xh);
  U_tmp=U(yl:yh,xl:xh);
  Qvar_tmp=Qvar(yl:yh,xl:xh);
  Uvar_tmp=Uvar(yl:yh,xl:xh);
  clear Q U Qvar Uvar
  % taper at the very edge
  taper=0.1*ones(length(m.y_tic),length(m.x_tic));
  taper(3:(length(m.y_tic)-2),3:(length(m.x_tic)-2))=1;
  taper=conv2(taper,ones(5,5),'same');
  taper=taper+flipud(taper); taper=taper+fliplr(taper);
  taper=taper/nanmedian(taper(:));
  Q=Q_tmp.*taper;
  U=U_tmp.*taper;
  Qvar=Qvar_tmp.*taper;
  Uvar=Uvar_tmp.*taper;
end

% initially adapted from code in reduc_plotcomap plot type 17

% downsample map - see comments above
m.x_tic=m.x_tic(1:downsamp:end);
m.y_tic=m.y_tic(1:downsamp:end);
Q=Q(1:downsamp:end,1:downsamp:end); U=U(1:downsamp:end,1:downsamp:end);
Qvar=Qvar(1:downsamp:end,1:downsamp:end); Uvar=Uvar(1:downsamp:end,1:downsamp:end);

% get expected noise as a function of position across map
pol_std=sqrt(Qvar+Uvar);

if ~cropwindow
  % instead of making a s/n cut just cut to center region of map where
  % noise is not much higher than minimum
  mask=pol_std<noicut*nanmedian(pol_std(:));
else
  % only trim the edges
  mask=true(size(Q));
  mask(1,:)=false; mask(end,:)=false; mask(:,1)=false; mask(:,end)=false;
  mask=mask(:);
end

% make grids of pol vector center locations
[xx,yy]=meshgrid(m.x_tic,m.y_tic);

% get magnitude and angle of pol signal
r=sqrt(Q.^2+U.^2);
th=atan2(U,-Q); % Jamie T code says -Q for IAU

% cut down to desired area (become vector lists)
xx=xx(mask); yy=yy(mask);
r=r(mask); th=th(mask);

if showlegend
  % append an extra length 1 vector in the corner as legend
  lx=m.x_tic(floor(length(m.x_tic)*0.12));
  ly=m.y_tic(floor(length(m.y_tic)*0.92));
  xx=[xx;lx];
  yy=[yy;ly];
  r=[r;1./scalefac]; th=[th;pi/2];
end

% calc delta x/y offsets for end of polvecs
dx=r.*cos(th/2)*scalefac;
dy=r.*sin(th/2)*scalefac;

% these offsets are correct for "axis equal" - below we will instead
% set the plot box aspect ratio to enforce square
% pixels - if we want the pol vector to be correctly oriented we need
% to rescale x by the ratio of the natural to forced box aspect ratios
dx=dx*((m.hx-m.lx)/(m.hy-m.ly))/(m.xdos/m.ydos);

% plot the pol vectors
lines = [];
if noplot
  lines = zeros(size(xx,1),2,20);
  for ii=1:size(xx,1)
    lines(ii,1,:)=linspace(xx(ii)-dx(ii),xx(ii)+dx(ii),20);
    lines(ii,2,:)=linspace(yy(ii)-dy(ii),yy(ii)+dy(ii),20);
  end
else
  plot([xx-dx xx+dx]', [yy-dy yy+dy]', 'k-');
  % set the axes such that the pixels - and grid of plotted pol vecs -
  % are square
  xlim([min(m.x_tic),max(m.x_tic)]); ylim([min(m.y_tic),max(m.y_tic)]); 
  axis xy; set(gca,'XDir','reverse');
  set(gca,'PlotBoxAspectRatio',[m.xdos,m.ydos,1])

  if showlegend
    % label the legend
  %    text(lx-m.xdos/40,m.y_tic(floor(length(m.y_tic)*0.94)),sprintf('%1.1f\\muK',1/scalefac),'FontSize',7,'FontWeight','normal');
    text(lx-m.xdos/40,m.y_tic(floor(length(m.y_tic)*0.94)),sprintf('%1.1f\\muK',1/scalefac),'FontWeight','normal');
  end
end

return
