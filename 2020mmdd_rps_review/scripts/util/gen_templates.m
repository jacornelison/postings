function T = gen_templates(hmap,d,p,ind,mapind,dp,residbeam,mapopt,simopt)
% gen_templates(hmap,d,p,ind,mapind,dp,residbeam,mapopt,simopt)
%
% hmap - healpix map to be interpolated from, loaded in variable or filename
% d - TOD data structure
% p,ind - output from get array info
% mapind - mapind=make_mapind(d,fs)
% dp - which templates to return:
%      1-relgain; 2-A/Bx; 3-A/By; 4-beamwid; 5-e_p; 6-e_c; 7-xtalk_up; 8-xtalk_down
% residbeam - used for constructing a residual beam mismatch template
%             structure containing unsmoothed skymaps and beammaps
% mapopt - several switches can be thrown:
%          deprojtest, curveskyrotbeam, curveskyrescale
% simopt - optional input; if deprojtest switch in
%          mapopt thrown use p.r and p.theta from simopt; necessary to use same
%          healpix pixel list for template generation as for simulation

disp('gen_templates')

if ischar(hmap)
  hmap=read_fits_map(hmap);
end

% are we doing crosstalk?
if any(dp==7 | dp==8)
  docrosstalk=1;
else
  docrosstalk=0;
end

n_templates = numel(dp);
% are we doing residual beam?
n_templates_sub = 6;
if any(dp==9)
  n_templates_sub = 9;  
end
if any(dp==10)
  n_templates_sub = 10;
end

% Only do the a's since the templates are per pair
n=length(ind.a);

% Cut down to only mapped points
ra = d.pointing.cel.ra(mapind);
dec = d.pointing.cel.dec(mapind);
dk = d.pointing.cel.dk(mapind);
dk_hor = d.pointing.hor.dk(mapind);

% initialize
T=zeros(length(d.pointing.cel.ra),n,n_templates);
temp=zeros(numel(ra),n,n_templates_sub);

clear d

% The derivatives of the healpix coordinate system (theta,phi) w.r.t. the focal plane
% coordinate system (x,y) only depends dec and dk, which are supposed to be
% constant. If the jitter in dec and dk are small enough, we can use their mean value
% and massively speed up the calculation of the template. This has the effect of making
% the derivatives dtheta/dx, dphi/dx, etc., scalars instead of vectors.
% This is turned off (though it does work well.)
%if std(dec)/mean(dec) < 1e-2 & std(dk)/mean(dk) <1e-2
%  disp('ignoring dk and dec jitter');
%  ignore_jitter=1;
%else
%  ignore_jitter=0;
%end
ignore_jitter=0;

% CONSTRUCT TEMPLATES
% order of temp will always be relgain, diffpoint-x, diffpoint-y, beamwid, ellip-plus,
% ellip-cross, xtalk-up, xtalk-down  
for i=1:n
  % Loop over channels. This is nearly as fast as doing it all at once and gets the
  % memory usage down to a reasonable level
  
  % expand the pointing
  ch=ind.a(i);
  if isfield(mapopt,'deprojtest') && mapopt.deprojtest;
    % Special test mode...
    % Force A/B to have common centroid in the manner of reduc_makesim
    chanind=[find(ind.b==ch),find(ind.a==ch)];
    chans0=[ind.a(chanind),ind.b(chanind)];
    n0=length(chans0);
    % find ra/dec trajectory for this channel
    [y,x]=reckon(repmat(dec,1,n0),repmat(ra,1,n0),...
      repmat(simopt.p.r(chans0)',length(dec),1),...
      repmat(simopt.p.theta(chans0)',length(dec),1)-repmat(dk,1,n0)-90);
    x=mean(x,2); y=mean(y,2);
  else
    % reduc_makepairmaps knows nothing about differential pointing, so we can use only
    % the a channels
    [y,x]=reckon(dec,ra,p.r(ch),repmat(p.theta(ch),size(dk))-dk-90);
  end

  x=cvec(x); y=cvec(y);

  % generate the standard 6 elliptical gaussian templates
  temp(:,i,1:6)=gen_templates_sub(hmap,ra,dec,dk,p,ind,ch,ignore_jitter,dp,x,y);
  
  % form the residual template if requested - with the cutdown of the
  % sky before conv2 below this does not add significantly to the run
  % time of the above
  if any(dp==9)
    % convolve the appropriate map with the appropriate rotated beam
    % this needs to follow the code in reduc_makesim subfuncs rotbeam and
    % gen_sig (under "case 'linear'")

    % take the diff beam
    db=residbeam.beammap.map(ind.a(i)).T-residbeam.beammap.map(ind.b(i)).T;
    
    % if beam maps are not available resid template will be zero
    if(~isempty(db))
      % rotate by mean dk angle - ignore any jitter!  Needs to be horizon dk
      % and not celestial, otherwise we get a 180 deg rotation
      if isfield(mapopt,'curveskyrotbeam') && mapopt.curveskyrotbeam;
        rotangle=curveskyrotbeam(mean(dec),mean(dk),p.r(ch),p.theta(ch));
      else
        rotangle=-mean(dk_hor);
      end    
      db=imrotate(db,rotangle,'bilinear','crop');

      % get the mean beam width - taking the mean is actually pointless at
      % the moment since these are nominal numbers
      fwhm=mean([p.fwhm_maj([ind.a(i),ind.b(i)]);p.fwhm_min([ind.a(i),ind.b(i)])]);
      
      % regress diffgain and pointing and remove
      %[c,t]=ffbm_deprojmap(db,residbeam.beammap.ad,fwhm/2.35,1:3);
      [c,t]=ffbm_deprojmap(db,residbeam.beammap.ad,fwhm/(sqrt(8*log(2))),2:3);
      %dbc=db-c(1).*t{1}-c(2).*t{2}-c(3).*t{3};
      dbc=db-c(1).*t{2}-c(2).*t{3}; % NB - only removing diff point
      
      % get the band index for this channel
      f=get_band_ind(p,ind,p.band(ch));
      
      if isfield(mapopt,'curveskyrescale') && mapopt.curveskyrescale
	% new better way - pull out the box around the scan
	
	% find the beam map pixel spacing and x/y span
	bmps=residbeam.beammap.ad.Field_size_deg/residbeam.beammap.ad.N_pix;
	bmsx=max(residbeam.beammap.ad.t_val_deg{1});
	bmsy=max(residbeam.beammap.ad.t_val_deg{2});
	% find the box on the sky containing the scanned trajectory
	xl=min(x(:)); xu=max(x(:));
	yl=min(y(:)); yu=max(y(:));

	% find the x direction curved sky scale factor -
	% make the pixels square at this dec angle
	xsf=1/cosd(mean(y(:)));
      
	% expand the box by the size of the beammaps
	xl=xl-bmsx*xsf; xu=xu+bmsx*xsf;
	yl=yl-bmsy;     yu=yu+bmsy;
	% make a grid of points
	m.x_tic=xl:bmps*xsf:xu;
	m.y_tic=yl:bmps:yu;
	[xx,yy]=meshgrid(m.x_tic,m.y_tic);
	% interpolate from healpix map
	tic; sm=healpix_interp(residbeam.skymap{f},xx(:),yy(:),'taylor'); toc
	% just keep the T map
	sm=reshape(sm,[size(xx),3]); sm=sm(:,:,1);
      else
	% old way - use the full sky map
	sm=residbeam.skymap(f).T;
	m=residbeam.skymap(f);
      end
	
      % convolve onto flat sky map
      mapc=conv2(sm,flipud(fliplr(dbc)),'same');
      
      % form the template
      temp(:,i,9)=interp2nan(m.x_tic,m.y_tic,mapc,x,y);

      % if wanted do it again using the uncleaned diff beam map
      if any(dp==10)
	map=conv2(sm,flipud(fliplr(db)),'same');
	temp(:,i,10)=interp2nan(m.x_tic,m.y_tic,map,x,y);
      end
      
      % diagnostic plot
      if(0)
	if(~isfield(m,'xdos'))
          m.xdos=diff(m.x_tic([1,end]))/xsf; m.ydos=diff(m.y_tic([1,end]));
	end
	setwinsize(gcf,800,1000); clf
	subplot(5,2,1); plot_healpix_map(hmap(f),1);
	hold on; plot(x,y,'k'); hold off
	colorbar; title('Healpix smoothed map from which Gauss modes derived');
	%line([xl,xu,xu,xl,xl],[yl,yl,yu,yu,yl]);
	subplot(5,2,2); plot(squeeze(temp(1:800,i,1:6)));
	title('Gaussian mode templates'); xlabel('time (samples in first 2 halfscans)');
	subplot(5,2,3); plot_map(m,sm);
	colorbar; title('Flat unsmoothed map which is convolved with residbeam');
	hold on; plot(x,y,'k');
	hold off;
	
	subplot(5,2,5);
	imagesc(residbeam.beammap.ad.t_val_deg{1},residbeam.beammap.ad.t_val_deg{2},db);
	axis xy; axis square; colorbar; title('Residbeam full');
	subplot(5,2,6);
	imagesc(residbeam.beammap.ad.t_val_deg{1},residbeam.beammap.ad.t_val_deg{2},dbc);
	axis xy; axis square; colorbar; title('Residbeam deproj');
	
	subplot(5,2,7); plot_map(m,map);
	hold on; plot(x,y,'k'); hold off
	colorbar; title('Flat map convolved w. residbeam full');
	subplot(5,2,8); plot_map(m,mapc);
	hold on; plot(x,y,'k'); hold off
	colorbar; title('Flat map convolved w. residbeam deproj');
	subplot(5,2,9); plot(squeeze(temp(1:800,i,10)));
	title('Residbeam dp0000 template'); xlabel('time (samples in first 2 halfscans)');
	subplot(5,2,10); plot(squeeze(temp(1:800,i,9)));
	title('Residbeam dp1100 template'); xlabel('time (samples in first 2 halfscans)');
	keyboard
      end
      
    end          
  end
end

if docrosstalk
  % Do crosstalk templates
  if size(temp,3)<7
    temp=cat(3,temp,zeros(size(temp,1),size(temp,2),2));
  end
  np=numel(find(p.mce_col(ind.a)==1)); % number of pairs per column
  for i=1:n
    % Each channel gets two additional templates: the T template from its upstream mux
    % neighbor and the T template from its downstream mux neighbor. Since the templates
    % are ordered by pair number, and pair number is monotonic with mux_row number, this
    % is just the adjacent templates. 
    downind=mod(i,16) + 1 + np*floor((i-1)/np);
    upind=mod(i-2,16) + 1 + np*floor((i-1)/np);
    
    temp(:,i,7)=temp(:,upind,1);
    temp(:,i,8)=temp(:,downind,1);
  end
end

% Output the templates in the requested order
for k=1:numel(dp)
  T(mapind,:,k)=temp(:,:,dp(k));
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function T=gen_templates_sub(map,ra,dec,dk,p,ind,ch,ignore_jitter,dp,x,y)

% Number of templates is always the 6 standard Gaussian modes.
% Modes not requested are returned as zeros
n=6;
T = zeros(length(dec),n);

% Polynomial templates, currently unused
if(0)
  azz=d.azoff;
  % Construct poly templates:
  for j=0:porder
    T(:,j+1)=azz.^j;
  end
end

% pick the right frequency index for this channel,
% this will deliberatly fail if there are more frequencies
% than deprojtion maps that have been handed in.
f = get_band_ind(p,ind,p.band(ch));
hmap = map(f);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Template 1: temperature map for monopole
T(:,1)=healpix_interp(hmap,x,y,'taylorT');

% Return if we don't want more templates
if max(dp(dp<n))==1
  return
end

% keep theta and phi around for later use
theta=cvec((90-y)*pi/180);
phi=cvec(x*pi/180);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Templates 2-3: differential pointing

% Trick healpix_interp into moving the 1st derivs around on the sky using the second
% derivs
dmap=hmap;
dmap=rmfield(dmap,'map');
dmap.map(:,1)=hmap.map(:,4); % dTdth -> T
dmap.map(:,4)=hmap.map(:,10); % dTdth2 -> dTdth
dmap.map(:,7)=hmap.map(:,13); %  dTdthdph -> dTdph
dTdth=healpix_interp(dmap,x,y,'taylor1T');

dmap.map(:,1)=hmap.map(:,7); % dTdph -> T
dmap.map(:,7)=hmap.map(:,16); % dTdph2 -> dTdph
dmap.map(:,4)=hmap.map(:,13); %  dTdphdth -> dTdth
dTdph=healpix_interp(dmap,x,y,'taylor1T');

% dT/dphi in a synfast map is really dT/dphi/sin(theta), so we'll multiply by
% sin(theta) up front to put it in terms of coordinate phi
dTdph=dTdph.*sin(theta);

% Using the chain rule, we'll construct the templates using d^2T/dx^2 and
% d^T/dy^2 for beamwidth and ellipticity. We'll need the second derivative maps, which
% will have to be nearest 
% neighbor interpolated for now. Trick healpix_interp into returning the two second
% derivatives and one cross derivative.
dmap=hmap;
dmap=rmfield(dmap,'map');
dmap.map(:,1)=hmap.map(:,10); % dTdth2 -> T
dmap.map(:,2)=hmap.map(:,13); % dTdthdph -> Q
dmap.map(:,3)=hmap.map(:,16); %  dTdph2 -> U
derivs=healpix_interp(dmap,x,y,'healpixnearest');

clear dmap

% Again, put phi derivatives in terms of coordinate phi
d2Tdth2 = derivs(:,1);
d2Tdthdph = derivs(:,2).*sin(theta);
d2Tdph2 = derivs(:,3).*sin(theta).^2;

clear derivs
clear x y

% The derivatives of the healpix coordinate system (theta,phi) w.r.t. the focal plane
% coordinate system (x,y) only depends dec and dk, which are supposed to be
% constant. If the jitter in dec and dk are small enough, we can use their mean value
% and massively speed up the calculation of the template. This has the effect of making
% the derivatives dtheta/dx, dphi/dx, etc., scalars instead of vectors.
if ignore_jitter
  theta=theta(1);
  dec=dec(1);
  dk=dk(1);
  phi=phi(1);
  ra=ra(1);
end

% Derivatives expressed in the healpix (theta,phi) coordinate system need to be
% expressed in the focal plane (x,y) coordinate system. This requires the chain
% rule. The parital derivatives of theta(x,y) and phi(x,y) are needed and could be
% calculated analytically, but we will determine them numerically by displacing x,y by
% a small amount and seeing how theta and phi vary. Here we displace by a small amount
% and then later scale the derivatives up to one degree.

% When calculating the numerical derivatives of theta,phi w.r.t. x,y, offset the beam
% by this many degrees
du=1e-3;

% First displace the beam centroids by +x, calculate dtheta and dphi
pd=displace_beamcen(p,du,0,ch);
[dthdx_1,dphdx_1,theta_plusx,phi_plusx]=calc_num_deriv(pd,theta,phi,dk,ra,dec,ch);

% Now displace the beam centroids by -x, calculate dtheta and dphi
pd=displace_beamcen(p,-du,0,ch);
[dthdx_2,dphdx_2,theta_minx,phi_minx]=calc_num_deriv(pd,theta,phi,dk,ra,dec,ch);
dthdx_2=-dthdx_2;
dphdx_2=-dphdx_2;

% Take the average of the derivative to the left and deriv to the right
dthdx=(dthdx_1+dthdx_2)/2;
dphdx=(dphdx_1+dphdx_2)/2;

% Now do the same for y
pd=displace_beamcen(p,0,du,ch);
[dthdy_1,dphdy_1,theta_plusy,phi_plusy]=calc_num_deriv(pd,theta,phi,dk,ra,dec,ch);

pd=displace_beamcen(p,0,-du,ch);
[dthdy_2,dphdy_2,theta_miny,phi_miny]=calc_num_deriv(pd,theta,phi,dk,ra,dec,ch);
dthdy_2=-dthdy_2;
dphdy_2=-dphdy_2;

% Take the average of the derivative above and below
dthdy=(dthdy_1+dthdy_2)/2;
dphdy=(dphdy_1+dphdy_2)/2;

% Apply chain rule for first derivatives in two dimensions
% dT/dx
T(:,2) = dTdth.*dthdx + dTdph.*dphdx;

%dT/dy
T(:,3) = dTdth.*dthdy + dTdph.*dphdy;

% Return if we don't want more templates
if max(dp(dp<n))<=3
  T=finalizeT(T,du,n);
  return
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Template 4: beamwidth

% The fit coefficient for beamwidth shouldn't depend on focal plane orientation, so we
% really don't need to convert to d^2T/dx^2 and d^2/dy^2, but we'll do it anyways
% because we need it for ellipticity below

d2thdx2 = dthdx_1-dthdx_2;
d2phdx2 = dphdx_1-dphdx_2;

d2thdy2 = dthdy_1-dthdy_2;
d2phdy2 = dphdy_1-dphdy_2;

% For cross derivatives (needed for ellipticity) we first calculate dtheta/dx and
% dphi/dx with zero y displacement (already done above) and then calculate dtheta/dx
% with with a non zero y-displacement (done below)

% offset trajectory at +y displacement
theta0=theta_plusy;
phi0=phi_plusy;

% dtheta/dx and dphi/dx at +y displacement
pd=displace_beamcen(p,du,du,ch);
[dthdx_1,dphdx_1]=calc_num_deriv(pd,theta0,phi0,dk,ra,dec,ch);

pd=displace_beamcen(p,-du,du,ch);
[dthdx_2,dphdx_2]=calc_num_deriv(pd,theta0,phi0,dk,ra,dec,ch);
dthdx_2=-dthdx_2;
dphdx_2=-dphdx_2;

dthdx_plusy=(dthdx_1+dthdx_2)/2;
dphdx_plusy=(dphdx_1+dphdx_2)/2;

% dtheta/dx and dphi/dx at -y displacement
theta0=theta_miny;
phi0=phi_miny;

pd=displace_beamcen(p,du,-du,ch);
[dthdx_1,dphdx_1]=calc_num_deriv(pd,theta0,phi0,dk,ra,dec,ch);

pd=displace_beamcen(p,-du,-du,ch);
[dthdx_2,dphdx_2]=calc_num_deriv(pd,theta0,phi0,dk,ra,dec,ch);
dthdx_2=-dthdx_2;
dphdx_2=-dphdx_2;

dthdx_miny=(dthdx_1+dthdx_2)/(2);
dphdx_miny=(dphdx_1+dphdx_2)/(2);

% How much did dtheta/dx and dphi/dx change from -y to +y displacement? This is the
% cross derivative. 
d2thdxdy = (dthdx_plusy - dthdx_miny)/2;
d2phdxdy = (dphdx_plusy - dphdx_miny)/2;

clear dthdx_1 dthdx_2 dphdx_1 dphdx_2 dthdx_plusy dphdx_plusy dthdx_miny dphdx_miny ...
    theta0 phi0 theta_plusy theta_miny

if(0)
  % Below we can calculate dth2dydx, but it is the the same as dth2dxdy, so we
  % don't. The code is left so we can verify this statement in the future.
  
  % offset trajectory at +x displacement
  theta0=theta_plusx;
  phi0=phi_plusx;

  % dtheta/dy and dphi/dy at +x displacement
  pd=displace_beamcen(p,du/2,du,ch);
  [dthdy_1,dphdy_1]=calc_num_deriv(pd,theta0,phi0,dk,ra,dec,ch);
  
  pd=displace_beamcen(p,du/2,-du,ch);
  [dthdy_2,dphdy_2]=calc_num_deriv(pd,theta0,phi0,dk,ra,dec,ch);
  dthdy_2=-dthdy_2;
  dphdy_2=-dphdy_2;
  
  dthdy_plusx=(dthdy_1+dthdy_2)/2;
  dphdy_plusx=(dphdy_1+dphdy_2)/2;
  
  % dtheta/dy and dphi/dy at -x displacement
  theta0=theta_minx;
  phi0=phi_minx;
  
  pd=displace_beamcen(p,-du/2,du,ch);
  [dthdy_1,dphdy_1]=calc_num_deriv(pd,theta0,phi0,dk,ra,dec,ch);
  
  pd=displace_beamcen(p,-du/2,-du,ch);
  [dthdy_2,dphdy_2]=calc_num_deriv(pd,theta0,phi0,dk,ra,dec,ch);
  dthdy_2=-dthdy_2;
  dphdy_2=-dphdy_2;
  
  dthdy_minx=(dthdy_1+dthdy_2)/2;
  dphdy_minx=(dphdy_1+dphdy_2)/2;
  
  % How much did dtheta/dy and dphi/dy change from -x to +x displacement? This is the
  % cross derivative
  d2thdydx = dthdy_plusx - dthdy_minx;
  d2phdydx = dphdy_plusx - dphdy_minx;
  
  clear dthdy_1 dthdy_2 dphdy_1 dphdy_2 dthdy_plusx dphdy_plusx dthdy_minx dphdy_minx ...
      theta0 phi0

end


d2Tdx2 = dTdth.*d2thdx2 + dTdph.*d2phdx2 + ...
         d2Tdth2.*(dthdx).^2 + d2Tdph2.*(dphdx).^2 + 2*d2Tdthdph.*dthdx.*dphdx;

d2Tdy2 = dTdth.*d2thdy2 + dTdph.*d2phdy2 + ...
         d2Tdth2.*(dthdy).^2 + d2Tdph2.*(dphdy).^2 + 2*d2Tdthdph.*dthdy.*dphdy;

d2Tdxdy = dTdth.*d2thdxdy + dTdph.*d2phdxdy + ...
          d2Tdth2.*dthdx.*dthdy + d2Tdph2.*dphdx.*dphdy + ...
          d2Tdthdph.*(dthdx.*dphdy + dphdx.*dthdy);

% Differential beamwidth
T(:,4) = d2Tdx2 + d2Tdy2;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Templates 5-6: ellipticity

% plus - d^2T/dx^2 - d^2T/dy^2
T(:,5) = d2Tdx2 - d2Tdy2;

% cross - 2(d^2T/dxdy)
T(:,6) = 2*d2Tdxdy;

%%%%%%%%%%%%%%%%%%%%%%%
% Scale to 1 degree offsets and reshape
T=finalizeT(T,du,n);

return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function T=finalizeT(T,du,n)

% Scale to 1 degree offsets

% First derivatives
T(:,2:3)=T(:,2:3)/du;

% Second derivatives
T(:,4:6)=T(:,4:6)/du^2;

% Reshape to final dims
T = reshape(T,size(T,1),1,n);

return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [dthdu,dphdu,theta_out,phi_out]=calc_num_deriv(pd,theta,phi,dk,ra,dec,ch)

% find ra/dec trajectory for the displaced channels
[y,x]=reckon(dec,ra,pd.r(ch),repmat(pd.theta(ch),size(dk))-dk-90);

theta_out = cvec((90-y)*pi/180);
phi_out = cvec(x*pi/180);

dthdu=theta_out - theta;
dphdu=phi_out - phi;

return
