function ffbm_makesetofegauss(beamopt,filename)
%function ffbm_makesetofegauss(beamopt,filename)
%
% beamopt structure has fields:
%
% mapsize: size of postage stamp map in degrees (default 8)
% stepsize: resolution of map in degrees (default 0.1)
% small: if true, use a 2x2 instead of 8x8 default map size
%
% epoch: a tag or date to control the epoch used in get_array_info
% experiment: almost certainly unused; might be used in 'modelcomposite'?
%
% Set one of the following to "true" to run some specific scenarios:
%   bigellip: multiply ellipticities by 5x
%   circ: ideal circular gaussians, 0.213 deg.  (should be redundant with 'ideal')
%   halfdiff: divide diffpoint by 2
%
% tqu: 1 to make T, Q, U beams; 0 to make a single beam
%
% signal: include mapped beams (0 for noise only)
% noise: noise level to be added to constructed maps
%        = 'none' or 'signal' noiseless constructed map
%        = 'uber' map with uberchopper noise level (default)
%        = 'zechoppa' map with ze choppa noise level
%        = 'midnoise' map with noise level between ze choppa and uber chopper
%        = 'modelcomposite' used for keck, the noise levels are
%        estimated from estimatenoise.m
%        = 'noiseonlyb2' map with only white noise at level of
%        bicep2 rev2 composite beam maps, stddev=7.3e-6;
%        = 'ubernoiseonly' map with only uberchopper noise levels, 
%        stddev=2.8e-5
%
% The following fields control the beam properties and
% act as in reduc_makesim.  Little buddies, xtalk not implemented.
%
% beamopt.beamcen  - determines which beam centers to get back from get_array_info
%                   = 'ideal' - feed offset angles design vals from fp_data (default)
%                   = 'obs'   - as observed and recorded in the beams files
%                   = 'zero'  - zero p.r and p.theta
% beamopt.beamwid   = 'zero'  - no map smoothing
%                   = 'ideal' - as returned by fp_data files
%                   = 'obs'   - read from aux_data file (default)
% beamopt.diffpoint = 'ideal' - no differential pointing (default)
%                   = 'obs'  - include differential pointing
% beamopt.chi       = 'ideal' - chi=design vals wrt theta. (except CP wrt theta=0) (default)
%                   = 'obs'   - chi as returned by get_array_info
% beamopt.epsilon   = 'ideal' - epsilon=0
%                   = 'obs'   - epsilon as returned by get_array_info
% beamopt.polofs    = 'ideal' - no "Roger effect" pol-dependent beams
%                   = 'obs'   - include "Roger effect" pol-dependent beams
% beamopt.lb        - not implemented
%                   = false (default) - do nothing
%                   = true - add 180 degrees to p.theta to simulate a little buddy. All
%                            dithered parameters are interpreted as occuring after the
%                            offset. This only occurs in reduc_makesim. reduc_makepairmaps
%                            will use the non-offset p.r and p.theta, and so accumulate
%                            the little buddy timestream assuming the main beam
%                            trajectory (unless mapopt.beamcen='assim', etc.). It is
%                            possibly useful to specify beamopt.abgain=ones(Ndet,Nreal)*x
%                            where x is the little buddy amplitude. 
%
% beamopt.xtalk  - not implemented
%                = false (default) do nothing
%                = level of inductive cross talk (i.e. .05 for 5 percent cross talk)
%
% beamopt.xtalk_relgain - not implemented
%                       = false to calculate xtalk in units of input map signal
%                       = true (default) to adjust by ADU/airmass scaling
%
%
% Example usage:
%   beamopt.beamwid='aux_data/beams/beams_bicep2_obs_rwa_20130607.csv';
%   beamopt.diffpoint='ideal';
%   beamopt.beamcen='ideal';
%   beamopt.epsilon='obs';
%   beamopt.mapsize=6;
%   beamopt.noise='none';
%   ffbm_makesetofegauss(beamopt,'bicep2_6deg');
%
%CLW 20140525
%RWO 20140813

if ~exist('beamopt','var')
  beamopt=[];
end
if ~exist('filename','var')
  filename=[];
end

if ~isfield(beamopt,'mapsize')
  beamopt.mapsize=[];
end
if ~isfield(beamopt,'stepsize')
  beamopt.stepsize=[];
end
if ~isfield(beamopt,'signal')
  beamopt.signal=[];
end
if ~isfield(beamopt,'epoch')
  beamopt.epoch=[];
elseif isnumeric(beamopt.epoch)
  beamopt.epoch=num2str(beamopt.epoch);
end
if ~isfield(beamopt,'experiment')
  beamopt.experiment=[];
end
if ~isfield(beamopt,'bigellip')
  beamopt.bigellip=[];
end
if ~isfield(beamopt,'halfdiff')
  beamopt.halfdiff=[];
end
if ~isfield(beamopt,'circ')
  beamopt.circ=[];
end
if ~isfield(beamopt,'small')
  beamopt.small=[];
end
if ~isfield(beamopt,'noise')
  beamopt.noise='uber';
end
if ~isfield(beamopt,'tqu')
  beamopt.tqu=0;
end

if isempty(beamopt.small)
  beamopt.small=false;
end
if isempty(beamopt.mapsize)
  beamopt.mapsize=8;
  if beamopt.small
    beamopt.mapsize=2;
  end
end
if isempty(beamopt.experiment)
  beamopt.experiment=get_experiment_name();
end
if isempty(beamopt.stepsize)
  beamopt.stepsize=0.1;
  if beamopt.small
    beamopt.stepsize=0.0625;
  end
end
if isempty(beamopt.bigellip)
  beamopt.bigellip=false;
end
if isempty(beamopt.halfdiff)
  beamopt.halfdiff=false;
end
if isempty(beamopt.circ)
  beamopt.circ=false;
end
if isempty(beamopt.signal)
  beamopt.signal=1;
end
if isempty(beamopt.noise)
  beamopt.noise='none';
end
if isempty(beamopt.tqu)
  beamopt.tqu=1;
end

% Make sure default beam props work exactly as in reduc_makesim
tmpopt=beamopt;
tmpopt.coord='C';
tmpopt=get_default_simopt(tmpopt);
flist={'beamcen','beamwid','diffpoint','chi','epsilon','lb','xtalk','polofs'};
for i=1:length(flist)
  beamopt.(flist{i})=tmpopt.(flist{i});
end

% Get array info and index arrays
if strcmp(beamopt.beamcen,'zero')
  beamcen=[];
else
  beamcen=beamopt.beamcen;
end
[p,ind]=get_array_info(beamopt.epoch,beamcen,beamopt.chi,beamopt.beamwid,beamopt.diffpoint,[],beamopt.epsilon,[],beamopt.polofs);
if beamopt.bigellip
  [tmps,tmpc,tmpp]=egauss2_mmt2scp(p.fwhm_maj,p.fwhm_min,p.alpha);
  tmpp=tmpp*5;
  tmpc=tmpc*5;
  [p.fwhm_maj,p.fwhm_min,p.alpha]=egauss2_scp2mmt(tmps,tmpc,tmpp);
end
p.sigma_maj=p.fwhm_maj./(2*sqrt(2*log(2)));
p.sigma_min=p.fwhm_min./(2*sqrt(2*log(2)));
% p.alpha is calculated from the p.theta line
% but p.alpha is calculated before p.theta is rotated
% by the drumangle
% but that's right, because if you rotate p.theta
% by the drumangle, you don't need to rotate p.alpha
% since p.alpha is calculated from the p.theta line
% so p.alpha2 from the x axis is as follows.
% Note we need pair center theta, not per-det theta with diffpoint included!
[p0,ind0]=get_array_info(beamopt.epoch,beamcen,beamopt.chi,beamopt.beamwid,'ideal',[],beamopt.epsilon,[],beamopt.polofs);
p.alpha2=p.alpha+p0.theta;

if (beamopt.circ)
  % Create circular gaussians, diff point included
  p.p=zeros(528,1);
  p.c=zeros(528,1);
  p.sigma=ones(528,1)*0.213;
  [fwhm_maj,fwhm_min,alpha2]=egauss2_scp2mmt(p.sigma,p.c,p.p);
  p.sigma_maj=fwhm_maj./(2*sqrt(2*log(2)));
  p.sigma_min=fwhm_min./(2*sqrt(2*log(2)));
  p.alpha2=alpha2;
  p.polofs_x=zeros(528,1);
  p.polofs_y=zeros(528,1);
end

% Recenter the maps.

% Use (new) standard pipeline helper functions.
% This correctly uses pair-centered coordinate system x', y'.
[p.rcent,p.thcent,xp,yp]=paircenter(p.r(ind.a),p.theta(ind.a),p.r(ind.b),p.theta(ind.b));
p.xp=0*p.r;
p.xp(ind.a)=xp{1};
p.yp(ind.a)=yp{1};
p.xp(ind.b)=xp{2};
p.yp(ind.b)=yp{2};

% Create x' and y' bins

% Note have to flip the y' direction
x_bin=[-beamopt.mapsize/2:beamopt.stepsize:beamopt.mapsize/2];
y_bin=[-beamopt.mapsize/2:beamopt.stepsize:beamopt.mapsize/2];
y_bin=fliplr(y_bin); % so that I can imagesc on sky at dk0

[X,Y]=meshgrid(x_bin,y_bin);
% Since we've called fliplr, the coordinates X and Y are now identical
% to x' and y' as defined in the coordinate system definition document
% at http://bicep0.caltech.edu/~spuder/analysis_logbook/analysis/20131005_sidelobe_coordinates/sidelobe_coordinates.pdf
% These coordinates are centered on the *pair* center and oriented so that
% y' is approximately oriented toward the zenith (south cel pole) at dk000.

%%%%%%%%%%%%%%%%%%%%%
% generate the ad structure with calc_ad2
% ad=calc_ad2(Field_size_deg,N_pix);

sizeindeg=length(x_bin)*beamopt.stepsize;
Field_size_deg=[sizeindeg sizeindeg];

N_pix=[length(x_bin) length(y_bin)];
ad=calc_ad2(Field_size_deg,N_pix);

%%%%%%%%%%%%%%%%%%%%%
% Construct the elliptical Gaussian beams for each channel

% Empty map structure -- now making a structure similar to pairmaps
% or coadded maps, not a single big matrix
map=repmat(struct(),length(p.gcp),1);

for ii=1:length(ind.a)
  for idx=[ind.a(ii) ind.b(ii)]    % A, B
    sigma_maj=p.sigma_maj(idx);
    sigma_min=p.sigma_min(idx);
    alpha=p.alpha2(idx);
    alpha=alpha*pi/180;
    % Note xcent, ycent are stored in B det only
    % x=p.x(idx)-p.xcent(ind.b(ii));
    % y=p.y(idx)-p.ycent(ind.b(ii));
    x=p.xp(idx);
    y=p.yp(idx);
    if beamopt.halfdiff
      x=x/2;
      y=y/2;
    end  
    beamparam=[1 x y sigma_maj sigma_min alpha];
   
    beammap = egauss2(beamparam,X,Y);
    map(idx).T=beammap;
    if beamopt.tqu
      % rotate for Roger effect
      rot=X*p.polofs_x(idx)+Y*p.polofs_y(idx);
      qbeam=beammap.*cosd(2*rot)*(1-p.epsilon(idx))/(1+p.epsilon(idx));
      ubeam=beammap.*sind(2*rot)*(1-p.epsilon(idx))/(1+p.epsilon(idx));
      % enforce rule that mean(umap)=0
      qmean=mean(qbeam(:));
      umean=mean(ubeam(:));
      th=atan2(umean,qmean);
      map(idx).Q=cos(th)*qbeam+sin(th)*ubeam;
      map(idx).U=cos(th)*ubeam-sin(th)*qbeam;
    end
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%
% nan out the dark squids and dark TESs
% can also nan out dets without well measured beams, if desired
map=nandarks(map,p,ind);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% inject noise
% the noise numbers here were taken from
% 1 map of 1 pixel using a low noise region
% detailed here:
% http://bmode.caltech.edu/~spuder/analysis_logbook/analysis/20130610_elipgausmaps/
% the midnoise level number is just a made up number
switch (beamopt.noise)
  case 'zechoppa'
    % ze choppa
    %noise in sky 
    %mu=3.0734e-05;
    %stddev=0.0032;
    mu=4.2e-05;
    stddev=0.0044;
  case 'uber'
    %for uberchopper
    mu=-2.4412e-06
    stddev=8.6043e-4
  case 'midnoise'
    mu=0;
    stddev=0.0022;
  case 'modelcomposite'
    if strcmp(beamopt.experiment,'keck')
      %the numbers here are the median of the noise
      %from estimatenoise.m using only rgls
      %{
      %script here:
      for ii=0:4
        mapfilename=['maps/rx' num2str(ii) '_all'];
        [noise(ii+1,:) mn(ii+1,:)]=estimatenoise(mapfilename,'keck',2012,ii,0);
        nn(ii+1)=nanmedian(noise(ii+1,:));
        mmn(ii+1)=nanmedian(mn(ii+1,:));
      end
      nanmedian(nn)
      %}
      mu=0;
      stddev=0.0010;
    else
      error(['Noise "modelcomposite" supported only for Keck Array.  Sort of.']);
    end
  case 'noiseonlyb2'
    mu=0;
    stddev=7.3e-6;
    beamopt.signal=false;
  case 'ubernoiseonly'
    mu=0;
    stddev=2.8e-5;
    beamopt.signal=false;
  case {'none','signal'}
    mu=0;
    stddev=0;
  otherwise
    error(['Unrecognized noise specification ' beamopt.noise]);
end  

%generate gaussian noise
noisemap=normrnd(mu,stddev,[length(y_bin),length(x_bin),length(p.gcp)]);
for ii=1:length(p.gcp)
  if isempty(map(ii).T)
    continue
  end
  if beamopt.signal
    map(ii).T=map(ii).T+noisemap(:,:,ii);
  else
    map(ii).T=noisemap(:,:,ii);
  end
end
if beamopt.tqu
  qnoisemap=normrnd(mu,stddev,[length(y_bin),length(x_bin),length(p.gcp)]);
  unoisemap=normrnd(mu,stddev,[length(y_bin),length(x_bin),length(p.gcp)]);
  for ii=1:length(p.gcp)
    if isempty(map(ii).Q)
      continue
    end
    if beamopt.signal
      map(ii).Q=map(ii).Q+qnoisemap(:,:,ii);
      map(ii).U=map(ii).U+unoisemap(:,:,ii);
    else
      map(ii).Q=qnoisemap(:,:,ii);
      map(ii).U=unoisemap(:,:,ii);
    end
  end
end

%%%%%%%%%%%%%%%%%%%%%
% renormalize map if signal is present
if beamopt.signal
  map=renormmap(map);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Save file

% The logic around filenames here needs to be cleaned up.
% Could either just use the supplied filename in all cases,
% or else construct one with various flags like in pairmaps & maps.

% If we're given a filename, just use it
if ~isempty(filename)
 fullfname=fullfile('beammaps/maps_constructed/',filename);
elseif strcmp(mapname,'signal')
  fullfname=['beammaps/maps_constructed/egaussmaps_signal'];
elseif beamopt.circ
  fullfname=['beammaps/maps_constructed/circgaussmaps_signal'];
elseif beamopt.halfdiff
  fullfname=['beammaps/maps_constructed/egaussmaps_halfdiff'];
elseif strcmp(beamopt.noise,'noiseonlyb2')
  fullfname=['beammaps/maps_constructed/B2_noiseonlymap'];
elseif strcmp(beamopt.noise,'ubernoiseonly')
  fullfname=['beammaps/maps_constructed/ubernoiseonlymap'];
else
  error(['Don''t know where to save the output.']);
end

% save
beamopt.p=p;
save(fullfname,'map','ad','beamopt');

return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function map=renormmap(map)
flist={'T','Q','U'};
for ii=1:size(map,1)
  % Norm Q and U (if present) from T
  mapnorm=sum(sum(map(ii).T));
  for jj=1:length(flist)
    if ~isfield(map(ii),flist{jj})
      continue
    end
    map(ii).(flist{jj})=map(ii).(flist{jj})/mapnorm;
  end
end
return

%%%%%%%%%%%%%%%%%%%%%%%
function map=nandarks(map,p,ind);
cdark=false(size(ind.e));
cdark(ind.o)=true;
cdark(ind.d)=true;

% Possibly NaN out the pixels that aren't well measured and have
% ideal beam params.  Currently not doing this.
% cnotmeas=false(size(p.r));
% cnotmeas(ind.a)=cnotmeas(ind.a) | (p.r(ind.a)==p.r(ind.b));
% cnotmeas(ind.b)=cnotmeas(ind.b) | cnotmeas(ind.a);
% cnotmeas(ind.l) = cnotmeas(ind.l) | (p.p(ind.l)==0 & p.c(ind.l)==0);
% pixelsmarkedasbad=find(cnotmeas);
% cdark = cdark | cnotmeas

for ii=find(cdark)
  map(ii).T=NaN*map(ii).T;
  if isfield(map(ii),'Q')
    map(ii).Q=NaN*map(ii).Q;
    map(ii).U=NaN*map(ii).U;
  end
end

return
