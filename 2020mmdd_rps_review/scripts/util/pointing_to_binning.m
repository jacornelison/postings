function point=pointing_to_binning(d,mapind,m,p,deta,detb)
% point=pointing_to_binning(d,mapind,m,p,deta,detb)
%
% Generates map pixel trajectory information given telescope motion and
% detector pointing.
%
% INPUTS
%   d         TOD date structure.
%
%   mapind    Boolean mask of time samples within the TOD which correspond to
%             half-scans to bin, as from make_mapind().
%
%   m         Map definition structure, as from get_map_defn().
%
%   p         Detector metadata structure, as from get_array_info().
%
%   deta      Index into fields of p for the A/pair-sum detector of interest.
%
%   detb      Index into fields of p for the B/pair-diff detector of interest.
%
% OUTPUTS
%   point     Structure containing map-pointing information:
%
%               .x0, .y0    Raw detector pointing on the map grid. (Different
%                           than .x,.y for HEALPix binning, and required for
%                           making Taylor-interpolating pointing operators.)
%
%               .nx, .ny    Number of histogram bins in the x,y directions.
%
%               .lx, .hx,   Low and high bounds of the histogram in each of
%               .ly, .hy    x,y directions.
%
%               .x, .y      Detector pointing on the map grid. In the case of
%                           binning to HEALPix maps, these will be the
%                           ring-ordered pixel numbers to bin to.
%
%               .c, .s      Corresponding detector cosines/sines for
%                           polarized signals.
%
% EXAMPLE
%
%   [tod,mapopt] = prepare_pairmap_tod(tag, mapopt);
%   liftvars(tod); clear tod
%   mapind = make_mapind(d, fs);
%   m = get_map_defn('bicep');
%   for ii=1:length(ind.rgla)
%     point = pointing_to_binning(d, mapind, m, p, ind.rgla(ii), ind.rglb(ii));
%     ...
%   end
%

% default assumption is to bin into RA/Dec maps. reduc_makepairmaps will
% normally set this explicitly with the value taken from mapopt.binning.
if ~isfield(m, 'binning') || isempty(m.binning)
  m.binning = 'ra';
end

% xaxis of map can be ra or az
switch m.binning
 case 'az'
  % Note that the simple replacement of ra with az makes the
  % detector reckoning below only approximately correct.  For
  % a telescope not at the South Pole it would be very wrong.
  x0=d.pointing.hor.az(mapind);
 case 'ra'
  x0=d.pointing.cel.ra(mapind);
end
y0=d.pointing.cel.dec(mapind);
dk=d.pointing.cel.dk(mapind);

% find ra/dec trajectory for this pair
[y,x]=reckon(y0,x0,p.r(deta),p.theta(deta)-90-dk);
if ~strncmpi(m.type,'healpix',7)
  if(m.lx>0)
    x(x<0)=x(x<0)+360;
  end
end

% alternate rectangular projections (rarely used)
switch m.proj
 case 'ortho'
  [x,y]=radec_to_sinproj(x,y,m.racen,m.deccen);
 case 'tan'
  [x,y]=radec_to_tanproj(x,y,m.racen,m.deccen);
 case 'arc'
  [x,y]=radec_to_arcproj(x,y,m.racen,m.deccen);
 case 'radec'
  % do nothing
end

% To support Taylor interpolation in the matrix, forward along the
% unmodified coordinates. (Healpix effectively bins in the next step.)
point.x0 = x;
point.y0 = y;

% Setup for binning
if ~strncmpi(m.type,'healpix',7)
  point.nx=m.nx; point.lx=m.lx; point.hx=m.hx;
  point.ny=m.ny; point.ly=m.ly; point.hy=m.hy;
  point.x=x; point.y=y;
else
  % determine mapping from x,y to healpix index
  point.pix=ang2pix(m.nside,x,y,m.proj);
  % find the list of pixels with hits with indices into upix array
  [point.upix,dum,point.x]=unique(point.pix);
  point.y=zeros(size(point.x));
  % ensure non-zero nx to avoid crash
  point.nx=max([1,length(point.upix)]);
  point.lx=0.5; point.hx=point.nx+0.5;
  point.ny=1; point.ly=-0.5; point.hy=0.5;
end

% calc the pol angle for the pair
th=p.theta(detb)-dk;
thref=p.chi_thetaref(detb)-dk;
alpha=chi2alpha(x0,y0,p.r(detb),th,p.chi_mean(detb),thref);

% Ganga alpha_h+beta_h
ap=(alpha+p.beta(detb))*pi/180;

% take sin and cos
point.c=cos(2*ap); point.s=sin(2*ap);

% attach the relevant map description structure for context
point.m=m;

end
