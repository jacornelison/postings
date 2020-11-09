function map=inject_pntsrc(m,map,ra,dec,flux,polfrac,polang)
% map=inject_pntsrc(m,map,ps)
%
% Inject point sources into maps according to input source list
%
% flux should be source fluxes in Jy
% map should be abs cal'd to uk
%
% e.g
% maps=inject_pntsrc(m,map,ra,dec,s,p)

% if pol frac scalar assume constant
if(isscalar(polfrac))
  polfrac=polfrac*ones(size(flux));
end

% if pol angle empty
if(~exist('polang'))
  polang=rand(size(flux))*pi;
end

% wmap lists unmeasured as negative flux
flux(flux<0)=0;

% bicep field wraps around RA=0
ra(ra>180)=ra(ra>180)-360;

% beam width in sigmas
beamwid=(31/60)/2.35;

% freq in Hz
nu=150e9;

% calc scale factor from Jy to uK peak in the map
jy2ukpeak=1./ukpeak_to_janskies(nu,beamwid*pi/180);

% only need to make lump out to 4 beam widths in pixels
n=ceil(beamwid*4/m.pixsize);

% for each source
for k=1:length(flux)
  
  x=interp1(m.x_tic,[1:m.nx],ra(k));
  y=interp1(m.y_tic,[1:m.ny],dec(k));
  
  xr=round(x); yr=round(y);
  
  if(isfinite(x)&isfinite(y))
    
      % find the map pixel boundaries for the source
      % making sure we don't go over the bounds
      lx=max([xr-n,1]); ux=min([xr+n,m.nx]);
      ly=max([yr-n,1]); uy=min([yr+n,m.ny]);
      
      % make a coordinate grid in degrees on sky
      % lx,ux,ly,uy are already limited to map bounds
      x=[lx-x:ux-x]*m.pixsize;
      y=[ly-y:uy-y]*m.pixsize;
      [x,y]=meshgrid(x,y);
      
      % generate gaussian lump, offset from nominal pixel by pixel offset
      g=egauss2([1,0,0,beamwid,beamwid,0],x,y);
      
      % scale gaussian blob such that 1Jy source has correct uK peak brightness
      g=g*jy2ukpeak;
      
      % inject the gaussian blob into T map
      map.T(ly:uy,lx:ux)=map.T(ly:uy,lx:ux)+g*flux(k);
      
      % inject the blob into Q/U maps
      % using IAU pol angle convention
      map.Q(ly:uy,lx:ux)=map.Q(ly:uy,lx:ux)+g*flux(k)*polfrac(k)*cos(2*polang(k));
      map.U(ly:uy,lx:ux)=map.U(ly:uy,lx:ux)+g*flux(k)*polfrac(k)*sin(2*polang(k));
  end
end

return
