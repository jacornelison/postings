function m=get_map_defn(type,resmult,proj)
% m=get_map_defn(type,resmult,proj)
%
% Fetch the map definition structure - used by reduc_makesim and
% reduc_makecomap
%
% type = map type
% resmult = resolution multiplier
% Positive odd integers will result in resolution supersets
% For negative integers the first row/col stay in the same
% position on the sky N additional row/cols are added between
% subsequent ones as occurs when padding an FT
% (See bkspt posting
% http://bicep.rc.fas.harvard.edu/bkspt/analysis/20160713_comb_bk_and_spt)

if(~exist('resmult','var'))
  resmult=[];
end
if(~exist('proj','var'))
  proj=[];
end

if(isempty(resmult))
  resmult=1;
end
if(isempty(proj))
  if(isempty(strfind(type,'healpix')))
    proj='radec';
  else
    proj='ring';
  end
end
  
if(resmult==0)
  % doesn't make sense and may happen due to old code
  resmult=1;
  warning('fed resmult=0 - forcing to 1 - fix the code');
end

% store these for ref
m.proj=proj;
m.type=type;
m.resmult=resmult;

switch lower(deblank(type))
  
  case {'healpix128','healpix256','healpix512'}
    % for healpix the only spec is the nside
    m.nside=str2double(type(8:end));
    return
   
  case 'quad05'

    % There is no point in trying to have each scan fall into a
    % specific row of map pixels - the feed offset angles are not
    % mutiples of 0.02 deg so it won't happen for any but the center feed.
    
    % as agreed with Michael
    m.lx=73.5; m.hx=91.5;
    m.ly=-51.5; m.hy=-42.5;
    m.pixsize=0.02;
    % square pixels half way up
    m.xdos=(m.hx-m.lx)*cos(-47*pi/180);
    m.ydos=m.hy-m.ly;
  
  case 'quad0506'    
    % to accomodate both years 
    m.lx=73.5; m.hx=91.5;
    m.ly=-54; m.hy=-42.5;
    m.pixsize=0.02;
    % square pixels half way up
    m.xdos=(m.hx-m.lx)*cos(-47*pi/180);
    m.ydos=m.hy-m.ly;
  
    % using a map area which accomodates 05 and 06 for 06 data alone
    % doesn't work well - m.xdos is used to set the nominal "plate
    % scale" in the 2D fft later - as such it should be the best
    % possible approx to the truth - the way I have it above with -47
    % is very bad as it is at edge of 06 map area leading to
    % significant distortion of ell scale in final spectra. Would be
    % better if it was the mean([m.ly,m.hy])=48.25 and I should have doen
    % that but better still would be to use the true mid point of the
    % actual 06 area which is more like -50
  
  case 'quad06'    
    % 06 area
    m.lx=73.5; m.hx=91.5;
    m.ly=-54; m.hy=-46;
    m.pixsize=0.02;
    % square pixels half way up
    m.xdos=(m.hx-m.lx)*cos(mean([m.ly,m.hy])*pi/180);
    m.ydos=m.hy-m.ly;
    
  case 'bicep'
    m.pixsize=0.25;
    switch proj
      case 'radec'
	m.lx=-55; m.hx=55;
	m.ly=-70; m.hy=-45;
	% square pixels half way up
	m.xdos=(m.hx-m.lx)*cos(mean([m.ly,m.hy])*pi/180);
	m.ydos=m.hy-m.ly;
      otherwise
	m.racen=0;
	m.deccen=-57;
	m.lx=-32; m.hx=32;
	m.ly=-17.5; m.hy=17.5;
	m.xdos=m.hx-m.lx;
	m.ydos=m.hy-m.ly;
    end
 
 case 'bicepext'
   % extend the traditional grid to accomodate BICEP3 2016 obs
   % while maintaining the el at which the pixels are square
   % - idea is to generate a superset of the traditional 'bicep'
   % defn
   % because 2016 B3 obs are shifted down by 2 deg versus
   % traditional el step range the up/down growth is assymetric -
   % we extend up by 3 deg and down by 7
   m.pixsize=0.25;
    switch proj
      case 'radec'
        %   % bixepext defn (approximate)
        %   lx1 = -60; hx1 = 60;
        %   xdos1 = (hx1-lx1)*cos(mean([-70,-45])*pi/180);
        %   sx1 = (hx1-lx1) / round(xdos1/m.pixsize);
        %   % bicep defn
        %   lx2 = -55; hx2 = 55;
        %   xdos2 = (hx2-lx2)*cos(mean([-70,-45])*pi/180);
        %   sx2 = (hx2-lx2) / round(xdos2/m.pixsize);
        %   % accumulate map size mismatch cause by pixel size diff
        %   dsx = (sx2-sx1) * round(xdos1/m.pixsize)/2;
        %   % eliminate acummulated difference by adjusting bounds.
        %   m.lx=lx1-dsx; m.hx=hx1+dsx;
        % constant retrieved by executing the above code:
        m.lx=-60.12711864406779;
        m.hx=+60.12711864406779;
        m.ly=-73; m.hy=-38;
        % square pixels in the traditional row
        m.xdos=(m.hx-m.lx)*cos(mean([-70,-45])*pi/180);
        m.ydos=m.hy-m.ly;
      otherwise
        error('not defined');
    end   
 
  case 'bicepext1'
   % extend the traditional grid down in el while maintaining the
   % el at which the pixels are square
   m.pixsize=0.25;
    switch proj
      case 'radec'
	m.lx=-55; m.hx=55;
	m.ly=-70; m.hy=-35;
	% square pixels in the traditional row
	m.xdos=(m.hx-m.lx)*cos(mean([-70,-45])*pi/180);
	m.ydos=m.hy-m.ly;
      otherwise
        error('not defined');
    end
    
  case 'bicepext2'
   % extend the traditional grid down in el while maintaining the
   % el at which the pixels are square
   m.pixsize=0.25;
    switch proj
      case 'radec'
	m.lx=-55; m.hx=55;
	m.ly=-70; m.hy=-40;
	% square pixels in the traditional row
	m.xdos=(m.hx-m.lx)*cos(mean([-70,-45])*pi/180);
	m.ydos=m.hy-m.ly;
      otherwise
        error('not defined');
    end

  % this is for simulation of spider (Odea) 'cleanest' patch 
  case 'spider'
    m.pixsize=0.25;
    switch proj
      case 'radec'
	m.lx=5; m.hx=115;
	m.ly=-70; m.hy=-45;
	% square pixels half way up
	m.xdos=(m.hx-m.lx)*cos(mean([m.ly,m.hy])*pi/180);
	m.ydos=m.hy-m.ly;
      otherwise
	m.racen=60;
	m.deccen=-57;
	m.lx=-32; m.hx=32;
	m.ly=-17.5; m.hy=17.5;
	m.xdos=m.hx-m.lx;
	m.ydos=m.hy-m.ly;
    end
    
  % Same as type 'bicep', but with 2x finer pixel size.
  % To be used for fitting beam centers
  % and radio pointing model.
  % Not quite the same as resmult=2
  case 'bicepfine'
    m.pixsize=0.125;
    switch proj
      case 'radec'
        m.lx=-55; m.hx=55;
        m.ly=-70; m.hy=-45;
        % square pixels half way up
        m.xdos=(m.hx-m.lx)*cos(mean([m.ly,m.hy])*pi/180);
        m.ydos=m.hy-m.ly;
      otherwise
        m.racen=0;
        m.deccen=-57;
        m.lx=-32; m.hx=32;
        m.ly=-17.5; m.hy=17.5;
        m.xdos=m.hx-m.lx;
        m.ydos=m.hy-m.ly;
    end

  % Same as type 'bicep', but with coarser pixel size.
  % To be used for matrix operations
  % Now redundant given resmult=0.5
  case 'bicepcoarse'
    m.pixsize=0.5;
    switch proj
      case 'radec'
         m.lx=-55; m.hx=55;
         m.ly=-70; m.hy=-45;
         % square pixels half way up
         m.xdos=(m.hx-m.lx)*cos(mean([m.ly,m.hy])*pi/180);
         m.ydos=m.hy-m.ly;
       otherwise
         m.racen=0;
         m.deccen=-57;
         m.lx=-32; m.hx=32;
         m.ly=-17.5; m.hy=17.5;
         m.xdos=m.hx-m.lx;
         m.ydos=m.hy-m.ly;
    end
         
  % galaxy faint arm
  case 'bicepgalf'
    m.pixsize=0.25;
    switch proj
      case 'radec'
	m.lx=80; m.hx=190;
	m.ly=-70; m.hy=-40;
	% square pixels half way up
	m.xdos=(m.hx-m.lx)*cos(mean([m.ly,m.hy])*pi/180);
	m.ydos=m.hy-m.ly;
      otherwise
	m.racen=9;
	m.deccen=-55;
	m.ly=-17.5; m.hy=17.5;
	m.xdos=m.hx-m.lx;
	m.ydos=m.hy-m.ly;
    end
 
  % galaxy bright arm
  case 'bicepgalb'
    m.pixsize=0.25;
    switch proj
      case 'radec'
	m.lx=180; m.hx=291;
	m.ly=-70; m.hy=-40;
	% square pixels half way up
	m.xdos=(m.hx-m.lx)*cos(mean([m.ly,m.hy])*pi/180);
	m.ydos=m.hy-m.ly;
      otherwise
	m.racen=15.712;
	m.deccen=-55;
	m.ly=-17.5; m.hy=17.5;
	m.xdos=m.hx-m.lx;
	m.ydos=m.hy-m.ly;
    end
      
  case 'bicepboomlmc'
    m.pixsize=0.25;
    switch proj
      case 'radec'
	m.lx=12.5; m.hx=122.5;
	m.ly=-75; m.hy=-40;
	% square pixels half way up
	m.xdos=(m.hx-m.lx)*cos(mean([m.ly,m.hy])*pi/180);
	m.ydos=m.hy-m.ly;
      otherwise
	m.racen=67.5;
	m.deccen=-55;
	m.ly=-17.5; m.hy=17.5;
	m.xdos=m.hx-m.lx;
	m.ydos=m.hy-m.ly;
    end  
    
  case 'bicepgalgap'
    m.pixsize=0.25;
    switch proj
      case 'radec'
	m.lx=140; m.hx=250;
	m.ly=-75; m.hy=-45;
	% square pixels half way up
	m.xdos=(m.hx-m.lx)*cos(mean([m.ly,m.hy])*pi/180);
	m.ydos=m.hy-m.ly;
      otherwise
	m.racen=200;
	m.deccen=-65;
	m.ly=-17.5; m.hy=17.5;
	m.xdos=m.hx-m.lx;
	m.ydos=m.hy-m.ly;
    end  
    
  case 'bicepfullgal'
    m.pixsize=0.25;
    switch proj
      case 'radec'
        m.lx=80; m.hx=290;
	m.ly=-75; m.hy=-40;
	% square pixels half way up
	m.xdos=(m.hx-m.lx)*cos(mean([m.ly,m.hy])*pi/180);
	m.ydos=m.hy-m.ly;
      otherwise
	m.racen=150;
	m.deccen=-65;
	m.ly=-17.5; m.hy=17.5;
	m.xdos=m.hx-m.lx;
	m.ydos=m.hy-m.ly;
    end  
    
  case 'biceplastgap'
    m.pixsize=0.25;
    switch proj
      case 'radec'
        m.lx=245; m.hx=355;
        m.ly=-70; m.hy=-44;
        % square pixels half way up
        m.xdos=(m.hx-m.lx)*cos(mean([m.ly,m.hy])*pi/180);
        m.ydos=m.hy-m.ly;
      otherwise
        m.racen=300;
        m.deccen=-55;
        m.ly=-17.5; m.hy=17.5;
        m.xdos=m.hx-m.lx;
        m.ydos=m.hy-m.ly;
      end  
       
  case 'bicepfull'
    m.pixsize=0.25;
    switch proj
      case 'radec'
	m.lx=0.01; m.hx=360;
	m.ly=-80; m.hy=-40;
	% square pixels half way up
	m.xdos=(m.hx-m.lx)*cos(mean([m.ly,m.hy])*pi/180);
	m.ydos=m.hy-m.ly;
      otherwise
	m.racen=150;
	m.deccen=-65;
	m.ly=-17.5; m.hy=17.5;
	m.xdos=m.hx-m.lx;
	m.ydos=m.hy-m.ly;
    end  
         
  %########################
  case 'spud'
    m.lx=-50; m.hx=50;
    m.ly=-70; m.hy=-45;
    m.pixsize=0.1;
    % square pixels half way up
    m.xdos=(m.hx-m.lx)*cos(mean([m.ly,m.hy])*pi/180);
    m.ydos=m.hy-m.ly;
          
  case 'galsurv'
    m.lx=230; m.hx=265;
    m.ly=-56; m.hy=-38;
    m.pixsize=0.02;
    % square pixels half way up
    m.xdos=(m.hx-m.lx)*cos(mean([m.ly,m.hy])*pi/180);
    m.ydos=m.hy-m.ly;
        
  case {'cena' 'carneb' 'galcen' 'g260.4' 'rcw38' 'mat6a' 'qso' 'j0538-440' 'g326.3' 'g328.4'}
    % this is for leakage simulations - they use the whole map
    [ra,dec]=src_coords(type);
    cosmd=cos(dec*pi/180);
    halfside=2.25; % degrees per half side of box
    m.lx=ra-halfside/cosmd; m.hx=ra+halfside/cosmd;
    m.ly=dec-halfside; m.hy=dec+halfside;
    m.pixsize=0.02;
    % square pixels half way up
    m.xdos=(m.hx-m.lx)*cosmd;
    m.ydos=m.hy-m.ly;
   
  case 'test'
    % this is for testing leakage simulations
    type='rcw38';
    [ra,dec]=src_coords(type);
    cosmd=cos(dec*pi/180);
    halfside=2.25; % degrees per half side of box
    m.lx=ra-halfside/cosmd; m.hx=ra+halfside/cosmd;
    m.ly=dec-halfside; m.hy=dec+halfside;
    m.pixsize=0.02;
    % square pixels half way up
    m.xdos=(m.hx-m.lx)*cosmd;
    m.ydos=m.hy-m.ly;
   
  case 'az'
    m.lx=-180; m.hx=180;
    m.ly=-70; m.hy=-45;
    m.pixsize=0.25;
    % square pixels half way up
    m.xdos=(m.hx-m.lx)*cos(mean([m.ly,m.hy])*pi/180);
    m.ydos=m.hy-m.ly;
    
  case 'azel'
    m.lx=-180; m.hx=180;
    m.ly=45; m.hy=70;
    m.pixsize=0.25;
    % square pixels half way up
    m.xdos=(m.hx-m.lx)*cos(mean([m.ly,m.hy])*pi/180);
    m.ydos=m.hy-m.ly;
    
  % test case for using galaxy maps for 2015 keck beamcenters
  case 'bicepgalbfine'
    m.pixsize=0.125;
    switch proj
      case 'radec'
        m.lx=180; m.hx=291;
        m.ly=-70; m.hy=-40;
        % square pixels half way up
        m.xdos=(m.hx-m.lx)*cos(mean([m.ly,m.hy])*pi/180);
        m.ydos=m.hy-m.ly;
      otherwise
        m.racen=15.712;
        m.deccen=-55;
        m.ly=-17.5; m.hy=17.5;
        m.xdos=m.hx-m.lx;
        m.ydos=m.hy-m.ly;
    end

  case 'fttest_even'
    % very low res for testing ft algorithms etc
    m.lx=-10; m.hx=10;
    m.ly=-10; m.hy=10;
    m.xdos=m.hx-m.lx;
    m.ydos=m.hy-m.ly;
    m.pixsize=m.xdos/4;
    
  case 'fttest_odd'
    % very low res for testing ft algorithms etc
    m.lx=-10; m.hx=10;
    m.ly=-10; m.hy=10;
    m.xdos=m.hx-m.lx;
    m.ydos=m.hy-m.ly;
    m.pixsize=m.xdos/5;
    
  otherwise
    error('unknown map defn type');
   
end

m.nx=round(m.xdos/m.pixsize);
m.ny=round(m.ydos/m.pixsize);

% do it this way to ensure pixel supersets for odd values
m.nx=m.nx*abs(resmult);
m.ny=m.ny*abs(resmult);
m.pixsize=m.pixsize/abs(resmult);

% calc pixel spacing
sx=(m.hx-m.lx)/m.nx;
sy=(m.hy-m.ly)/m.ny;

% if resmult negative use this as indication to provide a grid
% which is shifted so that FT is a superset - see bkspt posting
% http://bicep.rc.fas.harvard.edu/bkspt/analysis/20160713_comb_bk_and_spt
if(resmult<0)
  mf=abs(resmult)/2-0.5
  m.lx=m.lx+mf*sx; m.hx=m.hx+mf*sx;
  m.ly=m.ly+mf*sy; m.hy=m.hy+mf*sy;  
end

% calc the pixel centers in same way as hfill2 will later
m.x_tic=m.lx+sx/2:sx:m.hx-sx/2;
m.y_tic=m.ly+sy/2:sy:m.hy-sy/2;
return
