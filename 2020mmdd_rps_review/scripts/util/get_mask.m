function mask=get_mask(x,y,src)
% given ra and dec pairs, returns a logical array of same size where true
% denotes data to be excluded
%    mask=getmask(ra,dec,src)
%
%    src.name = 'galstrips' 
%                  masks data that lies within 0.5 deg in dec
%                  and within src.r(i) degrees in the az direction of each
%                  coordinate specified by src.ra(i) and src.dec(i)
%               'galplane' 
%                  masks points with galactic |b| < src.r
%               'fds'
%                  masks points that in the FDS map have T/max(T)>1/src.r
%               'pntsrc'
%                  other names default to pointsource masking within src.r(i)
%                  degrees of each coordinate specfied by src.ra and src.dec
%    src.r = see above
%

mp=false(size(x));

switch src.name
 case 'galstrips'
  for k=1:length(src.ra)
    % huge number of src for galsurv because we effectively have one
    % per dec pixel - a huge number given the dec range.
    % Presumably this makes it too slow to use spaceangle so Tom did
    % this instead...
    
    % select portion of x and y in reasonable proximity to this
    % slice of the galactic plane
    % xi are the indices into x and y for which the scan track
    % a) lies within 0.5 deg in dec of the feed dec 
    % b) lies within src.r of the galactic plane in the az direction
    xi=find(abs(src.ra(k)-x)<src.r(k)&abs(src.dec(k)-y)<0.5);
    
    mp(xi)=true;
  end
  
  
 case 'galplane'
  % make grid of l/b and interpolate to timestream
  [l,b]=euler(x,y,1,0);
  mp(abs(b)<src.r)=true;

  
 case 'fds'
  % read in FDS map and mask out regions with T > max(T)/r
  dy=0.1;
  decg=(min(y(:))-1):dy:(max(y(:))+1);
  dx=dy/abs(sind(mean([decg(1),decg(end)])))
  rag=(min(x(:))-1):dx:(max(x(:)+1));
  
  hmap=read_fits_map('~/bicep_analysis/skymaps/SFD_i100_healpix_256_fwhm55.2.fits');
  
  lmap=log10(healpix_to_longlat(hmap,rag,decg,1));
  msk=false(size(lmap));
  msk(lmap>max(lmap(:))*(src.r/100))=true;
  mp=interp2(rag,decg,msk,x,y,'nearest');
  
 case 'pntsrc'
  for k=1:length(src.ra)
    % find distance on sky from feed to source
    s=spaceangle(x,y,src.ra(k),src.dec(k),'deg');
    % where angle greater than threshold and current deck angle set the
          % mask to true for current channel - true means exclude data
    mp(s<src.r(k))=true;
  end

end

mask=mp;

return