function [mp,hmapcel,hmapgal]=plot_pntsrcmask(src,type,plottype)
% visualize a pointsource mask as a true/false ra/dec map
%
%  [mask,hmap,hmap_gal]=disp_pntsrcmask(src,type,plottype)
%
%   src is same as in reduc_pntsrcmask and get_mask
%
%   type is maptype input to get_map_defn, e.g. 'bicep', 'bicepgalp',
%   etc. 
%
%   if no output arguments are present, produces a plot with
%   plottype=1,2,3 for ra/dec, celestial healpix, galactic healpix 
%

if(~exist('plottype','var'))
  plottype=1;
end

m=get_map_defn(type,0,'radec');
% make a higher res lmap for interpolating onto a healpix grid
if(plottype==2|plottype==3)
  m.pixsize=0.1;
  m.y_tic=m.ly:m.pixsize:m.hy;
  m.x_tic=m.lx:(m.pixsize/cosd(mean([m.ly,m.hy]))):m.hx;
  m.nx=length(m.x_tic);
  m.ny=length(m.y_tic);
end

[x,y]=meshgrid(m.x_tic,m.y_tic);

% check which map coordinates are masked in the same way as
% reduc_pntsrcmask
mp=get_mask(x,y,src);
hmapcel=longlat_to_healpix(mp,x,y,256);

[l,b]=euler(x,y,1,0);
hmapgal=longlat_to_healpix(mp,l,b,256);

i=1;j=2;k=1;
if(nargout==0)
  % set up FDS healpix map
  hmapfds=read_fits_map('~/bicep_analysis/skymaps/SFD_i100_healpix_256_fwhm55.2.fits');
  lmapfds=healpix_to_longlat(hmapfds,m.x_tic,m.y_tic,1);
  switch plottype
   case 1
    subplot(2,2,1)
    imagesc(m.x_tic,m.y_tic,log10(lmapfds));colorbar;
    title('FDS map log10(T)');c=caxis;
    subplot(2,2,2)
    im=imagesc(m.x_tic,m.y_tic,log10(lmapfds).*abs(mp-1));colorbar;
    title('FDS masked map log10(T)');caxis(c);
    %set(im,'AlphaData',abs(mp*.5-1));
    i=3;j=2;k=2;    
    
    dir='/home/bicep0/dbarkats/bicep_analysis/maps/save/';

    switch type
     case 'bicepgalp'
      load([dir,'fullseason_galp_filtp3_weight2_jack0.mat'],'map');
     case 'bicepgalgap'
      load([dir,'fullseason_galgap_filtp3_weight2_jack0.mat'],'map');
     case 'bicepcmb3'
      load([dir,'fullseason_cmb3_filtp3_weight2_jack0.mat'],'map');
     case 'bicep'
      load([dir,'fullseason_filtp3_weight0_jack0.mat'],'map');
     otherwise
      map(1).T=ones(size(mp));
      map(2).T=map(1).T;
    end
    subplot(j,k,i)
    imagesc(m.x_tic,m.y_tic,map(2).T);caxis([-4,4]*1e-6);colorbar;
    xlabel('RA');ylabel('DEC');title('BICEP map (if available)');
    plotmask=ones(size(mp)); 
    %plotmask(mp==1)=NaN;
    plotmask=abs(mp-1);
    subplot(j,k,i+1);
    imagesc(m.x_tic,m.y_tic,map(2).T.*plotmask);caxis([-4,4]*1e-6);
    colorbar;xlabel('RA');ylabel('DEC');title('BICEP masked map');
    
   case 2
    % convert fds map to celestial
    [theta,phi]=pix2ang_ring(hmapfds.nside,hmapfds.pixel);
    ra=phi*180/pi; dec=-theta*180/pi+90;
    [l,b]=euler(ra,dec,1,0);
    phi=l*pi/180; theta=(-b+90)*pi/180;
    pix=ang2pix_ring(hmapfds.nside,theta,phi);
    temp=hmapfds.map;
    hmapfds.map(:)=NaN;
    [pix,ind]=unique(pix);
    hmapfds.map(pix+1)=temp(ind);
    while length(find(isnan(hmapfds.map)==1))>1
      ind=find(isnan(hmapfds.map)==1);
      ind(ind==1)=2;
      hmapfds.map(ind)=hmapfds.map(ind-1);
    end

    hmapfds.map=log10(hmapfds.map);
    figure(1)
    mollview(hmapfds);c=caxis;title('celestial log(T)');colorbar
    hmap=hmapcel;
    hmap.map=abs(hmap.map-1);    
    hmap.map(hmap.map==0)=NaN;
    hmap.map=hmap.map.*hmapfds.map;
    figure(2)
    mollview(hmap);title('celestial log(T)');caxis(c);colorbar

   case 3
    hmapfds.map=log10(hmapfds.map);
    figure(1)
    mollview(hmapfds);c=caxis;title('galactic log(T)');colorbar;
    hmap=hmapgal;
    hmap.map=abs(hmap.map-1);    
    hmap.map(hmap.map==0)=NaN;
    hmap.map=hmap.map.*hmapfds.map;
    figure(2)
    mollview(hmap);title('galactic log(T)');caxis(c);colorbar;
  end
end

return


function hmap=longlat_to_healpix(lmap,x,y,nside)

if(~exist('nside','var'))
  nside=128;
end
  
hmap.nside=nside;
hmap.ordering='RING';
npix=12*nside^2;
hmap.pixel=cvec(0:(npix-1));
hmap.map=ones(npix,1)*NaN;

phi=x*pi/180;
theta=(-y+90)*pi/180;
pix=ang2pix_ring(nside,theta,phi);
hmap.map(pix+1)=cvec(lmap);

return

  
