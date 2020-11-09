function ft=calc_map_fts(map,ad,w,t,pure_b,appsf)
% ft=calc_map_fts(map,ad,w,t,pure_b,appsf)
%
% Calculate the fourier components of T/Q/U maps using given weight
% mask and xform to E/B
%
% w = weight mask
% t = which ft's to calc, 1=T, 2=Q/U/E/B, 3=sum/diff, 4=kappa
% pure_b = normal or Kendrick estimator
% appsf = apply or not scale factor to normalize ft

if(~exist('pure_b','var'))
  pure_b='normal';
end
if(~exist('appsf','var'))
  appsf=true;
end

% force the mask to same size as map
% (this must be simply a resolution change)
if(~all(size(ad.t_r)==size(w)))
  disp(sprintf('Changing resolution of mask from %dx%d to %dx%d',size(w),size(ad.t_r)));
  w=interpft2(w,size(ad.t_r));
end

% calculate fourier grids
[u,v]=meshgrid(ad.u_val{1},ad.u_val{2});

if(any(t==1))
  ft.T=getft(ad,map.T,w,appsf);
end

if(any(t==2))  
  % always compute Q/U fts from base Q/U maps
  ft.Q=getft(ad,map.Q,w,appsf); ft.U=getft(ad,map.U,w,appsf);
  
  switch pure_b
   case 'normal'
    % compute E&B by rotating base Q/U fts
    [ft.E,ft.B]=qu2eb(u,v,ft.Q,ft.U,[],[],'iau');
    
   case 'kendrick'
    % compute E&B from base Q/U maps directly
    % (These E modes are identical to the 'normal' ones)
    [ft.E,ft.B]=qu2eb_pure(ad,map.Q,map.U,w,'iau',appsf); 
  end %switch
  
  % if matrix purified Q/U maps exist compute their F-modes and
  % replace B modes with these
  if(isfield(map,'QprojB'))
    if(~isempty(map.QprojB))
      ft.QprojB=getft(ad,map.QprojB,w,appsf); ft.UprojB=getft(ad,map.UprojB,w,appsf);
      [dum,ft.B]=qu2eb(u,v,ft.QprojB,ft.UprojB,[],[],'iau');
    end
  end
  
  % replace E modes with these
  if(isfield(map,'QprojE'))
    if(~isempty(map.QprojE))
      ft.QprojE=getft(ad,map.QprojE,w,appsf); ft.UprojE=getft(ad,map.UprojE,w,appsf);
      [ft.E,dum]=qu2eb(u,v,ft.QprojE,ft.UprojE,[],[],'iau');
    end
  end
  
end %if(any(t==2))

if(any(t==3))
  ft.Tsum=getft(ad,map.Tsum,w,appsf);
  ft.Tdif=getft(ad,map.Tdif,w,appsf);
end

if(any(t==4))
  ft.kappa=getft(ad,map.kappa,w,appsf);
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ft=getft(ad,im,w,appsf)

% de-sparse if not already so
if issparse(im)
  im=full(im);
end
if issparse(w)
  va=full(w);
end

% take product of map and mask
f=im.*w;
f(isnan(f))=0;

% take the ft
ft=i2f(ad,f);

if(appsf)
  % Apply scale factor to undo affect of mask.
  % This forces the application of the mask to leave total power
  % unaffected.
  % This is correct for a Gaussian random field where the map structure
  % is spread equally across the image plane - obviously if the power
  % is localized in a particular area it will not be correct.
  sf=prod(ad.N_pix)/nansum(rvec(w.^2));
  ft=ft*sqrt(sf);
end

if(0)
  figure(3); clf; setwinsize(gcf,1200,300);
  subplot(1,3,1); imagesc(im);
  axis image; axis xy; set(gca,'Xdir','reverse'); colorbar
  title('map')
  
  subplot(1,3,2); imagesc(w);
  axis image; axis xy; set(gca,'Xdir','reverse'); colorbar
  title('mask')
  
  subplot(1,3,3); imagesc(f);
  axis image; axis xy; set(gca,'Xdir','reverse')
  title('map x mask')
  %ca=caxis; caxis(ca/10);
  colorbar
  pause
  %drawnow
end

return
