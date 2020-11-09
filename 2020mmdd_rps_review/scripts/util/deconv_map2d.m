function [map,ad]=deconv_map2d(m,map,apsu2d,apsf2d,apsn2d,lpass,ap,deap,inv)
% [map,ad]=deconv_map2d(m,map,apsu2d,apsf2d,apsn2d,lpass,ap,deap,inv)
%
% Apodize, take map to Fourier space and, apply 2d supfac correction
% and 2d Wiener filter. Then go back to map space and deapodize.
%
% m = map defn
%
% map = mapfile name or preloaded data
%
% apsu2d/apsf2d/apsn2d = mean 2d F-modes from reduc_makeaps2d for
%   unfiltered signal, filtered signal and noise respectively
%
% lpass = highest ell modes which will be considered. default to
%   radius of F-plane.
% 
% ap = apply apodization before FT. default true. need to do
%   it unless maps are already apodized
%
% deap = remove apodization (and reimpose zero coverage regions)
%   after reverse FT. default false. May not wish to do this when
%   generating lensing template.
%
% inv = reverse filtering to undo action of a previous call.
%   default false. Filter is remembered in the map structure
%   so apsu2d,apsf2d,apsn2d,lpass args have no effect in this case
%
% Information about the applied filtering is passed back as
% additional fields of the map structure
%
% e.g.
% load maps/1459/real_aabd_filtp3_weight3_gs_dp1102_jack0.mat
% map=make_map(ac,m);
% plot_map(m,map(2).Q); caxis([-4,4]);
% apsu2d='aps2d/1459/0xx5_aabd_unfilt_nobeam_jack0';
% apsf2d='aps2d/1459/0xx5_aabd_filtp3_weight3_gs_dp1100_jack0';
% apsn2d='aps2d/1459/0xx6_aabd_filtp3_weight3_gs_dp1100_jack0';
% [mapd,ad]=deconv_map2d(m,map,apsu2d,apsf2d,apsn2d,700,1,0);
% plot_map(m,mapd(2).Q); caxis([-4,4]);
% mapdd=deconv_map2d(m,mapd,apsu2d,apsf2d,apsn2d,700,0,1,1);
% plot_map(m,mapdd(2).Q); caxis([-4,4]);

if(~exist('ap','var')) || isempty(ap)
  ap=true;
end  
if(~exist('deap','var')) || isempty(deap)
  deap=false;
end
if(~exist('inv','var')) || isempty(inv)
  inv=false;
end

% if map filename provided load it
if(ischar(map))
  load map
  % if it contains ac make the map
  if(exist(ac,'var'))
    map=make_map(ac,m);
  end
end

% if filenames load
if(ischar(apsu2d))
  load(apsu2d);
  U=aps2d;
end
if(ischar(apsf2d))
  load(apsf2d);
  F=aps2d;
end
if(ischar(apsn2d))
  load(apsn2d);
  N=aps2d;
end

% if ap masks are not already available make std ones using var
% fields of map
if(~isfield(map,'Pw'))
  map=add_masks(m,map);
end

% calculate the ad info
ad=calc_ad2([m.xdos,m.ydos],[m.nx,m.ny]);
ad.l=ad.u_r*2*pi;

% assign default value to lpass if needed (see below)
if(~exist('lpass','var')) || isempty(lpass)
  lpass=abs(ad.u_val{1}(1)*2*pi);
end

% set the peak of the weight maps to unity so the normalization in
% the central region does not change much - this is how it is done
% in reduc_plotcomap_pager - this may need further thought
for i=1:numel(map)
  map(i).Tw=map(i).Tw./max(map(i).Tw(:));
  map(i).Pw=map(i).Pw./max(map(i).Pw(:));
end

% get FT's
clear m
for i=1:numel(map)
  if(ap)  
    % get map Fourier modes (the 'false' arg means don't normalize the
    % modes to correct for masking)
    m(i)=calc_map_fts(map(i),ad,map(i).Pw,[2],'normal',false);
    mt(i)=calc_map_fts(map(i),ad,map(i).Tw,[1],'normal',false);
  else
    % do it without applying apodization
    m(i)=calc_map_fts(map(i),ad,1,[2],'normal',false);
    mt(i)=calc_map_fts(map(i),ad,1,[1],'normal',false);
  end
end
for i=1:numel(map)
  m(i).T=mt(i).T;
end

if(~inv)
  % The input map used to make the unfiltered sims needs to
  % have had an ell cutoff at the pixel scale to prevent aliasing
  % - hence the modes "in the corners" of the F-plane can have non-zero
  % values only due to mixing due to ap mask. Therefore we force
  % them to zero
  for i=1:numel(map)
    F(i).T(ad.l>lpass)=0;
    F(i).E(ad.l>lpass)=0;
  end
  
  % calc the 2d supfacs
  for i=1:numel(map)
    map(i).RT=F(i).T./U(i).T;
    map(i).RE=F(i).E./U(i).E;
  end
  
  % where less than 1% in amplitude survives filtering
  % assume the mode is not recoverable
  for i=1:numel(map)
    map(i).RT(map(i).RT<1e-4)=0;
    map(i).RE(map(i).RE<1e-4)=0;
  end
  
  % correct the 2d noise
  for i=1:numel(map)
    % divide straight out of 2d noise aps
    Np(i).T=N(i).T./map(i).RT;
    Np(i).Q=N(i).Q./map(i).RE;
    Np(i).U=N(i).U./map(i).RE;
    Np(i).E=N(i).E./map(i).RE;
  end
  
  % make the Weiner filter
  for i=1:numel(map) 
    map(i).WT=U(i).T./(U(i).T+Np(i).T);
    map(i).WE=U(i).E./(U(i).E+Np(i).E);
  end
  
  % make the overall filter
  for i=1:numel(map)
    map(i).OT=map(i).WT./sqrt(map(i).RT);
    map(i).OE=map(i).WE./sqrt(map(i).RE);
  end
  
  % where overall filter is NaN make it zero
  % (1% amplitude threshold above leads to 0/0=NaN)
  for i=1:numel(map)
    map(i).OT(isnan(map(i).OT))=0;
    map(i).OE(isnan(map(i).OE))=0;
  end
  
else
  % invert the overall filter previously calculated
  for i=1:numel(map)
    map(i).OT=1./map(i).OT;
    map(i).OE=1./map(i).OE;
  end
  
  % remove infinite values resulting from 1/0
  for i=1:numel(map)
    map(i).OT(isinf(map(i).OT))=0;
    map(i).OE(isinf(map(i).OE))=0;
  end
end

% apply the overall filter
for i=1:numel(map)
  mpp(i).T=m(i).T.*map(i).OT;
  mpp(i).Q=m(i).Q.*map(i).OE;
  mpp(i).U=m(i).U.*map(i).OE;
  mpp(i).E=m(i).E.*map(i).OE;
end

% go back to image plane
for i=1:numel(map)
  map(i).T=f2i(ad,mpp(i).T);
  map(i).Q=f2i(ad,mpp(i).Q);
  map(i).U=f2i(ad,mpp(i).U);
end

% if requested deapodize
if(deap)
  map(i).T=map(i).T./map(i).Tw;
  map(i).Q=map(i).Q./map(i).Pw;
  map(i).U=map(i).U./map(i).Pw;
  
  % and reimpose same coverage as initial
  map(i).T(isnan(map(i).Tvar))=NaN;
  map(i).Q(isnan(map(i).Qvar))=NaN;
  map(i).U(isnan(map(i).Uvar))=NaN;
end

return
