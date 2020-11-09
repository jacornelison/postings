function [map,ad]=deconv_map(m,map,beam,lpass,apss,apsn,ap,deap,inv,wienereb)
% [map,ad]=deconv_map(m,map,beam,lpass,apss,apsn,ap,deap,inv)
%
% Apodize, take map to Fourier space and, as requested, apply any
% or all of beam deconv, Wiener filter and low pass filter. Then go
% back to map space and, if requested, deapodize.
%
% m = map defn
%
% map = mapfile name or preloaded data
%
% beam = beamfile name or preloaded data, set to 'none' for no
% correction. elements<1e-2 are set to zero - beam correction >100x
% doesn't make any practical sense.
%
% lpass = 2 element vector specifying lower and upper limits of
% cos^2 rolloff low pass filter in units of ell. Empty means
% default value. Default value is upper edge of roll at the lowest
% of 650, top ell bin of the aps provided for Wiener filter, or the
% point at which the beam hits zero (actually 1e-2 given note above).
% For B2/Keck 650 is the highest value that would ever make sense
% to allow through given the 3Hz low pass timestream filter (see
% Chris post 20150331_delensing)
%
% apss,apsn = apsset simset names or preloaded data to perform
% Wiener filter
%
% ap = apodize map prior to FFTing? Default true, but might want false if input maps is
%   already apodized, as would be the case for a second call to this function to re-apply beam smoothing after
%   delensing the output of a first call to this function
%
% deap = undo apodization after return to image plane (this is
% probably not wanted for practical purposes)
%
% inv = apply inverse of filter (i.e. reconvolve instead of deconvolve) default = false
%
% Full information about the applied filtering is passed back as
% additional fields of the ad structure
%
% e.g.
% load maps/1459/real_aabd_filtp3_weight3_gs_dp1102_jack0.mat
% map=make_map(ac,m);
% load aps/1459x1615/xxx6_aabd_filtp3_weight3_gs_dp1100_jack0_xxx6_aabd_filtp3_weight3_gs_dp1100_jack0_pureB_matrix_cm.mat
% n=aps;
% load aps/1459x1615/xxx5_aabd_filtp3_weight3_gs_dp1100_jack0_xxx5_aabd_filtp3_weight3_gs_dp1100_jack0_pureB_matrix_cm.mat
% s=aps;
% mapd=deconv_map(m,map(2),[],[],s(2),n(2))
% plot_map(m,mapd.Q); caxis([-4,4])

% fill in args which have defaults
if(~exist('beam','var')) || isempty(beam)
  % Standard B2 beam
  beam='aux_data/beams/beamfile_20130222_sum.fits';
end
if(~exist('lpass','var'))
  lpass=[];
end
if(~exist('apss','var'))
  apss=[];
end
if(~exist('apsn','var'))
  apsn=[];
end
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
  load(map)
end

% calculate the 2d ell values
ad=calc_ad2([m.xdos,m.ydos],[m.nx,m.ny]);
ad.l=ad.u_r*2*pi;

% initial filters are unity - we will multiply in some or all of the
% three filters
ad.tf.Tf2d=ones(size(ad.l));
ad.tf.Pf2d=ones(size(ad.l));

% if requested unfold beam
if(~strcmp(beam,'none'))
  
  % get beam
  if(ischar(beam))
    beam=fitsread(beam,'bintable');
    ad.beam.f=beam{1}; ad.beam.l=0:(length(ad.beam.f)-1);
  else
    ad.beam.l=0:1500;
    [l,bl]=get_bl(beam,ad.beam.l);
    ad.beam.f=bl;
  end
  
  % don't allow the beam to go below 1e-2 - beam correct of greater
  % than factor 100 doesn't make much sense
  ad.beam.f(ad.beam.f<1e-2)=0;

  % interp to 2d
  ad.beam.f2d=interp1(ad.beam.l,ad.beam.f,ad.l);

  % mult into total filter
  ad.tf.Tf2d=ad.tf.Tf2d.*(1./ad.beam.f2d);
  ad.tf.Pf2d=ad.tf.Pf2d.*(1./ad.beam.f2d);  
end

% if requested apply Wiener filter - a filter which varies
% between zero and one depending on the s/n of the Fourier modes -
% note this means that the output map will be a biased
% representation of the true sky
if(~isempty(apss))
  
  % if passed filenames load
  if(ischar(apss))
    load apss
    apss=aps;
  end
  if(ischar(apsn))
    load apsn
    apsn=aps;
  end
  
  % get the mean of signal and noise sims for T and E - we use E to
  % construct the pol filter
  ad.wf.l=apss.l;
  ad.wf.Ts=mean(apss.Cs_l(:,1,:),3);
  ad.wf.Tn=mean(apsn.Cs_l(:,1,:),3);
  ad.wf.Ps=mean(apss.Cs_l(:,3,:),3);
  ad.wf.Pn=mean(apsn.Cs_l(:,3,:),3);
  
  % construct Wiener filters - note that since we are working with
  % power spectra there is no need to square
  ad.wf.Tf=ad.wf.Ts./(ad.wf.Ts+ad.wf.Tn);
  ad.wf.Pf=ad.wf.Ps./(ad.wf.Ps+ad.wf.Pn);
  
  % expand Wiener filters to 2d
  ad.wf.Tf2d=interp1(ad.wf.l,ad.wf.Tf,ad.l,[],0);
  ad.wf.Pf2d=interp1(ad.wf.l,ad.wf.Pf,ad.l,[],0);
  
  % mult into total filter
  ad.tf.Tf2d=ad.tf.Tf2d.*ad.wf.Tf2d;
  ad.tf.Pf2d=ad.tf.Pf2d.*ad.wf.Pf2d;
end

if inv
  % Invert Wiener and beam filters before low pass 
  ad.tf.Tf2d=1./ad.tf.Tf2d;
  ad.tf.Pf2d=1./ad.tf.Pf2d;
end

% roll off smoothly at high ell
if(~strcmp(lpass,'none'))

  if(isempty(lpass))
    % default to lower of beam=0 or wiener last bin locations
    llb=Inf; llw=Inf;
    if(isfield(ad,'beam'))
      llb=min(find(ad.beam.f==0));
    end
    if(isfield(ad,'wf'))
      llw=floor(apss.l(end));
    end
    ul=min([650,llb,llw]);
    lpass(1)=ul-50;
    lpass(2)=ul;
  end
  
  rr=diff(lpass);
  % tukeywin is cos^2 roll off
  h=tukeywin(2*rr,1); h=h(rr:end);
  ad.lpass.f=[ones(lpass(1),1);h];
  ad.lpass.l=0:(length(ad.lpass.f)-1);
  ad.lpass.f2d=interp1(ad.lpass.l,ad.lpass.f,ad.l,[],0);

  % mult into total filter
  ad.tf.Tf2d=ad.tf.Tf2d.*ad.lpass.f2d;
  ad.tf.Pf2d=ad.tf.Pf2d.*ad.lpass.f2d;
end

% if ap masks are not already available make std ones using var
% fields of map
if(~isfield(map,'Pw'))
  map=add_masks(m,map);
end

% set the peak of the weight maps to unity so the normalization in
% the central region does not change much - this is how it is done
% in reduc_plotcomap_pager - this may need further thought
if ap
  map.Tw=map.Tw./max(map.Tw(:));
  map.Pw=map.Pw./max(map.Pw(:));
  Tw=map.Tw; Pw=map.Pw;
else
  Tw=1; Pw=1;
end

% get map Fourier modes (the 'false' arg means don't normalize the
% modes to correct for masking) 
ftt=calc_map_fts(map,ad,Tw,[1],'normal',false);
ftp=calc_map_fts(map,ad,Pw,[2],'normal',false);

% set non finite total filter values to zero
ad.tf.Tf2d(~isfinite(ad.tf.Tf2d))=0;
ad.tf.Pf2d(~isfinite(ad.tf.Pf2d))=0;

% apply total filter
ftt.T=ftt.T.*ad.tf.Tf2d;
if wienereb
  disp('wienereb on')
  [Qf, Uf] = apply_wiener_eb(ftp.Q,ftp.U, ad, ad.tf, 1);
  ftp.Q = Qf;
  ftp.U = Uf;
else
  ftp.Q=ftp.Q.*ad.tf.Pf2d;
  ftp.U=ftp.U.*ad.tf.Pf2d;
end



% keep coverage masks
tm=isnan(map.T);
qm=isnan(map.Q);
um=isnan(map.U);

% go back to image plane
map.T=f2i(ad,ftt.T);
map.Q=f2i(ad,ftp.Q);
map.U=f2i(ad,ftp.U);

% deapodize if requested
if(deap)
  map.T=map.T./map.Tw;
  map.Q=map.Q./map.Pw;
  map.U=map.U./map.Pw;
end

% reimpose same coverage as initial
map.T(tm)=NaN;
map.Q(qm)=NaN;
map.U(um)=NaN;

% extract the 1d total filter for plotting purposes
i=ad.N_pix(1)/2+1; j=ad.N_pix(1); k=ad.N_pix(2)/2+1;
ad.tf.l=ad.l(k,i:j);
ad.tf.Tf=ad.tf.Tf2d(k,i:j);
ad.tf.Pf=ad.tf.Pf2d(k,i:j);

return
