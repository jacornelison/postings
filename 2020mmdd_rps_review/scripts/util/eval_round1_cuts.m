function [c,cp]=eval_round1_cuts(cp,cut,p,ind)
% [c,cp]=eval_round1_cuts(cp,cut,p,ind)
%
% Evaluate cuts which are applied in reduc_makepairmaps -
% i.e. to remove a subset of the halfscans on
% each channel before they are binned into the pixels.
%
% Using the cut parameters in cp and the values and options in
% structure "cut" generate a series of cut masks - while the cut-params can
% be any size the output cut size must be something that expand_cuts
% knows how to make into nhalf-scan x nchannels.
% Creation of each cut mask is triggered by the existence of a field
% in structure "cut" which specifies a range/threshold etc. Each cut
% can be evaluated using one or more cut parameters according to
% arbitrary rules as coded below.
%
% Once the cuts have been evaluated we expand them to full size and
% then combine to an overall mask. From this additional cut-on-cut
% parameters are generated and cut upon.
%
% p is used to tell which channels are in which rx
%
% Ideally we wouldn't use the g and rg fields of ind here because we
% don't want to decide what they are until reduc_coaddpairmaps, but it
% seems unavoidable

disp('evaluating round1 cuts...');

% copy "utility" parameters
c.nhs=cp.nhs;
c.nch=cp.nch;
c.nrx=cp.nrx;
% nmce added in BICEP3-era, so ensure compatibility with older data products
if isfield(cp,'nmce')
  c.nmce=cp.nmce;
else
  c.nmce=cp.nrx; % True for BICEP2 and Keck
end

% number of nans in each half-scan must be less than threshold
if(isfield(cut,'fb_nancount'))
  c.fb_nancount=cp.fb_nancount<cut.fb_nancount;
end

% field scan std with p0 filtering must be in range
if(isfield(cut,'fb_std_p0'))
  c.fb_std_p0=cp.fb_std_p0>cut.fb_std_p0(1)&...
              cp.fb_std_p0<cut.fb_std_p0(2);
end

% when darks are anomolous cut the whole receiver
if(isfield(cut,'fb_std_p0_darks'))
  for r=unique(p.rx)'
    % find the good dark channels in this rx
    rp=intersect(find((p.rx)==r),ind.gd);
    % find the median parameter value for these channels
    x=nanmedian(cp.fb_std_p0(:,rp),2);
    % cut on parameter out of range
    c.fb_std_p0_darks(:,r+1)=x>=cut.fb_std_p0_darks(1)&...
                             x<cut.fb_std_p0_darks(2);
  end
end		 

% field scan std with p3 filtering must be in range
if(isfield(cut,'fb_std_p3'))
  c.fb_std_p3=cp.fb_std_p3>cut.fb_std_p3(1)&...
              cp.fb_std_p3<cut.fb_std_p3(2);
end

% uncalibrated field scan per-channel std with p3 filtering must be in range
if(isfield(cut,'fb_std_uncal'))
  c.fb_std_uncal=cp.fb_std_uncal>cut.fb_std_uncal(1)&...
                 cp.fb_std_uncal<cut.fb_std_uncal(2);
end

% field scan pair-diff std with p0 filtering must be in range
if(isfield(cut,'fb_std_sd_p0'))
  % this cut applies only to paired channels - true by default
  c.fb_std_sd_p0=true(c.nhs,c.nch);
  % take cut from diff
  val=cp.fb_std_sd_p0(ind.b)>cut.fb_std_sd_p0(1)&...
      cp.fb_std_sd_p0(ind.b)<cut.fb_std_sd_p0(2);
  % and apply to sum and diff
  c.fb_std_sd_p0(ind.a)=val; c.fb_std_sd_p0(ind.b)=val;
end

% field scan pair-diff std with p3 filtering must be in range
if(isfield(cut,'fb_std_sd_p3'))
  % this cut applies only to paired channels - true by default
  c.fb_std_sd_p3=true(c.nhs,c.nch);
  % take cut from diff
  val=cp.fb_std_sd_p3(ind.b)>cut.fb_std_sd_p3(1)&...
      cp.fb_std_sd_p3(ind.b)<cut.fb_std_sd_p3(2);
  % and apply to sum and diff
  c.fb_std_sd_p3(ind.a)=val; c.fb_std_sd_p3(ind.b)=val;
end

% is any pixel in colomn/row flux jumping?
if(isfield(cut,'is_fj_row'))
  c.is_fj_row=cp.is_fj_row<cut.is_fj_row;
end
if(isfield(cut,'is_fj_col'))
  c.is_fj_col=cp.is_fj_col<cut.is_fj_col;
end

% timing error cuts
if(isfield(cut,'syncsampnum_diff1'))
  c.syncsampnum_diff1=abs(cp.syncsampnum_diff1)<=cut.syncsampnum_diff1;
end
if(isfield(cut,'syncsampnum_diff2'))
  c.syncsampnum_diff2=abs(cp.syncsampnum_diff2)<=cut.syncsampnum_diff2;
end

% do the expand here because we need to do a combine below
c=expand_cuts(c,p);


% DO PASSFRAC COL/CHAN CUTS

% calculate overall first round mask
cm=combine_cuts(c);

% find the fraction of nominally good channels passing in each column
seq_col_list=unique(p.sequential_mce_col(isfinite(p.sequential_mce_col)));
for seq_col=seq_col_list'
  % find the mask values for the good light channels in this
  % column
  x=cm(:,intersect(find((p.sequential_mce_col)==seq_col),ind.gl));
  if ~isempty(x)
    % find the fraction passing
    cp.passfrac_col(:,seq_col)=sum(x,2)/size(x,2);
  end
end

% find the fraction of nominally good channels passing in each receiver
for r=unique(p.rx)'
  % find the mask values for the good light channels in this rx
  x=cm(:,intersect(find((p.rx)==r),ind.gl));
  if ~isempty(x)
    % find the fraction passing
    cp.passfrac_chan(:,r+1)=sum(x,2)/size(x,2);
  end
end
  
% apply the column fraction cut
if(isfield(cut,'passfrac_col'))
  c.passfrac_col=true(size(cm));
  seq_col_list=unique(p.sequential_mce_col(isfinite(p.sequential_mce_col)));
  for seq_col=seq_col_list'
    i=find(p.sequential_mce_col==seq_col);
    if ~isempty(i)
      v=cp.passfrac_col(:,seq_col)>cut.passfrac_col;
      c.passfrac_col(:,i)=repmat(v,[1,length(i)]);
    end
  end
end
% cm=cm&c.passfrac_col; % add into overall mask... this does nothing

% apply the channel fraction cut
if(isfield(cut,'passfrac_chan'))
  c.passfrac_chan=true(size(cm));
  for r=unique(p.rx)'
    i=find(p.rx==r);
    v=cp.passfrac_chan(:,r+1)>cut.passfrac_chan;
    c.passfrac_chan(:,i)=repmat(v,[1,length(i)]);
  end
end

return
