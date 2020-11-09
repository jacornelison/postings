function r=rescale_sim_dust(r,bsnn,bpwf,mb,ms,mc,md,alpha)
% r=rescale_sim(r,bpwf)
%
% construct an r.sim field in r containing realizations of an
% arbitrary cosmology+dust model
%
% r - reduc_final output used for db and rwf fields and as
% container for output
%
% bsnn - auto and cross aps from reduc_makeaps for b x s x n1 x n2
% x ... nN where b is original r=0.1 sims (being repurposed as
% dust), s is original lensed-LCDM sims (being repurposed as
% lensed-LCDM+r) and n1 and n2 are noise in two respective
% channels. These terms are needed to construct the realizations
%
% bpwf - the bpwf's for each element of r
%
% mb - the model used to make the standard r=0.1 sims
% ms - the model corresponding to the lensed-LCDM sims
%
% mc - the desired output "cosmology" component - presumably
% lensed-LCDM+r
% md - the desired output "dust" spectral shape - presumably a simple
% power law normalized to unity at l=80 or similar
% alpha - the dust power at ell=80 in CMB units in each
% frequency channel
%
% e.g:
% load final/0751x1613/real_a_filtp3_weight3_gs_dp1102_jack0_real_a_filtp3_weight3_gs_dp1100_jack0_pureB_matrix_directbpwf
% load aps/0751x1613/xxx4_a_filtp3_weight3_gs_dp1100_jack0_xxx5_a_filtp3_weight3_gs_dp1100_jack0_xxx6_a_filtp3_weight3_gs_dp1100_jack0_xxx6_a_filtp3_weight3_gs_dp1100_jack0_matrix.mat
%               
% ms=load_cmbfast('input_maps/official_cl/camb_planck2013_r0_lensing.fits');
% mb=load_cmbfast('input_maps/official_cl/camb_planck2013_r0p1_noE.fits');
%
% pl=powlaw([1,-0.42],ms.l); pl=pl./pl(find(ms.l==80));
% md=ms; md.Cs_l=zeros(size(md.Cs_l)); md.Cs_l(:,4)=pl; md.Cs_l(:,3)=2*pl;
% mc=ms;
%
% dustrelpow=[1,0,0,0,[0.430,0.901,3.096,24.63].^2];
%
% rn=rescale_sim_dust(r,aps,bpwf,mb,ms,mc,md,alpha);

% number of underlying maps (freq channels)
n=-0.5+sqrt(2*length(r)+0.25);

% strip off the alt crosses - I think they could in principle be
% dealt with correctly but as things stand even EE doesn't work
% because the r=0.1 sims have zero E-mode. We are primarily
% interested in BB...
for i=1:length(bsnn)
  bsnn(i).Cs_l=bsnn(i).Cs_l(:,1:6,:);
end
  
% loop over all spectra
for i=1:n
  for j=i:n

    % find the index of this auto or cross spectrum
    k=aps_getxind(r,i,j);
    
    % calc the expv for the r=0.1 model
    rtmp=calc_expvals(r(k),mb,bpwf(k));
    expvb=rtmp.expv;
    % calc the expv for the lensed-LCDM model
    rtmp=calc_expvals(r(k),ms,bpwf(k));
    expvs=rtmp.expv;
    
    % calc the expv for the new cosmology model
    rtmp=calc_expvals(r(k),mc,bpwf(k));
    expvc=rtmp.expv;
    % calc the expv for the new "dust" model
    rtmp=calc_expvals(r(k),md,bpwf(k));
    expvd=rtmp.expv;
    
    % per bandpower scale factors are just the ratio
    a1=expvc./expvs;
    a2=expvd./expvb;
    
    % fix when expv is nan
    a1(isnan(a1))=1; a2(isnan(a2))=1;
  
    % get rid of alt cross when present
    a1=a1(:,1:6); a2=a2(:,1:6);
    
    % expand c over realz
    nrlz=size(bsnn(1).Cs_l,3);
    a1=repmat(a1,[1,1,nrlz]);
    a2=repmat(a2,[1,1,nrlz]);
    
    % add up the terms to make scaled sims
    % bsnn contains bb, ss, n1n1, n2n2, n3n3, n4n4, ...
    % c=sqrt(a1)*s, d=sqrt(a2)*b
    % terms are (c+alpha_1*d+n_1)(c+alpha_2*d+n_2)
    t1=a1.*bsnn(2).Cs_l;
    t2=sqrt(alpha(j)).*sqrt(a1).*sqrt(a2).*bsnn(aps_getxind(bsnn,1,2)).Cs_l;
    t3=sqrt(a1).*bsnn(aps_getxind(bsnn,2,2+j)).Cs_l;

    t4=sqrt(alpha(i)).*sqrt(a1).*sqrt(a2).*bsnn(aps_getxind(bsnn,1,2)).Cs_l;
    t5=sqrt(alpha(i)).*sqrt(alpha(j)).*a2.*bsnn(1).Cs_l;
    t6=sqrt(alpha(i)).*sqrt(a1).*bsnn(aps_getxind(bsnn,1,2+j)).Cs_l;

    t7=sqrt(a1).*bsnn(aps_getxind(bsnn,2,2+i)).Cs_l;
    t8=sqrt(alpha(j)).*sqrt(a2).*bsnn(aps_getxind(bsnn,1,2+i)).Cs_l;
    t9=bsnn(aps_getxind(bsnn,2+i,2+j)).Cs_l;

    sim=t1+t2+t3+t4+t5+t6+t7+t8+t9;
    %sim=t5;
    
    % get debias and supfac
    db=r(k).db; rwf=r(k).rwf;
    
    % again get rid of alt cross
    db=db(:,1:6); rwf=rwf(:,1:6);
    
    % apply
    r(k).sim=sim.*repmat(rwf,[1,1,nrlz])-repmat(db,[1,1,nrlz]);
    %r(k).sim=sim;
  end
end

return
