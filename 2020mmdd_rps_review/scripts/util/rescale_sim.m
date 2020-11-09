function r=rescale_sim(r,bpwf,mo,mn)
% r=rescale_sim(r,bpwf,mo,mn)
%
% rescale the r.sim field (lensed-LCDM plus noise) to an arbitrary
% model using the r.sim, r.noi, and r.sigl (type5 lensed-LCDM only)
% fields
%
% e.g.
% mo=load_cmbfast('input_maps/official_cl/camb_planck2013_r0_lensing.fits');
% mr=load_cmbfast('input_maps/official_cl/camb_planck2013_r0p1_noE.fits');
% mn=mo; mn.Cs_l(:,4)=mo.Cs_l(:,4)+2*mr.Cs_l(:,4);
% rn=rescale_sim(r,bpwf,mo,mn)

for i=1:length(r)

  % calc the expv for the old model
  rtmp=calc_expvals(r(i),mo,bpwf(i));
  expvo=rtmp.expv;
  % calc the expv for the new model
  rtmp=calc_expvals(r(i),mn,bpwf(i));
  expvn=rtmp.expv;
  % per bandpower scale factors are just the ratio
  c=expvn./expvo;
  
  % fix when expv is zero
  c(isnan(c))=1;
  
  % expand c over realz
  c=repmat(c,[1,1,size(r(i).sim,3)]);
  
  % back out the sn term
  db=repmat(r(i).db,[1,1,size(r(i).sim,3)]);
  sn=(r(i).sim+db-r(i).sigl-r(i).noi)/2;
  
  % apply rescale and store back to sim
  r(i).sim=c.*r(i).sigl+r(i).noi+2*sqrt(c).*sn-db;
end

return
