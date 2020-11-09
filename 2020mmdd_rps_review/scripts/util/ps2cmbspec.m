function [l,Cs_l,C_l]=ps2cmbspec(u,del_u,C_u)
% [l,Cs_l,C_l]=ps2cmbspec(u,C_u)
%
% Convert output of powspec to Cscript_l form
%
% e.g.: ad=calc_ad(20,256);
%       [l,Cs_l,C_l] = get_cmb_model(0.05,0.35,0.60,65,1,0);
%       cmb=cmb_skysim(ad,l,C_l);
%       [C_u,u]=powspec(ad,i2f(ad,cmb));
%       [lo,Cs_lo]=ps2cmbspec(u,ad.del_u,C_u);
%       plot(l,Cs_l,'r');
%       hold on; plot(lo,Cs_lo,'b'); hold off

if(length(del_u)==1)
  del_u=[del_u,del_u];
end

l=u*2*pi;

C_l=C_u*prod(del_u);
Cs_l=C_l.*l.*(l+1)/(2*pi);

return
