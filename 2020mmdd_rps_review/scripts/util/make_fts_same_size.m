function [ft1,ft2,ad]=make_fts_same_size(ad1,ad2,ft1,ft2)
% [ft1,ft2,ad]=make_fts_same_size(ad1,ad2,ft1,ft2)
%
% When one ft is bigger than the other pad the smaller one
% to get same size
%
% This is utility function used by reduc_makeaps.m and util/make_aps2d.m

if(size(ad2.t_r,1)==size(ad1.t_r,1))
  ad=ad1;
  return
end

if(size(ad2.t_r,1)>size(ad1.t_r,1))
  disp(sprintf('warning: padding ft1 from %dx%d to %dx%d',size(ad1.t_r),size(ad2.t_r)));
  ad=ad2;
  [dum,ft1]=pad_ft([],ft1,size(ad2.t_r));
else
  disp(sprintf('warning: padding ft2 from %dx%d to %dx%d',size(ad2.t_r),size(ad1.t_r)));
  ad=ad1;
  [dum,ft2]=pad_ft([],ft2,size(ad1.t_r));
end

return
