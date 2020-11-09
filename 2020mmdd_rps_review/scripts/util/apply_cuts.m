function d=apply_cuts(d,fs,cm)
% d=apply_cuts(d,fs,cm)
%
% nan out any half-scan/channel which is marked false in mask cm
%
% d=apply_cuts(d,fs,cm)

disp('applying cuts...');

% Apply to tod
for ihs=1:length(fs.sf)
  d.mce0.data.fb(fs.sf(ihs):fs.ef(ihs),cm(ihs,:)==0)=NaN;
end

return
