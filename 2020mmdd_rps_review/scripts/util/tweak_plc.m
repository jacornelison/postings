function lc=tweak_plc(d,lc,sampratio)
% lc=tweak_plc(d,lc,sampratio)
%
% Select sample rate=1 and not part of field scans.
% The feature bit is not entirely reliable since it may be
% turned off as soon as the PLC script completes on one MCE.
% Here, expand the PLC region to the end of the PLC script
% on the last MCE to finish.
nsnap=double(d.array.frame.nsnap);
fbits=double(d.array.frame.features);
plc_good=(nsnap==1) & (bitand(fbits,2^14) | fbits==0);
% Expand to fast rate
plc_good=reshape(repmat(reshape(plc_good,1,[]),sampratio,1),[],1);
% Extend each PLC to just before first non-good-plc sample
for ii=1:length(lc.sf)
  j=find(~plc_good(lc.ef(ii):end),1,'first');
  lc.ef(ii)=lc.ef(ii)+j-1;
end

return

