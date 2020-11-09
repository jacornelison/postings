function caxis_symzero()
% caxis_symzero()
%
% Force the color stretch to be symmetric about zero

ca=caxis;
ca=max(abs(ca));
caxis([-ca,ca]);

return
