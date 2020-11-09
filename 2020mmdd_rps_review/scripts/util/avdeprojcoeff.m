function b=avdeprojcoeff(b,bw)
% b=avdeprojcoeff(b,bw)
%
% average coefficients over scan dir and time

sz=size(b);
nch=sz(1);
nsc=sz(2);
nb=sz(3);
nph=sz(4);

bw=reshape(bw,[nch,nsc,1,nph]);
bw=repmat(bw,[1,1,nb,1]);

% average over left/right
b=nansum(b.*bw,2)./nansum(bw,2);
bw=nansum(bw,2);

% average over time
b=squeeze(nansum(b.*bw,4)./nansum(bw,4));

return
