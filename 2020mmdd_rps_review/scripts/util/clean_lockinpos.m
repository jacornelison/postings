function d=clean_lockinpos(d)
% d=clean_lockinpos(d)
%
% lockin.azPos etc contain glitch values in 2005,6 data due to word
% tearing in PMAC readout. We can recognize these due to an apparent
% step in the timestream which goes back after one sample.
% Find these values and replace them with the previous value.
%
% This algorithm does not get 100% because there are a very few "double
% glitches" at least in asPos where the value changes for 2 samples.

% find suspiciously large +ve and -ve steps
indp=find(diff(d.lockin.azPos)>+0.04);
indn=find(diff(d.lockin.azPos)<-0.04);
% find +ve step followed by -ve and -ve followed by +ve
ind=[intersect(indp,indn-1);intersect(indn,indp-1)];
% copy values from pre-glitch location into glitch
d.lockin.azPos(ind+1)=d.lockin.azPos(ind);

% repeat for el and dk

indp=find(diff(d.lockin.elPos)>+0.04);
indn=find(diff(d.lockin.elPos)<-0.04);
ind=[intersect(indp,indn-1);intersect(indn,indp-1)];
d.lockin.elPos(ind+1)=d.lockin.elPos(ind);

indp=find(diff(d.lockin.dkPos)>+0.04);
indn=find(diff(d.lockin.dkPos)<-0.04);
ind=[intersect(indp,indn-1);intersect(indn,indp-1)];
d.lockin.dkPos(ind+1)=d.lockin.dkPos(ind);

return
