function [lc en t] = get_cals(t1,t2)
% [lc t] = get_cals(t1,t2)
% Accumulate load curve data over the interval defined by t1 and t2
% Returns an n x m structure with elements .pj_av, .pj_end, .rnorm and
% .rdet. n is 2x the number of tags and m is the number of channels

% NOTE: Does not currently accumulate anything other than time and load
% curve data.  Need to expand to include field scans and elnods

tags=list_tags(t1,t2,'has_tod',1);

pj_av=[];
pj_end=[];
rnorm=[];
rdet=[];
t=[];
g1=[]; g2=[]; g3=[]; g4=[];

for ii=1:length(tags)
  load(['data/real/' tags{ii}(1:6) '/' tags{ii} '_calval'])
  t = [t; lc.t];
  pj_av = [pj_av; lc.g(:,:,1)];
  pj_end = [pj_end; lc.g(:,:,2)];
  rnorm = [rnorm; lc.g(:,:,3)];
  rdet = [rdet; lc.g(:,:,4)];
  g1 = [g1; en.g(:,:,1)];
  g2 = [g2; en.g(:,:,2)];
  g3 = [g3; en.g(:,:,3)];
  g4 = [g4; en.g(:,:,4)]; 
end
clear lc en

lc.pj_av = pj_av;
lc.pj_end = pj_end;
lc.rnorm = rnorm;
lc.rdet = rdet;

en.g(:,:,1) = g1;
en.g(:,:,2) = g2;
en.g(:,:,3) = g3;
en.g(:,:,4) = g4;
