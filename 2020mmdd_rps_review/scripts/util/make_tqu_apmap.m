function map=make_tqu_apmap(m,map,pure_b,ellrng,purifmat)
% map=make_tqu_apmap(m,map,pure_b,ellrng,purifmat)
%
% Take Q/U to F-plane, filter to ellrng, and go back to image
% code derived from make_ebmap.m

if(~exist('pure_b','var'))
  pure_b=[];
end
if(~exist('ellrng','var'))
  ellrng=[];
end
if(~exist('purifmat','var'))
  purifmat=[];
end

if(isempty(pure_b))
  pure_b='normal';
end
if(isempty(ellrng))
  ellrng=[50,120];
end


ad=calc_ad2([m.xdos,m.ydos],[m.nx,m.ny]);
l=ad.u_r*2*pi;

% average the Q/U var maps
for i=1:numel(map)
  map(i).Qvar=(map(i).Qvar+map(i).Uvar)/2;
  map(i).Uvar=map(i).Qvar;
end

for i=1:numel(map)
  
  % make weight mask
  win=1./map(i).Qvar; win=win./prctile(win(:),90);

  % do ft's using same code as reduc_makeaps
  ft(i)=calc_map_fts(map(i),ad,win,[1,2],pure_b,false);
end

% simple f plane mask
mask=l>ellrng(1)&l<ellrng(2);
for j=1:length(map)
  ft(j).Tf=ft(j).T;
  ft(j).Tf(~mask)=0;
  
  ft(j).Qf=ft(j).Q;
  ft(j).Qf(~mask)=0;
  
  ft(j).Uf=ft(j).U;
  ft(j).Uf(~mask)=0;
end

% go back to image plane
for j=1:length(map)
  map(j).Tap=f2i(ad,ft(j).Tf);
  map(j).Qap=f2i(ad,ft(j).Qf);
  map(j).Uap=f2i(ad,ft(j).Uf);
end

% apply mask to make extent same as Q/U
for j=1:length(map)
  mask=isnan(map(j).Qvar) | isnan(map(j).Uvar);
  map(j).Tap(mask)=NaN;
  map(j).Qap(mask)=NaN;
  map(j).Uap(mask)=NaN;
end

return
