function map=make_ebmap(m,map,pure_b,ellrng,purifmat)
% map=make_ebmap(m,map,pure_b,ellrng,purifmat)
%
% Take Q/U to F-plane, convert to E/B, filter to ellrng, and go
% back to image

if ~exist('pure_b','var') || isempty(pure_b)
  pure_b='normal';
end
if ~exist('ellrng','var') || isempty(ellrng)
  ellrng=[50,120];
end
if ~exist('purifmat','var')
  purifmat=[];
end

if ischar(purifmat)
  purifmat=repmat({purifmat},size(map));
end

ad=calc_ad2([m.xdos,m.ydos],[m.nx,m.ny]);
l=ad.u_r*2*pi;

% apply matrix purification to make Q/U maps containing only
% B-modes
if iscell(purifmat)
  % adding .QprojB and .UprojB fields on first iteration causes later entries
  % to gain empty fields; send original through projection
  map_orig = map;
  for i=1:numel(map)
    if ~isempty(purifmat{i})
      [mapE,mapB]=do_projection(m,map_orig(i),purifmat{i});
      map(i).QprojB=mapB.Q; map(i).UprojB=mapB.U;
    end
  end
  clear map_orig
end

% if pol ap mask is not already available make it
if(~isfield(map,'Pw'))
  map=add_masks(m,map);
end

for i=1:numel(map)
  % normalize the weight mask so E/B maps have roughly the same
  % normalization as Q/U
  win=map(i).Pw; win=win./prctile(win(:),90);

  % generate E and B modes exactly as is done for power spectrum estimation
  ft{i}=calc_map_fts(map(i),ad,win,[2],pure_b,false);  
end

% apply simple annular f plane mask
mask=l>ellrng(1)&l<ellrng(2);
for j=1:length(map)
  ft{j}.Ef=ft{j}.E;
  ft{j}.Ef(~mask)=0;
  
  ft{j}.Bf=ft{j}.B;
  ft{j}.Bf(~mask)=0;
end

% go back to image plane
for j=1:length(map)
  map(j).E=f2i(ad,ft{j}.Ef);
  map(j).B=f2i(ad,ft{j}.Bf);
end

% make "vector components" of E&B maps
[u,v]=meshgrid(ad.u_val{1},ad.u_val{2});
for j=1:length(map)
  [Qf,Uf]=eb2qu(u,v,ft{j}.Ef,zeros(size(ft{j}.Bf)));
  map(j).EQ=f2i(ad,Qf);
  map(j).EU=f2i(ad,Uf);
  [Qf,Uf]=eb2qu(u,v,zeros(size(ft{j}.Bf)),ft{j}.Bf);
  map(j).BQ=f2i(ad,Qf);
  map(j).BU=f2i(ad,Uf);
end

% apply mask to make extent same as Q/U
for j=1:length(map)
  mask=isnan(map(j).Pw);

  if(all(size(mask)==size(map(j).E)))
    % only do the re-mask when an upsample has not occurred
    map(j).E(mask)=NaN;
    map(j).B(mask)=NaN;
    
    map(j).EQ(mask)=NaN; map(j).EU(mask)=NaN;
    map(j).BQ(mask)=NaN; map(j).BU(mask)=NaN;
  end
end

return
