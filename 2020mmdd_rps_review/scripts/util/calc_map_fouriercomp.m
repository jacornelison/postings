function map=calc_map_fouriercomp(map,ad,pure_b)
% map=calc_map_fouriercomp(map,ad,pure_b)
%
% This is legacy function which wraps new calc_map_fts so it looks
% like old calc_map_fouriercomp

if(~exist('pure_b','var'))
  pure_b='normal';
end

% calculate fourier grids
[u,v]=meshgrid(ad.u_val{1},ad.u_val{2});

% do ft's
for i=1:numel(map)
  
  if isfield(map,'T')
    ft=calc_map_fts(map(i),ad,1./map(i).Tvar,[1]);
    map(i).Tft=ft.T;
    % are these fields actually used anywhere?...
    map(i).Tw =  nansum(rvec(1./map(i).Tvar));
  end

  if isfield(map,'Tsum')
    ft=calc_map_fts(map(i),ad,1./map(i).Tvar,[3]);
    map(i).Tsumft=ft.Tsum; map(i).Tdifft=ft.Tdif;
  end
  
  if(isfield(map,'Q'))
    % average the Q/U var maps
    map(i).Qvar=(map(i).Qvar+map(i).Uvar)/2;
    map(i).Uvar=map(i).Qvar;

    ft=calc_map_fts(map(i),ad,1./map(i).Qvar,[2],pure_b);
    map(i).Qft=ft.Q; map(i).Uft=ft.U;
    map(i).Eft=ft.E; map(i).Bft=ft.B;

    map(i).Ew = nansum(rvec(1./map(i).Qvar));
    map(i).Bw = map(i).Ew;

  end
end

return
