function [map,coaddopt,m,map_n,map_lcdm,map_r,map_d,map_s]=construct_b2planck_sims(expt,rlzn,rlzs,r,d)
% [map,coaddopt,m,map_n,map_lcdm,map_r,map_d]=construct_b2planck_sims(f,Plband,rlzn,rlzs,r,d)
%
% Outputs simulated Planck or B2 map in uK_CMB containing lensed LCDM, CMB tensors,
% dust, and statistical noise:
% 
% Outputs:
%
% map,coaddopt,m   - all maps added together
% map_n            - noise map(s)
% map_lcdm         - lensed lcdm signal only map
% map_r            - tensors map
% map_d            - dust map
% map_s            - total signal (map_lcdm+map_r+map_d)
%
% Inputs:
%
% expt -  determines what noise realizations are used and frequency scaling of dust map
%         'b150'
%         'p100'
%         'p143'
%         'p217'
%         'p353'
% rlzn - noise realization (if vector, returns an array of noise maps and total
%        signal+noise maps)
% rlzs - signal realization, constrols map_lcdm and map_r
% r - tensor/scalar ratio
% d - scaling factor for Planck PSM, sim mode, pmax=17.5%. Set to one to use the
%     straight PSM map from Jamie Tolan

% Get maps
map_lcdm = get_lcdm_map(rlzs); % lensed LCDM
map_r    = get_bmode_map(rlzs); % r=0.1 tensors
[map_d,coaddopt,m] = get_dust_map(); % PSM dust

% Scale tensor map
map_r = cal_coadd_maps(map_r,r/0.1);

% Scale dust map using numbers from Ken Ganga
% http://b2p.planck.fr/index.php/AnalysisLogbook/2014-08-14DustConversion
f=str2num(expt(2:end));
switch f
 case 100
  fsc=0.424885357233;
 case 143
  fsc=0.900287214698;
 case 150
  fsc=1;
 case 217
  fsc=3.12985190361;
 case 353
  fsc=25.1198188496;
end
map_d = cal_coadd_maps(map_d,d*fsc);

% Add signal maps together
map_s = map_lcdm;
map_s = addmaps(map_s,map_r);
map_s = addmaps(map_s,map_d);

% Get noise maps and add to signal map (add systematics map here if necessary)
for k=1:numel(rlzn)
  map_n(k) = get_noise_map(expt,rlzn(k));
  map(k) = addmaps(map_s,map_n(k));
end

% Only return certain fields, necessary for array concatenation
map=stripmap(map);
map_n=stripmap(map_n);
map_lcdm=stripmap(map_lcdm);
map_r=stripmap(map_r);
map_d=stripmap(map_d);
map_s=stripmap(map_s);

return


%%%%%%%%%%%%%%%%%%%%%%%%
function map=addmaps(map1,map2)

map=map1;
map.T=map1.T+map2.T;
map.Q=map1.Q+map2.Q;
map.U=map1.U+map2.U;

return


%%%%%%%%%%%%%%%%%%%%%%%%%%
function map=get_noise_map(expt,rlz)

switch expt(1)
 
 case 'b'
  load(sprintf('/n/bicepfs2/bicep2/pipeline/maps/0751/%03d6_a_filtp3_weight3_gs_dp1100_jack0.mat',rlz));
  map=cal_coadd_maps(make_map(ac,m),3150);
  
 case 'p'
  load(sprintf('/n/bicepfs2/b2planck/pipeline/maps/1613/%03d6_a_filtp3_weight3_gs_dp1100_jack0.mat',rlz));
  fac=[30,44,70,100,143,217,353];
  ind=find(fac==str2num(expt(2:end)));
  map=cal_coadd_maps(make_map(ac(ind),m),3150);

end

return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function map=get_lcdm_map(rlz)

load(sprintf('/n/bicepfs2/bicep2/pipeline/maps/0751/%03d5_a_filtp3_weight3_gs_dp1100_jack0.mat',rlz));
map=cal_coadd_maps(make_map(ac,m),coaddopt.mapopt{1}.simopt.ukpervolt);

return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function map=get_bmode_map(rlz)

load(sprintf('/n/bicepfs2/bicep2/pipeline/maps/0751/%03d4_a_filtp3_weight3_gs_dp1100_jack0.mat',rlz));
map=cal_coadd_maps(make_map(ac,m),coaddopt.mapopt{1}.simopt.ukpervolt);

return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [map,coaddopt,m]=get_dust_map

load('/n/bicepfs2/b2planck/pipeline/maps/psmdust/0012_a_filtp3_weight3_gs_dp1100_jack0.mat');
map=cal_coadd_maps(make_map(ac,m),3150);

return


%%%%%%%%%%%%%%%%%%%%%%%%%%5
function d=stripmap(b)

for k=1:numel(b)

  a=b(k);
  
  c.x_tic=a.x_tic;
  c.y_tic=a.y_tic;
  c.T=a.T;
  c.Tvar=a.Tvar;
  c.Titime=a.Titime;
  c.Q=a.Q;
  c.U=a.U;
  c.Qvar=a.Qvar;
  c.Uvar=a.Uvar;
  c.QUcovar=a.QUcovar;
  c.Pitime=a.Pitime;

  d(k)=c;

end

return
