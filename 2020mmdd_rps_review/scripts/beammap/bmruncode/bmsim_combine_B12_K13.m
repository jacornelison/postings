function bmsim_combine_B12_K13()
%  bmsim_combine_B12_K13()
%  this makes the beam map sims for the 1456/????_aab_ maps (BK13) 

%{
% CLW's B2 map corresponding to sim 0751
% and coadded Keck map corresponding to 1351_ab
for j='0'
  mapB12 = (['2063/0002_b_filtp3_weight3_gs_dp1102_jack',j,'.mat']);
  mapK13 = (['2086/0002_ab_filtp3_weight3_gs_dp1102_jack',j,'.mat']);
  mapBK13= (['2549/0002_aab_filtp3_weight3_gs_dp1102_jack',j,'.mat']);
  load_and_combine(mapB12,mapK13,mapBK13)
  %farmit('farmfiles/coaddcoadd/','load_and_combine(mapK13,mapB12,mapK13)','func',{@load_and_combine,@minimize_coaddopt},'var',{'mapK13','mapB12','mapK13'},'queue',queue,'mem',10000,'maxtime',15,'submit',0)    
end
%}    

% KSK's B2 map corresponding to sim 0751
% and coadded Keck map corresponding to 1351_ab
jack = '0'; %123456789abcde';
deproj = {'0000','1000','0100','1100','1110','1101','1111','1102','1112'};
% Pre- 201612: B12 = 2568/0002_a, K13 = 2589/0001_ab, BK13 =
% 2590/0001_aab

% Now we combine 1) composite maps, 2) split maps, 3) floor maps
% for r = 1.2, 2, 4 degrees
% However, the BICEP2 split maps are kind of screwed up, so we will forgo
% them for now.  In the BK13 + K2014, use K13 instead of BK13 splits.

% composite 1.2, composite 2, composite 4, composite 6, composite 8
% floor 1.2, floor 2, floor 4, floor 6, floor 8
B12 = {'3625/0001_a','3626/0001_a','3627/0001_a','3628/0001_a','3629/0001_a',...
       '3645/0001_a','3646/0001_a','3647/0001_a','3648/0001_a','3649/0001_a'};
K13 = {'3630/0001_ab','3631/0001_ab','3632/0001_ab','3633/0001_ab','3634/0001_ab',...
       '3650/0001_ab','3651/0001_ab','3652/0001_ab','3653/0001_ab','3654/0001_ab'};
BK13 = {'3630/0001_aab','3631/0001_aab','3632/0001_aab','3633/0001_aab','3634/0001_aab',...
        '3650/0001_aab','3651/0001_aab','3652/0001_aab','3653/0001_aab','3654/0001_aab'};

for jj = 1:length(jack)
  for dd = 1:length(deproj)
    for mm = 1:length(B12)
      mapB12 = ([B12{mm} '_filtp3_weight3_gs_dp' deproj{dd} '_jack',jack(jj),'.mat']);
      mapK13 = ([K13{mm} '_filtp3_weight3_gs_dp' deproj{dd} '_jack',jack(jj),'.mat']);
      mapBK13= ([BK13{mm} '_filtp3_weight3_gs_dp' deproj{dd} '_jack',jack(jj),'.mat']);
      load_and_combine(mapB12,mapK13,mapBK13)
    end
  end
end

return


function load_and_combine(mapB12,mapK13,mapBK13)
%  this needs to be adjusted year to year
  mapB12
  mapK13
  mapBK13
  
  mB12  = load(['maps/',mapB12]);
  mK13  = load(['maps/',mapK13]);
    
  % make sure you got the right ukpervolt for that year here:
  mB12.coaddopt.ukpv_applied = 3150;
  mB12.ac=cal_coadd_ac(mB12.ac,mB12.coaddopt.ukpv_applied);
  
  % for K13 nothing needs to be done...
  %
  
  ac = struct_merge(mB12.ac,mK13.ac);
  ac = coadd_ac_overrx(ac);
  
  mB12.coaddopt = minimize_coaddopt(mB12.coaddopt);
  
  mB12.coaddopt.mapname = mapB12;
  
  coaddopt={mB12.coaddopt,mK13.coaddopt{:}};
  m = mK13.m;
  
  saveandtest(['maps/',mapBK13],'ac','coaddopt','m','-v7.3');
return

function coaddopt = minimize_coaddopt(coaddopt)
  try coaddopt = rmfield(coaddopt,'b'); end
  try coaddopt = rmfield(coaddopt,'bi'); end
  try coaddopt = rmfield(coaddopt,'bw'); end
  try coaddopt = rmfield(coaddopt,'hsmax'); end
  try coaddopt = rmfield(coaddopt,'whist'); end
  try coaddopt = rmfield(coaddopt,'devhist'); end
  try coaddopt = rmfield(coaddopt,'traj'); end
  try coaddopt = rmfield(coaddopt,'c'); end
  try coaddopt = rmfield(coaddopt,'jackmask'); end
return
