function bmsim_combine_K2012_K2013()
%  bmsim_combine_K2012_K2013()
%  this makes the beam map sims for the 1351/????_ab_ maps (K13) 

%{    
% CLW's individual year maps corresponding to
% the real maps 1351 dp1102
for j='0'
  mapK2012 = (['2086/0002_a_filtp3_weight3_gs_dp1102_jack',j,'1.mat']);
  mapK2013 = (['2086/0002_b_filtp3_weight3_gs_dp1102_jack',j,'1.mat']);
  mapK13   = (['2086/0002_ab_filtp3_weight3_gs_dp1102_jack',j,'.mat']);
  load_and_combine(mapK2012,mapK2013,mapK13)
  %farmit('farmfiles/','load_and_combine(mapK2012,mapK2013,mapK13)','func',{@load_and_combine,@minimize_coaddopt},'var',{'mapK2012','mapK2013','mapK13'},'queue',queue,'mem',10000,'maxtime',15,'submit',0)    

end
%}
    
% KSK's individual year maps corresponding to
% the real maps 1351 dp1102
% These were made with "x-year" beam maps incorporating all available
% data
jack = '0'; %123456789abcde';
deproj = {'0000','1000','0100','1100','1110','1101','1111','1102','1112'};
% Pre- 201612: K2012 = 2587/0001_a, K2013 = 2588/0001_a, K13 =
% 2589/0001_ab

% Now we combine 1) composite maps, 2) split maps, 3) floor maps
% for r = 1.2, 2, 4, 6, 8 degrees
%
% composite 1.2, composite 2, composite 4, composite 6, composite 8
% split 1.2, split 2, split 4, split 6, split 8
% floor 1.2, floor 2, floor 4, floor 6, floor 8
K2012 = {'3630/0001_a','3631/0001_a','3632/0001_a','3633/0001_a','3634/0001_a'...
         '3640/0001_a','3641/0001_a','3642/0001_a','3643/0001_a','3644/0001_a'...
         '3650/0001_a','3651/0001_a','3652/0001_a','3653/0001_a','3654/0001_a'};
K2013 = {'3630/0001_b','3631/0001_b','3632/0001_b','3633/0001_b','3634/0001_b'...
         '3640/0001_b','3641/0001_b','3642/0001_b','3643/0001_b','3644/0001_b'...
         '3650/0001_b','3651/0001_b','3652/0001_b','3653/0001_b','3654/0001_b'};
K13 =   {'3630/0001_ab','3631/0001_ab','3632/0001_ab','3633/0001_ab','3634/0001_ab'...
         '3640/0001_ab','3641/0001_ab','3642/0001_ab','3643/0001_ab','3644/0001_ab'...
         '3650/0001_ab','3651/0001_ab','3652/0001_ab','3653/0001_ab','3654/0001_ab'};

for jj = 1:length(jack)
  for dd = 1:length(deproj)
    for mm = 1:length(K2012)
      mapK2012 = ([K2012{mm} '_filtp3_weight3_gs_dp' deproj{dd} '_jack',jack(jj),'1.mat']);
      mapK2013 = ([K2013{mm} '_filtp3_weight3_gs_dp' deproj{dd} '_jack',jack(jj),'1.mat']);
      mapK13   = ([K13{mm} '_filtp3_weight3_gs_dp' deproj{dd} '_jack',jack(jj),'.mat']);
      load_and_combine(mapK2012,mapK2013,mapK13)
    end
  end
end


%babysitjobs('farmfiles/comb/*','wait5')
      
return


function load_and_combine(mapK2012,mapK2013,mapK13)
%  this needs to be adjusted year to year
  mapK2012
  mapK2013
  mapK13
  
  mK2012  = load(['maps/',mapK2012]);
  mK2013  = load(['maps/',mapK2013]);
  
  % make sure you got the right ukpervolt for that year here:
  mK2012.coaddopt.ukpv_applied = get_ukpervolt('2012');
  mK2012.ac = cal_coadd_ac(mK2012.ac,mK2012.coaddopt.ukpv_applied);
  mK2012.ac = coadd_ac_overfreq(mK2012.ac,mK2012.coaddopt);
  mK2012.coaddopt.coaddtype = 0; % needed to keep reduc_makaps from choking
  
  mK2013.coaddopt.ukpv_applied = get_ukpervolt('2013');
  mK2013.ac = cal_coadd_ac(mK2013.ac,mK2013.coaddopt.ukpv_applied);
  mK2013.ac = coadd_ac_overfreq(mK2013.ac,mK2013.coaddopt);
  mK2013.coaddopt.coaddtype = 0;
  
  ac = struct_merge(mK2012.ac,mK2013.ac);
  ac = coadd_ac_overrx(ac);
  
  mK2012.coaddopt = minimize_coaddopt(mK2012.coaddopt);
  mK2013.coaddopt = minimize_coaddopt(mK2013.coaddopt);  
  
  mK2012.coaddopt.mapname = mapK2012;
  mK2013.coaddopt.mapname = mapK2013;
  
  coaddopt = {mK2012.coaddopt;mK2013.coaddopt};
  m = mK2013.m;
  
  saveandtest(['maps/',mapK13],'ac','coaddopt','m','-v7.3');
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
