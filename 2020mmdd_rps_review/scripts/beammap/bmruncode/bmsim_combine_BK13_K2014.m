function bmsim_combine_BK13_K2014()
%  bmsim_combine_BK13_K2014()
%  this makes the beam map sims for the 1459/????_aabd_ maps (BK14)

% Stuff from real maps combination
%queue='serial_requeue,itc_cluster';
%deprojsreal={'1100','1102'};
% this the latest addition to the set: K2014
%nbase1 = 1351
%daughter1 = 'd'
% this is the BK13, the previous map deepest combination
%nbaseP = 1456
% this is the BK14 = BK13 + K2014
%nbase2 = 1459
    
%{ 
% This was done for 2014 beams analysis
deprojsreal={'1102'};
deprojsreal2={'110200'};
jacks = get_default_coaddopt();
jacks = '0';
nbase1 = 2544
daughter1 = 'a'
nbaseP = 2549
daughterP = 'aab'
nbase2 = 2549
daughter2 = 'aabd'
%  real
for deproj = deprojsreal
  for j=jacks
    
    mBK = [num2str(nbaseP),'/0002_',daughterP,'_filtp3_weight3_gs_dp',deproj{:},'_jack',j,'.mat']
    mK  = [num2str(nbase1),'/0001_',daughter1,'_filtp3_weight3_gs_dp',deprojsreal2{:},'_jack',j,'1.mat']
    mS  = [num2str(nbase2),'/0001_',daughter2,'_filtp3_weight3_gs_dp',deproj{:},'_jack',j,'.mat']
    %if ~exist(['maps/',mS],'file')
      %farmit('farmfiles/coaddcoadd/','load_and_combine(mBK,mK,mS)','func',{@load_and_combine,@minimize_coaddopt},'var',{'mBK','mK','mS'},'queue',queue,'mem',10000,'maxtime',15,'submit',0)    
      load_and_combine(mBK,mK,mS);
    %end
  end
end
%}    
    
% 2015 analysis
% 2590_aab = BK13
% 2562_a = K2014
% 2591_aabd = BK14
jack = '0'; %123456789abcde'; 
deproj = {'0000','1000','0100','1100','1110','1101','1111','1102','1112'};

% Now we combine 1) composite maps, 2) split maps, 3) floor maps
% for r = 1.2, 2, 4 degrees
% However, the BICEP2 split maps are kind of screwed up, so here I will
% use K13 instead of BK13 splits.

% composite 1.2, composite 2, composite 4, composite 6, composite 8
% split 1.2, split 2, split 4, split 6, split 8
% floor 1.2, floor 2, floor 4, floor 6, floor 8
BK = {'3630/0001_aab','3631/0001_aab','3632/0001_aab','3633/0001_aab','3634/0001_aab',...
      '3640/0001_ab','3641/0001_ab','3642/0001_ab','3643/0001_ab','3644/0001_ab',...
      '3650/0001_aab','3651/0001_aab','3652/0001_aab','3653/0001_aab','3654/0001_aab'};
K = {'3630/0001_d','3631/0001_d','3632/0001_d','3633/0001_d','3634/0001_d',...
     '3640/0001_d','3641/0001_d','3642/0001_d','3643/0001_d','3644/0001_d',...
     '3650/0001_d','3651/0001_d','3652/0001_d','3653/0001_d','3654/0001_d'};
S = {'3630/0001_aabd','3631/0001_aabd','3632/0001_aabd','3633/0001_aabd','3634/0001_aabd',...
     '3640/0001_aabd','3641/0001_aabd','3642/0001_aabd','3643/0001_aabd','3644/0001_aabd',...
     '3650/0001_aabd','3651/0001_aabd','3652/0001_aabd','3653/0001_aabd','3654/0001_aabd'};

for jj = 1:length(jack)
  for dd = 1:length(deproj)
    for mm = 1:length(BK)
      mBK = ([BK{mm} '_filtp3_weight3_gs_dp' deproj{dd} '_jack' jack(jj) '.mat']);
      mK =  ([K{mm} '_filtp3_weight3_gs_dp' deproj{dd} '_jack' jack(jj) '1.mat']);
      mS =  ([S{mm} '_filtp3_weight3_gs_dp' deproj{dd} '_jack' jack(jj) '.mat']);
      load_and_combine(mBK,mK,mS);
    end
  end
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function load_and_combine(mBK,mK,mS)
%  this needs to be adjusted year to year
  mBK
  mK
  mS

  mapBK  = load(['maps/',mBK]);
  mapK   = load(['maps/',mK]);

  % make sure you got the right ukpervolt for that year and apply it:
  mapK.coaddopt.ukpv_applied = get_ukpervolt('2014');
  mapK.ac=cal_coadd_ac(mapK.ac,mapK.coaddopt.ukpv_applied);
  mapK.ac=coadd_ac_overfreq(mapK.ac,mapK.coaddopt);
  
  mapK.ac = rmfield(mapK.ac,'wsd');  
  mapK.ac = rmfield(mapK.ac,'wcd');
  actemp=struct_merge(mapK.ac,mapBK.ac);
  ac(1,1) = actemp(1);
  % this is handpicking the right freqs to make the deepest map
  % at each frequency:
  ac(2,1)=coadd_ac_overrx(actemp([2,3]));

  mapBK.coaddopt = minimize_coaddopt(mapBK.coaddopt);
  mapK.coaddopt  = minimize_coaddopt(mapK.coaddopt);
  
  mapK.coaddopt.mapname = mK;
  mapBK.coaddopt.mapname = mBK;
  
  % Above line kills all information in coaddopt since it's a 3-element cell
  % Below allows us to make the aps!
  mapBK.coaddopt.coaddtype = 0;
  mapBK.coaddopt.jacktype = 0;
  mapBK.coaddopt.ukpv_applied = [1]; % field needs to exist so we don't
                                     % cal during reduc_makeaps
  
  coaddopt={mapBK.coaddopt;mapK.coaddopt};
  m = mapBK.m;
  
  saveandtest(['maps/',mS],'ac','coaddopt','m','-v7.3');
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
