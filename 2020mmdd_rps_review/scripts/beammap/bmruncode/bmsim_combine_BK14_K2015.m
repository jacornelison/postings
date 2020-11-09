function bmsim_combine_BK14_K2015()
%  bmsim_combine_BK14_K2015()
%  this makes the beam map sims for the 1459/????_aabde_ maps (BK15)

%{    
queue='serial_requeue,itc_cluster';

deprojsreal={'1100','1102'};
jacks = '0345';

% this the latest addition to the set: K2014
nbase1 = 1351
daughter1 = 'e'

% this is BK14, the previous map deepest combination
nbaseP = 1459
daughterP = 'aabd'

% this will be BK15 = BK14 + K2015
nbase2 = 1459
daughter2 = 'aabde'
%}

% 2015 analysis
% 2591_aabd = BK14
% 2559_a = K2015
% 2592_aabde = BK15
jack = '0'; 
deproj = {'0000','1000','0100','1100','1110','1101','1111','1102','1112'};
    
% Now we combine 1) composite maps, 2) split maps, 3) floor maps
% for r = 1.2, 2, 4, 6, 8 degrees

% composite 1.2, composite 2, composite 4, composite 6, composite 8
% split 1.2, split 2, split 4, split 6, split 8
% floor 1.2, floor 2, floor 4, floor 6, floor 8
BK = {'3630/0001_aabd','3631/0001_aabd','3632/0001_aabd','3633/0001_aabd','3634/0001_aabd',...
      '3640/0001_aabd','3641/0001_aabd','3642/0001_aabd','3643/0001_aabd','3644/0001_aabd',...
      '3650/0001_aabd','3651/0001_aabd','3652/0001_aabd','3653/0001_aabd','3654/0001_aabd'};
K = {'3630/0001_e','3631/0001_e','3632/0001_e','3633/0001_e','3634/0001_e',...
     '3640/0001_e','3641/0001_e','3642/0001_e','3643/0001_e','3644/0001_e',...
     '3650/0001_e','3651/0001_e','3652/0001_e','3653/0001_e','3654/0001_e'};
S = {'3630/0001_aabde','3631/0001_aabde','3632/0001_aabde','3633/0001_aabde','3634/0001_aabde',...
     '3640/0001_aabde','3641/0001_aabde','3642/0001_aabde','3643/0001_aabde','3644/0001_aabde',...
     '3650/0001_aabde','3651/0001_aabde','3652/0001_aabde','3653/0001_aabde','3654/0001_aabde'};

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
  mapK.coaddopt.ukpv_applied = get_ukpervolt('2015');
  mapK.ac=cal_coadd_ac(mapK.ac,mapK.coaddopt.ukpv_applied);
  mapK.ac=coadd_ac_overfreq(mapK.ac,mapK.coaddopt);
  
  mapK.ac = rmfield(mapK.ac,'wsd');  
  mapK.ac = rmfield(mapK.ac,'wcd');
  % Turn ac(3,1)+ac(2,1) into ac(5,1) by concatenating
  actemp=struct_merge(mapK.ac,mapBK.ac);
  % This is handpicking the right freqs to make the deepest map at each
  % frequency:
  ac(1,:) = coadd_ac_overrx(actemp([1,4],:)); % 100GHz
  ac(2,:) = coadd_ac_overrx(actemp([2,5],:)); % 150GHz
  ac(3,:) = actemp(3,:); % 220GHz

  % For BK=1459_aabd, already done when combining BK13 + K2014
  %mapBK.coaddopt = minimize_coaddopt(mapBK.coaddopt);
  %mapBK.coaddopt.mapname = mBK;

  mapK.coaddopt  = minimize_coaddopt(mapK.coaddopt);
  mapK.coaddopt.mapname = mK;
  
  % Above line kills all information in coaddopt since it's a 3-element cell
  % Below allows us to make the aps!
  mapBK.coaddopt.coaddtype = 0;
  mapBK.coaddopt.jacktype = 0;
  mapBK.coaddopt.ukpv_applied = [1]; % field needs to exist so we don't
                                     % cal during reduc_makeaps

  coaddopt= [mapBK.coaddopt];%; mapK.coaddopt];
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

