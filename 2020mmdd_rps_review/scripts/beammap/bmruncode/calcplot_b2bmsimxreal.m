function calcplot_b2bmsimxreal(simnums,calcopt)
% function calcplot_b2bmsimxreal(simnums,calcopt)
% Calculate and plot beam map sim aps, maps pagers, and final spectra
%
% INPUTS
%   
%   simnums: struct specifying the sim numbers to calc/plot 
%              simser, simrlz, simdau, realser, realrlz, realdau
%              "sim" refers to beam map sim - one realization
%              "real" refers to whatever is being crossed with "sim"
%              This could be real data (realrlz = 'real') or sims, e.g.
%              realrlz = '???8' for type 8 dust sims
%   calcopt: struct specifying spectrum calculation option and what to do
%              calcaps: 0 (default), 1 - generate aps
%              updateaps: 0 (default) remakes all, 1 only makes missing
%              farmaps: 0 (default) interactively, 1 to farm
%              makesimset: 0 (default), 1 - make simset from realrlz
%              mappager: 0 (default), 1 - generate maps pager
%              plotaps: 0 (default), 1 - generate aps spectrum plots
%              est: 'pureB' (default), 'matrix'
%              autocross: 'auto' (default), 'cross'
%              crosstype: 'real' (default) toonly cross with 1 real rlz
%                         'sim' to cross the real bm sim with all sim rlz
%              simstocross: realizations of sims to cross (only works
%                           with crosstype = 'sim')

% Parse input options
if ~isfield(simnums,'simser')
  simser = '3626';
else
  simser = simnums.simser;
end
if ~isfield(simnums,'simrlz')
  simrlz = '0001';
else
  simrlz = simnums.simrlz;
end
if ~isfield(simnums,'simdau')
  simdau = 'a';
else
  simdau = simnums.simdau;
end
if ~isfield(simnums,'realser')
  realser = '1450';
else
  realser = simnums.realser;
end
if ~isfield(simnums,'realrlz')
  realrlz = 'real';
else
  realrlz = simnums.realrlz;
end
if ~isfield(simnums,'realdau')
  realdau = 'a';
else
  realdau = simnums.realdau;
end
if ~isfield(calcopt,'calcaps')
  calcaps = 0;
else
  calcaps = calcopt.calcaps;
end
if ~isfield(calcopt,'updateaps')
  updateaps = 0;
else
  updateaps = calcopt.updateaps;
end
if ~isfield(calcopt,'farmaps')
  farmaps = 0;
else
  farmaps = calcopt.farmaps;
end
if ~isfield(calcopt,'makesimset')
  makesimset = 0;
else
  makesimset = calcopt.makesimset;
end
if ~isfield(calcopt,'mappager')
  mappager = 0;
else
  mappager = calcopt.mappager;
end
if ~isfield(calcopt,'plotaps')
  plotaps = 0;
else
  plotaps = calcopt.plotaps;
end
if ~isfield(calcopt,'est')
  est = 'pureB';
else
  est = calcopt.est;
end
if ~isfield(calcopt,'autocross')
  autocross = 'auto';
else
  autocross = calcopt.autocross;
end
if ~isfield(calcopt,'crosstype')
  crosstype = 'real';
else
  crosstype = calcopt.crosstype;
end
if ~isfield(calcopt,'simstocross')
  simstocross = 1:499;
else
  simstocross = calcopt.simstocross;
end

% General setup

% BICEP2 matrix estimator
xb12 = ['/n/panlfs2/bicep/bicep2/pipeline/matrixdata/c_t/0704/' ...
        'healpix_red_spectrum_lmax700_beamB2bbns_reob0704_proj.mat']; 

% apsopt template
switch autocross
  case 'auto'
    load('aps/0751/real_a_filtp3_weight3_gs_dp1102_jack0_matrix');
  case 'cross'
    load(['aps/0751x1351/' ...
          'real_a_filtp3_weight3_gs_dp1102_jack0_'...
          'real_a_filtp3_weight3_gs_dp1102_jack01_matrix_overrx']);
end
apsopt = rmfield(apsopt,'purifmatname');
apsopt.makebpwf = 0;

switch est
  case 'matrix'
    apsopt.purifmatname = xb12;
  case 'pureB'
    apsopt.pure_b = 'kendrick';
end

if updateaps
  apsopt.update = 1;
else
  apsopt.update = 0;
end

apsopt.overrideukpervolt = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate aps

if calcaps
  switch autocross
    case 'auto'
      apsopt.ukpervolt = {get_ukpervolt()};
      apsopt.polrot = [0];
      mapn = {[simser '/' simrlz '_' simdau ...
               '_filtp3_weight3_gs_dp*_jack*.mat']};
      if farmaps
        mkdir(['farmfiles/' simser]);
        cmd = 'reduc_makeaps(mapn,apsopt)';
        farmfile = ['farmfiles/' simser '/' simrlz '_' simdau ...
                    '_auto_' est];
        farmit(farmfile,cmd,'var',{'mapn','apsopt'},'mem',40000, ...
               'maxtime',1440);
      else
        reduc_makeaps(mapn,apsopt);
      end
    case 'cross'
      apsopt.ukpervolt = {get_ukpervolt(),get_ukpervolt()};
      switch est
        case 'matrix'
          apsopt.purifmatname = {xb12,xb12};
      end
      apsopt.commonmask = 2;
      
      switch crosstype
        case 'real' 
          apsopt.polrot = [0,-1.1];
          % do all jacks/dps in one go since it's only one rlz
          % Forget about those extra dp digits...
          mapn = {[simser '/' simrlz '_' simdau '_filtp3_weight3_gs_dp????_jack*.mat'],...
                  [simser '/' simrlz '_' simdau],[realser '/' realrlz '_' realdau]};
          if farmaps
            mkdir(['farmfiles/' simser]);
            cmd = 'reduc_makeaps(mapn,apsopt)';
            farmfile = ['farmfiles/' simser '/' simrlz '_' simdau ...
                        '_crossreal_' est];
            farmit(farmfile,cmd,'var',{'mapn','apsopt'},'mem',40000, ...
                   'maxtime',1440);
          else
            reduc_makeaps(mapn,apsopt);
          end
        case 'sim' 
          apsopt.polrot = [0,0];
          type = realrlz(4); % sim type
          % just do dp1100 / jack0 and farm out all rlz
          for ii = 1:length(simstocross)
            mapreplacenum = sprintf('%03i',simstocross(ii)); % sim number
            mapn = {[simser '/' simrlz '_' simdau '_filtp3_weight3_gs_dp1100_jack0.mat'],...
                    [simser '/' simrlz '_' simdau],...
                    [realser '/' mapreplacenum type '_' realdau]};
            if farmaps
              mkdir(['farmfiles/' simser]);
              cmd = 'reduc_makeaps(mapn,apsopt)';
              farmfile = ['farmfiles/' simser '/' simrlz '_' simdau '_x_' ...
                          realser '_' mapreplacenum type '_' realdau '_' est];
              farmit(farmfile,cmd,'var',{'mapn','apsopt'},'mem',40000, ...
                     'maxtime',1440);
            else
              reduc_makeaps(mapn,apsopt);
            end 
          end % simstocross
      end % crosstype

  end
  
end 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Make simset - only for bm sim X sim realizations, dp1100/jack0

if makesimset
  apsname = [simser 'x' realser '/' simrlz '_' simdau ...
             '_filtp3_weight3_gs_dp1100_jack0_???' ...
             realrlz(4) '_' realdau '_filtp3_weight3_gs_dp1100_jack0_' ...
             est '_cm_overrx.mat'];
  reduc_makesimset(apsname);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Make maps pager

mapn = {[simser '/' simrlz '_' simdau '_filtp3_weight3_gs_dp*_jack*.mat']};
apsopt.mapname = {simser};
apsopt.polrot = 0;
apsopt.ukpervolt = {get_ukpervolt()};

if mappager  % free, then fixed color scale
  reduc_plotcomap_pager(mapn,apsopt,0,[],[],[],[],[],[],0,[],1); 
  reduc_plotcomap_pager(mapn,apsopt,1,[],[],[],[],[],[],0,[],1);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot aps

nbase = simser;
nbase2 = [simser 'x' realser];

% Deprojections - all made for auto/cross
depj{1} = '0000'; depj{2} = '1000'; depj{3} = '0100'; depj{4} = '1100'; 
depj{5} = '1110'; depj{6} = '1101'; depj{7} = '1111'; depj{8} = '1102';
depj{9} = '1112';

% Same for jacks
jack{1} = '0'; jack{2} = '1'; jack{3} = '2'; jack{4} = '3';
jack{5} = '4'; jack{6} = '5'; jack{7} = '6'; jack{8} = '7';
jack{9} = '8'; jack{10} = '9'; jack{11} = 'a'; jack{12} = 'b';
jack{13} = 'c'; jack{14} = 'd'; jack{15} = 'e'; 

lmax{1} = '200'; lmax{2} = '500';

% Get the beam map sim algorithmic floor from noiseless constructed
% Planck sky dp1111 beam map sim
floor = load(['/n/bicepfs1/users/clwong/aps/2027/'...
              '0002_a_filtp3_weight3_gs_dp1111_jack0_pureB.mat']);
% Supfac correct it
flc = load(['final/1450x1350/real_a_filtp3_weight3_gs_dp1111_jack0_'...
            'real_a_filtp3_weight3_gs_dp1111_jack01_pureB_overrx.mat']);
floor.aps.Cs_l = floor.aps.Cs_l.*(flc.r(1).rwf);
% So of course this must be removed from real aps AFTER supfac correction

if plotaps

  set(0,'defaultlinelinewidth',1.0);
  set(0,'DefaultAxesFontSize',14);
  
  figure(1); clf; 
  setwinsize(gcf,1000,600);
  set(gcf,'Visible','off');
  
  for dd = 1:length(depj)

    % Load up the appropriate bias calculated by Chris
    % This needs to be removed AFTER supfac correction
    % We only have 0000, 1100, 1101, 1111 on disk...
    switch depj{dd}
      case {'1000','0100','1110'}
        bias = load(['beammaps/aux_products/bias/' ...
                     'beammap_bias_and_sigma_dp1100']);
      case {'1102'}
        bias = load(['beammaps/aux_products/bias/' ...
                     'beammap_bias_and_sigma_dp1101']);
      otherwise
        bias = load(['beammaps/aux_products/bias/' ...
                     'beammap_bias_and_sigma_dp' depj{dd}]);
    end
    
    % Load up the appropriate supfac - in principle this should have
    % been done for all jacks too (so it would go below), but we only
    % have jack0, so to save time load it here
    switch est
      case 'pureB'
        % See paper_plots_b2_respap14, subfunc make_beammapsimsplot
        final = load(['final/1450x1350/real_a_filtp3_weight3_gs_dp'...
                      depj{dd} '_jack0_real_a_filtp3_weight3_gs_dp'...
                      depj{dd} '_jack01_pureB_overrx.mat']);
      case 'matrix'
        % See paper_plots_b2_respap14, subfunc get_b2xb2
        % This supfac was only made for dp1100/1102!  Just use 1102 
        final = load(['final/0751/real_a_filtp3_weight3_gs_dp1102_'...
                      'jack0_matrix_directbpwf_rbc.mat']);
    end

    for jj = 1:length(jack)
        
      % Load in aps, apply supfac, prepare for plotting fuction
      switch autocross
        case 'auto'
          aps = load(['aps/' simser '/' simrlz '_' simdau ...
                      '_filtp3_weight3_gs_dp' depj{dd} '_jack' ...
                      jack{jj} '_' est '.mat']); 
          % pureB supfac loaded above has 3 fields (B2XKeck), choose first
          % Matrix only has 1 so final.supfac(1) is fine
          aps = apply_supfac(aps.aps,final.supfac(1));
          leg = nbase;
        case 'cross'
          aps = load(['aps/' simser 'x' realser '/' simrlz '_' simdau ...
                      '_filtp3_weight3_gs_dp' depj{dd} '_jack' ...
                      jack{jj} '_' realrlz '_' realdau ...
                      '_filtp3_weight3_gs_dp' depj{dd} '_jack' ...
                      jack{jj} '_' est '_cm_overrx.mat']);
          % aps(3) contains alt crosses which cause apply_supfac to fail
          % zero them out
          aps.aps(3).Cs_l = aps.aps(3).Cs_l(:,1:6);
          % pureB supfac loaded above is B2, Keck, B2xKeck - so for our
          % BM sim result we actually want B2, B2, B2 for all fields!
          % For matrix there is only 1, and this extends final.supfac to
          % auto and cross!
          final.supfac(2) = final.supfac(1);
          final.supfac(3) = final.supfac(1);
          aps = apply_supfac(aps.aps,final.supfac); % Apply to all 3
          % Do we want to plot the real spectra too?  If not, reduce
          % dimensionality of aps 
          % [sim auto, real auto, sim x real] ->
          % [sim auto, sim x real]
          aps = aps([1 3]); 
          leg = {nbase,nbase2};
      end

      % Apply beam map floor, bias correction, and calculate rho - 
      % Is this appropriate for anything but jack0?  
      % Only apply to first element of aps (the auto) since the second
      % will (always?) be a cross
      switch jack{jj}
        case {'0'}
          %aps(1).Cs_l(:,4) = aps(1).Cs_l(:,4) - bias.bias;
          %aps(1).Cs_l(:,4) = aps(1).Cs_l(:,4) - floor.aps.Cs_l(:,4);
          rho = calc_rho(aps(1).Cs_l(:,4),final);
        otherwise
          rho = [];
      end
    
      % Loop over plot ranges
      for ll = 1:length(lmax)
      
        % Plot all 6 spectra
        subplot(2,3,1)
        plot_bmspec(aps,leg,'TT',lmax{ll},'150',jack{jj},depj{dd},[])
        subplot(2,3,2)
        plot_bmspec(aps,leg,'TE',lmax{ll},'150',jack{jj},depj{dd},[])
        subplot(2,3,3)
        plot_bmspec(aps,leg,'EE',lmax{ll},'150',jack{jj},depj{dd},[])
        subplot(2,3,4)
        plot_bmspec(aps,leg,'BB',lmax{ll},'150',jack{jj},depj{dd},rho)
        subplot(2,3,5)
        plot_bmspec(aps,leg,'TB',lmax{ll},'150',jack{jj},depj{dd},[])
        subplot(2,3,6)
        plot_bmspec(aps,leg,'EB',lmax{ll},'150',jack{jj},depj{dd},[])
        

        % Even if we just have auto, save in subfolder nbase x real cause
        % we'll populate  it eventually...
        mkdir(['plots/' nbase2]);
        savename = ['plots/' nbase2 '/' simrlz '_' simdau '_filtp3_weight3' ...
                    '_gs_dp' depj{dd} '_jack' jack{jj} '_150_' est '_' lmax{ll} ...
                    '.eps'];
        print('-depsc2',savename);
        clf()
        
      end % plot range
    end % jack
  end % depj

end

return % Main function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function rho = calc_rho(bb,final)
% Calculate rho statistic - specific for B2, could probably generalize
% Depends on the final file for r=0.1 sim means and std
% aps is already supfac-corrected

rhobins = 2:6;

% This assumes that the FIRST field in the final struct is B2 auto
Brp1 = mean(final.r(1).sig(:,4,:),3);
Bstd = std(final.r(1).sim(:,4,:),[],3);
w = Brp1./Bstd.^2;
Brp1sum = sum(Brp1(rhobins).*w(rhobins))/sum(w(rhobins));

Bcontam = sum(bb(rhobins).*w(rhobins))/sum(w(rhobins));
rho = 0.1*Bcontam/Brp1sum; % Because model is r = 0.1

return