function calcplot_kbmsimxreal(simnums,calcopt)
% function calcplot_kbmsimxreal(simnums,calcopt)
% Calculate and plot beam map sim aps, maps pagers, and final spectra
% This only does individual-year Keck maps -- refer to bmsim_combine
% scripts for multi-year coadds
%
% INPUTS
%   
%   simnums: struct specifying the sim numbers to calc/plot 
%              simser, simrlz, simdau, realser, realrlz, realdau
%              For Keck we also require a year - make sure this matches
%              the sim!
%              k2012, k2013, k2014, k2015
%              Special cases: 2012+2013 = k13, k13 + b12 = bk13
%              bk13 + k2014 = bk14, bk14 + k2015 = bk15
%              These cases change the name of the map file read
%              "sim" refers to beam map sim - one realization
%              "real" refers to whatever is being crossed with "sim"
%              This could be real data (realrlz = 'real') or sims, e.g.
%              realrlz = '???8' for type 8 dust sims
%   calcopt: struct specifying spectrum calculation option and what to do
%              calcaps: 0 (default), 1 - generate aps
%              updateaps: 0 (default) remakes all aps, 1 only makes missing
%              farmaps: 0 (default) interactively, 1 to farm
%              makesimset: 0 (default), 1 - makesimset from realrlz
%              mappager: 0 (default), 1 - generate maps pager
%              plotaps: 0 (default), 1 - generate aps spectrum plots
%              est: 'pureB' (default), 'matrix'
%              autocross: 'auto' (default), 'cross'
%              crosstype: 'real' (default) to only cross with 1 real rlz
%                         'sim' to cross the real bm sim with all sim rlz
%              simstocross: realizations of sims to cross (only works
%                           with crosstype = 'sim')
%              perwhat: 'perfreq' (default), 'perrx'

% Parse input options
if ~isfield(simnums,'year')
  year = 'k2012';
else
  year = simnums.year;
end
% If we are in a cross-year map coadd, we probably are not coadding by rx
% so the map files will not have jack?1 - modify filename accordingly
if ismember(year,{'k2012','k2013','k2014','k2015'})
  coaddrx = '1';
else
  coaddrx = '';
end
if ~isfield(simnums,'simser')
  simser = '2553';
else
  simser = simnums.simser;
end
if ~isfield(simnums,'simrlz')
  simrlz = '0001';
else
  simrlz = simnums.simrlz;
end
if ~isfield(simnums,'simdau')
  simdau = 'b';
else
  simdau = simnums.simdau;
end
if ~isfield(simnums,'realser')
  realser = '1351';
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
if ~isfield(calcopt,'perwhat')
  perwhat = 'perfreq';
else
  perwhat = calcopt.perwhat;
end

% Matrices for each frequency/year - not all were made per-year or
% per-rx, so choose the closest for each case

% 95
bk14_95 = ['/n/panlfs2/bicep/keck/pipeline/matrixdata/c_t/1351/' ...
           'healpix_red_spectrum_lmax700_beamKuber100rev1_reob1351_' ...
           'd_100GHz_QQQUUU_proj.mat'];
bk15_95 = ['/n/panlfs2/bicep/keck/pipeline/matrixdata/c_t/1459/' ...
           'healpix_red_spectrum_lmax700_beamKuber100rev1_reob1459_' ...
           'aabde_100GHz_QQQUUU_proj.mat'];

% 150
% BICEP2 
b12_150 = ['/n/panlfs2/bicep/bicep2/pipeline/matrixdata/c_t/0704/' ...
           'healpix_red_spectrum_lmax700_beamB2bbns_reob0704_proj.mat']; 
% K13 - this was the original tag-subset 0706, but later on K2012/K2013
% observing matrices were made with flavor-reweight coadds and the
% subsequent BK14/BK15 150 matrices used these instead.
k13_150 = ['/n/panlfs2/bicep/keck/pipeline/matrixdata/c_t/0706/' ...
           'healpix_red_spectrum_lmax700_beamB2bbns_reob0706_rxa_proj' ...
           '.mat'];
% BK14_150
bk14_150 = ['/n/panlfs2/bicep/keck/pipeline/matrixdata/c_t/1459/' ...
            'healpix_red_spectrum_lmax700_beamB2bbns_reob1459_' ...
            'aabd_150GHz_QQQUUU_proj.mat'];
% BK15_150
bk15_150 = ['/n/panlfs2/bicep/keck/pipeline/matrixdata/c_t/1459/' ...
            'healpix_red_spectrum_lmax700_beamB2bbns_reob1459_' ...
            'aabde_150GHz_QQQUUU_proj.mat'];

% 220
bk15_220 = ['/n/panlfs2/bicep/keck/pipeline/matrixdata/c_t/1459/' ...
            'healpix_red_spectrum_lmax700_beamKuber220_reob1459_' ...
            'aabde_220GHz_QQQUUU_proj.mat'];

% apsopt template
load(['aps/0751x1351/' ...
      'real_a_filtp3_weight3_gs_dp1102_jack0_' ...
      'real_a_filtp3_weight3_gs_dp1102_jack01_matrix']);
apsopt = rmfield(apsopt,'purifmatname');
apsopt.makebpwf = 0;

% Assemble matrix if necessary, by default suitable for auto
switch est
  case 'pureB'
    apsopt.pure_b = 'kendrick';
  case 'matrix'
    switch year
      case {'k2012','k2013','k13','bk13'} % Just have coadded K2012+K2013 = K13
        switch perwhat
          case 'perfreq'
            apsopt.purifmatname = {k13_150};
          case 'perrx' % 
            apsopt.purifmatname = {k13_150,k13_150,k13_150,...
                                   k13_150,k13_150};
        end
      case {'k2014','bk14'} % Again, no per-year 150
        switch perwhat
          case 'perfreq'
            apsopt.purifmatname = {bk14_95,bk14_150};
          case 'perrx'
            apsopt.purifmatname = {bk14_95,bk14_150,bk14_95,...
                                   bk14_150,bk14_150};
        end
      case {'k2015','bk15'}
        switch perwhat
          case 'perfreq'
            apsopt.purifmatname = {bk15_95,bk15_150,bk15_220};
          case 'perrx'
            apsopt.purifmatname = {bk15_95,bk15_220,bk15_95,...
                                   bk15_220,bk15_150};
        end
    end
end

% Set up other apsopt - by default suitable for auto and maps pager
switch perwhat
  case 'perfreq'
    apsopt.coaddrx = 2;
    switch year
      case {'k2012','k2013'}
        apsopt.polrot = [0];
      case 'k2014'
        apsopt.polrot = [0,0];
      case 'k2015'
        apsopt.polrot = [0,0,0];
      case {'k13','bk13'}
        apsopt.polrot = [0];
        apsopt.coaddrx = 0; % switch back, it's already a single map
      case {'bk14'}
        apsopt.polrot = [0,0];
        apsopt.coaddrx = 0;
      case {'bk15'}
        apsopt.polrot = [0,0,0];
        apsopt.coaddrx = 0;
    end
  case 'perrx'
    apsopt.coaddrx = 0;
    apsopt.polrot = [0,0,0,0,0];
end % perwhat

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
      if coaddrx == '1'
        apsopt.ukpervolt = {get_ukpervolt(year(end-3:end))};
      else
        apsopt = rmfield(apsopt,'ukpervolt');
      end
      mapn = {[simser '/' simrlz '_' simdau '_filtp3_weight3_gs_dp*' ...
               '_jack?' coaddrx '.mat']};
      if farmaps
        mkdir(['farmfiles/' simser]);
        cmd = 'reduc_makeaps(mapn,apsopt)';
        farmfile = ['farmfiles/' simser '/' simrlz '_' simdau ...
                    '_auto_' est '_' perwhat];
        farmit(farmfile,cmd,'var',{'mapn','apsopt'},'mem',40000, ...
               'maxtime',1440);
      else
        reduc_makeaps(mapn,apsopt);
      end
    case 'cross' 
      apsopt.commonmask = 2;
      if coaddrx == '1'
        apsopt.ukpervolt = {get_ukpervolt(year(end-3:end)),get_ukpervolt(year(end-3:end))};
      else
        apsopt = rmfield(apsopt,'ukpervolt');
      end
      % Reassemble polrots and matrices for cross
      switch perwhat
        case 'perfreq' % There are n_freqs + n_freqs maps
          switch year
            case {'k2012','k2013','k13','bk13'}
              apsopt.polrot = [0,0];
              apsopt.purifmatname = {k13_150,k13_150};
            case {'k2014','bk14'}
              apsopt.polrot = [0,0,-0.5,0];
              apsopt.purifmatname = {bk14_95,bk14_150,bk14_95,bk14_150};
            case {'k2015','bk15'}
              apsopt.polrot = [0,0,0,-0.5,0,0];
              apsopt.purifmatname = {bk15_95,bk15_150,bk15_220,...
                                     bk15_95,bk15_150,bk15_220};
          end
        case 'perrx' % Here we always have 5 rxs + 5 rxs
          % Don't know if we have per-rx polrots...
          % If only we had per-rx matrices
          apsopt.polrot = [0,0,0,0,0,0,0,0,0,0]; % lol
          switch year
            case {'k2012','k2013'}
              apsopt.purifmatname = {k13_150,k13_150,k13_150,k13_150,k13_150,...
                                     k13_150,k13_150,k13_150,k13_150,k13_150};
            case 'k2014'
              apsopt.purifmatname = {bk14_95,bk14_150,bk14_95,bk14_150,bk14_150,...
                                     bk14_95,bk14_150,bk14_95,bk14_150,bk14_150};
            case 'k2015'
              apsopt.purifmatname = {bk15_95,bk15_220,bk15_95,bk15_220,bk15_150,...
                                     bk15_95,bk15_220,bk15_95,bk15_220,bk15_150};
          end
      end

      % Janky pivot back to pureB if required - but now we have the polrots!
      switch est
        case 'pureB';
          apsopt.pure_b = 'kendrick';
          apsopt = rmfield(apsopt,'purifmatname');
      end
          
      switch crosstype
        case 'real'
          % Not all dps were made for real data, but in general dp0000,
          % dp1100, and dp1102 should exist
          mapn1 = {[simser '/' simrlz '_' simdau '_filtp3_weight3_gs_' ...
                   'dp0000_jack0' coaddrx '.mat'],...
                  [simser '/' simrlz '_' simdau],...
                  [realser '/' realrlz '_' realdau]};
          mapn2 = {[simser '/' simrlz '_' simdau '_filtp3_weight3_gs_' ...
                   'dp1100_jack0' coaddrx '.mat'],...
                  [simser '/' simrlz '_' simdau],...
                  [realser '/' realrlz '_' realdau]};
          mapn3 = {[simser '/' simrlz '_' simdau '_filtp3_weight3_gs_' ...
                   'dp1102_jack0' coaddrx '.mat'],...
                  [simser '/' simrlz '_' simdau],...
                  [realser '/' realrlz '_' realdau]};
          if farmaps
            mkdir(['farmfiles/' simser]);
            cmd = 'reduc_makeaps(mapn1,apsopt)';
            farmfile = ['farmfiles/' simser '/' simrlz '_' simdau ...
                        '_crossreal_dp0000_' est '_' perwhat];
            farmit(farmfile,cmd,'var',{'mapn1','apsopt'},'mem',40000, ...
                   'maxtime',1440);
            cmd = 'reduc_makeaps(mapn2,apsopt)';
            farmfile = ['farmfiles/' simser '/' simrlz '_' simdau ...
                        '_crossreal_dp1100_' est '_' perwhat];
            farmit(farmfile,cmd,'var',{'mapn2','apsopt'},'mem',40000, ...
                   'maxtime',1440);
            cmd = 'reduc_makeaps(mapn3,apsopt)';
            farmfile = ['farmfiles/' simser '/' simrlz '_' simdau ...
                        '_crossreal_dp1102_' est '_' perwhat];
            farmit(farmfile,cmd,'var',{'mapn3','apsopt'},'mem',40000, ...
                   'maxtime',1440);
          else
            reduc_makeaps(mapn1,apsopt);
            reduc_makeaps(mapn2,apsopt);
            reduc_makeaps(mapn3,apsopt);
          end
          
        case 'sim'
          % For sims, just do dp1100 / jack0
          type = realrlz(4); % sim type
          for ii = 1:length(simstocross)
            mapreplacenum = sprintf('%03i',simstocross(ii)); % sim number
            mapn = {[simser '/' simrlz '_' simdau '_filtp3_weight3_gs_dp1100_jack0' ...
                     coaddrx '.mat'],...
                    [simser '/' simrlz '_' simdau],...
                    [realser '/' mapreplacenum type '_' realdau]};
            if farmaps
              mkdir(['farmfiles/' simser]);
              cmd = 'reduc_makeaps(mapn,apsopt)';
              farmfile = ['farmfiles/' simser '/' simrlz '_' simdau '_x_' ...
                          realser '_' mapreplacenum type '_' realdau '_' ...
                          est '_' perwhat];
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

switch perwhat
  case 'perfreq'
    switch year
      case {'k13','bk13','bk14','bk15'}
        coaddsuffix = '';
      otherwise
        coaddsuffix = '_overfreq';
    end
  case 'perrx'
    coaddsuffix = '';
end

if makesimset
  apsname = [simser 'x' realser '/' simrlz '_' simdau ...
             '_filtp3_weight3_gs_dp1100_jack0' coaddrx '_???' ...
             realrlz(4) '_' realdau '_filtp3_weight3_gs_dp1100_jack0' ...
             coaddrx '_' est '_cm' coaddsuffix '.mat'];
  reduc_makesimset(apsname);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Make maps pager - here always auto

mapn = {[simser '/' simrlz '_' simdau '_filtp3_weight3_gs_dp*' ...
         '_jack?' coaddrx '.mat']};
switch perwhat
  case 'perfreq'
    switch year
      case {'k2012','k2013','k13','bk13'}
        apsopt.mapname = {'150'};
      case {'k2014','bk14'}
        apsopt.mapname = {'95','150'};
      case {'k2015','bk15'}
        apsopt.mapname = {'95','150','220'};
    end
end

if coaddrx == '1'
  apsopt.ukpervolt = {get_ukpervolt(year(end-3:end))};
else
  if isfield(apsopt,'ukpervolt')
    apsopt = rmfield(apsopt,'ukpervolt');
  end
end

if mappager % free, then fixed color scale
  reduc_plotcomap_pager(mapn,apsopt,0,[],[],[],[],[],[],0,[],1); 
  reduc_plotcomap_pager(mapn,apsopt,1,[],[],[],[],[],[],0,[],1);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot aps

nbase = simser;
nbase2 = [simser 'x' realser];

% Deprojections - all made for auto
depj1{1} = '0000'; depj1{2} = '1000'; depj1{3} = '0100'; 
depj1{4} = '1100'; depj1{5} = '1110'; depj1{6} = '1101'; 
depj1{7} = '1111'; depj1{8} = '1102';

% Only a few for cross
depj2{1} = '0000'; depj2{2} = '1100'; depj2{3} = '1102';

% Same for jacks
jack{1} = '0'; jack{2} = '1'; jack{3} = '2'; jack{4} = '3';
jack{5} = '4'; jack{6} = '5'; jack{7} = '6'; jack{8} = '7';
jack{9} = '8'; jack{10} = '9'; jack{11} = 'a'; jack{12} = 'b';
jack{13} = 'c'; jack{14} = 'd'; jack{15} = 'e'; 

lmax{1} = '200'; lmax{2} = '500';

% Define suffixes for loading in aps
switch perwhat
  case 'perfreq'
    pwsuffix = '_overfreq';
    % But of course if we have K13 cross-year, the aps does NOT have this
    % cause it's coaddtype = 0...
    switch year
      case {'k13','bk13','bk14','bk15'}
        pwsuffix = '';
    end
  case 'perrx'
    pwsuffix = '';
end

% And the length/names of the aps to plot
% This can get way more complicated with crosses...shit
switch perwhat
  case 'perrx' % Def not doing all inter-rx crosses
    freq{1} = 'rx0'; freq{2} = 'rx1'; freq{3} = 'rx2';
    freq{4} = 'rx3'; freq{5} = 'rx4';
  case 'perfreq'
    switch year
      case {'k2012','k2013','k13','bk13'}
        freq{1} = '150';
      case {'k2014','bk14'}
        freq{1} = '100'; freq{2} = '150'; freq{3} = '100x150';
      case {'k2015','bk15'} % get these from aps_get_xind
        freq{1} = '100'; freq{2} = '150'; freq{3} = '220';
        freq{4} = '100x150'; freq{5} = '100x220'; freq{6} = '150x220';
    end
end

% BICEP2!!! Change when we have same for Keck
% Get the beam map sim algorithmic floor from noiseless constructed
% Planck sky dp1111 beam map sim
floor = load(['/n/bicepfs1/users/clwong/aps/2027/'...
              '0002_a_filtp3_weight3_gs_dp1111_jack0_pureB.mat']);
% Supfac correct it
flc = load(['final/1450x1350/real_a_filtp3_weight3_gs_dp1111_jack0_'...
            'real_a_filtp3_weight3_gs_dp1111_jack01_pureB_overrx.mat']);
floor.aps.Cs_l = floor.aps.Cs_l.*(flc.r(1).rwf);
% So of course this must be removed from real aps AFTER supfac correction

% NOPE, no floor is better than the wrong one...
floor.aps.Cs_l = zeros(size(floor.aps.Cs_l));

% Load up the sim-dependent algorithmic floor and beam map noise bias
% This should be supfac corrected
%[floor bias] = get_floorbias(expt);

if plotaps

  set(0,'defaultlinelinewidth',1.0);
  set(0,'DefaultAxesFontSize',14);
  
  figure(1); clf; 
  setwinsize(gcf,1000,600);
  set(gcf,'Visible','off');
  
  for dd = 1:length(depj1)

    % Load up bias for Keck once calculated
    
    % Load up the appropriate final file for supfac - check if multiple
    % dp options are available!
    % For K13 just use 2013
    switch est
      case 'pureB'
        switch year
          case 'k2012' 
            switch perwhat
              case 'perfreq' % B2, Keck, B2xKeck
                final = load(['final/0751x1351/real_a_filtp3_weight3_gs_dp' ...
                              '1102_jack0_real_a_filtp3_weight3_gs_dp1102_' ...
                              'jack01_pureB_overrx']);
                final.supfac = final.supfac([2]); % Just choose Keck
                final.r = final.r(2);
              case 'perrx' % This bad boy is B2 X all 5 rxs!
                final = load(['final/0751x1351/real_a_filtp3_weight3_gs_dp' ...
                              '1102_jack0_real_a_filtp3_weight3_gs_dp1102_' ...
                              'jack01_pureB']);
                final.supfac = final.supfac([2,3,4,5,6]); % All 5 rxs
                final.r = final.r([2,3,4,5,6]);
            end
          case {'k2013','k13'}
            switch perwhat
              case 'perfreq' % B2, Keck, B2xKeck
                final = load(['final/0751x1351/real_a_filtp3_weight3_gs_dp' ...
                              '1102_jack0_real_b_filtp3_weight3_gs_dp1102_' ...
                              'jack01_pureB_overrx']);
                final.supfac = final.supfac([2]); 
                final.r = final.r(2);
              case 'perrx' % This bad boy is B2 X all 5 rxs!
                final = load(['final/0751x1351/real_a_filtp3_weight3_gs_dp' ...
                              '1102_jack0_real_b_filtp3_weight3_gs_dp1102_' ...
                              'jack01_pureB']);
                final.supfac = final.supfac([2,3,4,5,6]);
                final.r = final.r([2,3,4,5,6]);
            end
          case {'k2014','bk14'}
            % No per-rx versions of this were made, so we have to cobble
            % it out of the overfreq versions
            % Apparently directbpwf version of this is available...
            final = load(['final/1351/real_d_filtp3_weight3_gs_dp' ...
                          '1102_jack01_pureB_overfreq.mat']);
            switch perwhat
              case 'perfreq' % 100, 150, 100x150
                final.supfac = final.supfac([1,2,3]); 
              case 'perrx'
                % rearrange to a 5-struct - 100, 150, 100, 150, 150
                final.supfac = [final.supfac(1),final.supfac(2),final.supfac(1),...
                                final.supfac(2),final.supfac(2)];
                % Do the same for the sim bandpowers
                final.r = [final.r(1),final.r(2),final.r(1),...
                           final.r(2),final.r(2)];
            end % perwhat
          case 'k2015'
            final = load(['final/1351/real_e_filtp3_weight3_gs_dp' ...
                          '1102_jack01_pureB_overfreq.mat']);
            switch perwhat
              case 'perfreq' % 100, 150, 220, 100x150, 100x220, 150x220
                final.supfac = final.supfac([1,2,3,4,5,6]);
              case 'perrx'
                final.supfac = [final.supfac(1),final.supfac(3),final.supfac(1),...
                                final.supfac(3),final.supfac(2)];
                final.r = [final.r(1),final.r(3),final.r(1),...
                           final.r(3),final.r(2)];
            end
        end % year

      case 'matrix'
        switch year
          case 'k2012' 
            switch perwhat
              case 'perfreq' % B2, Keck, B2xKeck
                final = load(['final/0751x1351/real_a_filtp3_weight3_gs_dp' ...
                              '1102_jack0_real_a_filtp3_weight3_gs_dp1102_' ...
                              'jack01_matrix_overrx']);
                final.supfac = final.supfac([2]); 
                final.r = final.r(2);
              case 'perrx' % This bad boy is B2 X all 5 rxs!
                final = load(['final/0751x1351/real_a_filtp3_weight3_gs_dp' ...
                              '1102_jack0_real_a_filtp3_weight3_gs_dp1102_' ...
                              'jack01_matrix']);
                final.supfac = final.supfac([2,3,4,5,6]);
                final.r = final.r([2,3,4,5,6]);
            end
          case {'k2013','k13'}
            switch perwhat
              case 'perfreq' % B2, Keck, B2xKeck
                final = load(['final/0751x1351/real_a_filtp3_weight3_gs_dp' ...
                              '1102_jack0_real_b_filtp3_weight3_gs_dp1102_' ...
                              'jack01_matrix_overrx']);
                final.supfac = final.supfac([2]); 
                final.r = final.r(2);
              case 'perrx' % This bad boy is B2 X all 5 rxs!
                final = load(['final/0751x1351/real_a_filtp3_weight3_gs_dp' ...
                              '1102_jack0_real_b_filtp3_weight3_gs_dp1102_' ...
                              'jack01_matrix']);
                final.supfac = final.supfac([2,3,4,5,6]);
                final.r = final.r([2,3,4,5,6]);
            end
          case {'k2014','bk14'}
            % STF told me to use this last year, and it's probably fine
            % for 95 -- but why not use the equivalent of the pureB
            % above?  Presumably the supfac will change for 3 Rxs worth
            % of 150 in 2014, and the below is coadded aabd...
            final = load(['final/1459x1615/real_aabd_filtp3_weight3_gs_' ...
                          'dp1102_jack0_real_aabd_filtp3_weight3_gs_' ...
                          'dp1100_jack0_pureB_matrix_cm_directbpwf.mat']);
            switch perwhat
              case 'perfreq' % 100, 150, 100x150
                final.supfac = final.supfac([1,2,8]);
              case 'perrx'
                % rearrange to a 5-struct - 100, 150, 100, 150, 150
                final.supfac = [final.supfac(1),final.supfac(2),final.supfac(1),...
                                final.supfac(2),final.supfac(2)];
                % Do the same for the sim bandpowers
                final.r = [final.r(1),final.r(2),final.r(1),...
                           final.r(2),final.r(2)];
            end % perwhat
          case 'k2015'
            final = load(['final/1459x1615/real_aabde_filtp3_weight3_gs_' ...
                          'dp1102_jack0_real_aabde_filtp3_weight3_gs_' ...
                          'dp1100_jack0_pureB_matrix_cm_directbpwf.mat']);
            switch perwhat
              case 'perfreq' % 100, 150, 220, 100x150, 100x220, 150x220
                final.supfac = final.supfac([1,2,3,13,14,24]);
              case 'perrx'
                % 100, 220, 100, 220, 150
                final.supfac = [final.supfac(1),final.supfac(3),final.supfac(1),...
                                final.supfac(3),final.supfac(2)];
                final.r = [final.r(1),final.r(3),final.r(1),...
                           final.r(3),final.r(2)];
            end % perwhat
        end % year
    end % est
    
    for jj = 1:length(jack)
        
      % Load in aps, apply supfac, prepare for plotting function
      switch autocross
        case 'auto'
          aps = load(['aps/' simser '/' simrlz '_' simdau ...
                      '_filtp3_weight3_gs_dp' depj1{dd} '_jack' ...
                      jack{jj} coaddrx '_' est pwsuffix '.mat']); 
          switch year
            case {'k2012','k2013','k13'}
              switch perwhat
                case 'perfreq'
                  aps = apply_supfac(aps.aps,final.supfac);
                case 'perrx'
                  aps.aps = aps.aps([1,2,3,4,5]);
                  aps = apply_supfac(aps.aps,final.supfac);
              end
            case {'k2014','bk14'}
              switch perwhat
                case 'perfreq'
                  aps.aps = aps.aps([1,2,3]); % to match supfac
                  aps.aps(3).Cs_l = aps.aps(3).Cs_l(:,1:6); % kill alt
                  final.supfac(3).rwf = final.supfac(3).rwf(:,1:6);
                  aps = apply_supfac(aps.aps,final.supfac);
                case 'perrx'
                  aps.aps = aps.aps([1,2,3,4,5]); % to match supfac
                  aps = apply_supfac(aps.aps,final.supfac);
              end
            case 'k2015'
              switch perwhat
                case 'perfreq'
                  aps.aps = aps.aps([1,2,3,4,5,6]);
                  aps.aps(4).Cs_l = aps.aps(4).Cs_l(:,1:6);
                  aps.aps(5).Cs_l = aps.aps(5).Cs_l(:,1:6);
                  aps.aps(6).Cs_l = aps.aps(6).Cs_l(:,1:6);
                  final.supfac(4).rwf = final.supfac(4).rwf(:,1:6);
                  final.supfac(5).rwf = final.supfac(5).rwf(:,1:6);
                  final.supfac(6).rwf = final.supfac(6).rwf(:,1:6);
                  aps = apply_supfac(aps.aps,final.supfac);
                case 'perrx'
                  aps.aps = aps.aps([1,2,3,4,5]); % to match supfac
                  aps = apply_supfac(aps.aps,final.supfac);
              end
          end % year
          leg = nbase;
        case 'cross'
          % Dear god this is going to be complicated
          leg = {nbase,nbase2};
      end % autocross
        
      % Apply beam map floor, bias correction, and calculate rho
      % Is this appropriate for anything but jack0?
      switch jack{jj}
        case '0'
          % No bias yet
          
          % What indices do we apply the floor to?  I suppose anything
          % that is bmsim x bmsim, but not bmsim x real
          switch year
            case {'k2012','k2013','k13'}
              switch perwhat
                case 'perfreq'
                  f_ind = 1;
                case 'perrx'
                  f_ind = [1,2,3,4,5];
              end
            case 'k2014'
              switch perwhat
                case 'perfreq'
                  f_ind = [1,2,3];
                case 'perrx'
                  f_ind = [1,2,3,4,5];
              end
            case 'k2015'
              switch perwhat
                case 'perfreq'
                  f_ind = [1,2,3,4,5,6];
                case 'perrx'
                  f_ind = [1,2,3,4,5];
              end
          end
          
          % For each relevant spectrum, remove the floor
          for ii = 1:length(f_ind)
            aps(f_ind(ii)).Cs_l(:,4) = aps(f_ind(ii)).Cs_l(:,4) - ...
                                       floor.aps.Cs_l(:,4);
          end
          
          % Calculate rho - make it the same length as aps
          rho = calc_rho(aps,final);
          
        otherwise
          rho = [];
                  
      end % jack0 special stuff
      
      % Plot away!
      
      % Loop over plot ranges
      for ff = 1:length(freq)
        for ll = 1:length(lmax)
      
          % Plot all 6 spectra
          subplot(2,3,1)
          plot_bmspec(aps(ff),leg,'TT',lmax{ll},freq{ff},jack{jj},depj1{dd},[])
          subplot(2,3,2)
          plot_bmspec(aps(ff),leg,'TE',lmax{ll},freq{ff},jack{jj},depj1{dd},[])
          subplot(2,3,3)
          plot_bmspec(aps(ff),leg,'EE',lmax{ll},freq{ff},jack{jj},depj1{dd},[])
          subplot(2,3,4)
          %if ~isempty(rho)
          switch jack{jj}
            case '0'
              plot_bmspec(aps(ff),leg,'BB',lmax{ll},freq{ff},jack{jj}, ...
                          depj1{dd},rho(ff))
            otherwise
              plot_bmspec(aps(ff),leg,'BB',lmax{ll},freq{ff},jack{jj}, ...
                          depj1{dd},[])
          end
            %else
            %plot_bmspec(aps(ff),leg,'BB',lmax{ll},freq{ff},jack{jj},depj1{dd})
            %end
          subplot(2,3,5)
          plot_bmspec(aps(ff),leg,'TB',lmax{ll},freq{ff},jack{jj},depj1{dd},[])
          subplot(2,3,6)
          plot_bmspec(aps(ff),leg,'EB',lmax{ll},freq{ff},jack{jj},depj1{dd},[])
          

          % Even if we just have auto, save in subfolder nbase x real cause
          % we'll populate  it eventually...
          mkdir(['plots/' nbase2]);
          savename = ['plots/' nbase2 '/' simrlz '_' simdau '_filtp3_weight3' ...
                      '_gs_dp' depj1{dd} '_jack' jack{jj} '_' freq{ff} '_' ...
                      est '_' lmax{ll} '.eps'];
          print('-depsc2',savename);
          clf()
          
        end % plot range  
      end % freq
    end % jack
  end % dp
end % plotaps


return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function rho = calc_rho(aps,final)
% Calculate rho statistic - 
% Depends on the final file for r=0.1 sim means and std
% aps is already supfac-corrected

rhobins = 2:6;

for ii = 1:length(aps)
  % This assumes that fields of bb (i.e. aps struct) and final correspond
  % to the same spectra
  % See documentation in final.r.doc
  Brp1 = mean(final.r(ii).sig(:,4,:),3);
  Bstd = std(final.r(ii).sim(:,4,:),[],3);
  w = Brp1./Bstd.^2;
  Brp1sum = sum(Brp1(rhobins).*w(rhobins))/sum(w(rhobins));

  bb = aps(ii).Cs_l(:,4);
  Bcontam = sum(bb(rhobins).*w(rhobins))/sum(w(rhobins));
  rho(ii) = 0.1*Bcontam/Brp1sum; % Because model is r = 0.1
end
  
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [floor bias] = get_floorbias(year)
% Each sim/year has a slightly different floor and bias


return