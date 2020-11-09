function ffbm_makebeamparamfiles(compopt)
% function ffbm_makebeamparamfiles(compopt)
%
% Function to read fits from windowed maps, apply cuts generated in
% ffbm_makecuts, and save all/some of
% a) .mat file with all good measurements (good for plotting)
% b) .csv: saves sigma, p, c, corresponding errors, and 
%          a/b offset errors.  Should be committed to aux_data/beams
% c) fpdata_beam: .csv files made per-receiver which are read by
%                 get_array_info and used to cut on bad a/b offsets and beam
%                 shapes  
%
% As always, relies on bmrunlist_201?.csv
%
% INPUTS (should all be sent in with compopt)
%        Works best if you use the compopt saved with the cut struct
%
%   expt:            'keck','bicep3'
%   year:            keck: 2012,...,2018
%                    bicep3: 2016,...,2018
%
% OPTIONAL INPUTS
%   
%   mapcoord:        coordinate system of windowed maps to load
%                    'xpyp' (default), 'azel_ap'
%   demodtype:       'square' (default), 'choplet'
%   mapsize:         windowed map size in deg (8 default)
%   suffix:          optional string to add to map name, '' default
%   bmdir:           directory from which to load beammaps
%                    (default beammaps/maps)
%   mapstoload:      vector of beam map run numbers to look at, i.e. 1:10
%                    (default is everything in bmrunlist.csv)
%   cutdir:          directory from which to load the cut structure 
%                    (default beammaps/cuts)
%   cutfile:         file in cutdir to load
%                    (default ffbm_cuts_year)
%   cutlist:         list of things to cut on, i.e.
%                    {'mirror','nanfit','peak'} (default all)
%   beamparamdir:    directory in which to save output files
%                    (default beammaps/beamparams)
%   beamparamfile:   filename to save to
%                    (default beamparams_year)
%   makemat:         make .mat file with all good measurements
%                    (default 1)
%   makecsv:         make .csv file with summary numbers
%                    (default 0)
%   makefpdata:      make 'aboffsetflags' and 'beamshapeflags' files
%                    (default 0)

% Parse compopt
expt = compopt.expt;
year = compopt.year;
if ~isfield(compopt,'mapcoord');
  mapcoord = 'xpyp';
  compopt.mapcoord = mapcoord;
else
  mapcoord = compopt.mapcoord;
end
if ~isfield(compopt,'demodtype');
  demodtype = 'square';
else
  demodtype = compopt.demodtype;
end
if ~isfield(compopt,'mapsize');
  sizestr = '8deg';
else
  sizestr = [strrep(num2str(compopt.mapsize),'.','p') 'deg'];
end
if ~isfield(compopt,'suffix')
  suffix = '';
else
  suffix = compopt.suffix;
end
if ~isfield(compopt,'bmdir')
  bmdir = 'beammaps/maps_windowed/';
  compopt.bmdir = bmdir;
else
  bmdir = compopt.bmdir;
end
% Which maps do we want?
% Get beam map run info
bm = get_bm_info(year);
if ~isfield(compopt,'mapstoload')
  compopt.mapstoload = bm.number;
else
  bm = structcut(bm,compopt.mapstoload);
end
if ~isfield(compopt,'cutdir')
  compopt.cutdir = 'beammaps/cuts/';
  cutdir = compopt.cutdir;
else
  cutdir = compopt.cutdir;
end
if ~isfield(compopt,'cutfile')
  compopt.cutfile = ['ffbm_cuts_' num2str(year)];
  cutfile = compopt.cutfile;
else
  cutfile = compopt.cutfile;
end
if ~isfield(compopt,'cutlist')
  compopt.cutlist = {'mirror','notlight','nanfit','peak',...
      'sigma','ellip','abspoint','hand'};
  cutlist = compopt.cutlist;
else
  cutlist = compopt.cutlist;
end
if ~isfield(compopt,'beamparamdir')
  compopt.beamparamdir = 'beammaps/beamparams';
  beamparamdir = compopt.beamparamdir;
else
  beamparamdir = compopt.beamparamdir;
end
if ~isfield(compopt,'beamparamfile')
  compopt.beamparamfile = ['beamparams_' num2str(year)];
  beamparamfile = compopt.beamparamfile;
else
  beamparamfile = compopt.beamparamfile;
end
if ~isfield(compopt,'makemat')
  makemat = 1;
else
  makemat = compopt.makemat;
end
if ~isfield(compopt,'makecsv')
  makecsv = 0;
else
  makecsv = compopt.makecsv;
end
if ~isfield(compopt,'makefpdata')
  makefpdata = 0;
else
  makefpdata = compopt.makefpdata;
end

[p ind] = get_array_info([num2str(year) '0201']);
p = rmfield(p,'expt'); % for structcut to work

% Get the right dk for each detector (rotation option in findingparam)
p.theta = p.theta + p.drumangle;

% Load up cuts
cutname = [cutdir '/' cutfile];
load(cutname);
cut = cuts.cuts;

% Check that cut struct and mapstoload are consistent
if length(bm.number) ~= length(cut)
  error('Cut struct does not match length of requested maps');
end

% If we're in xp/yp, we don't need to rotate when calculating
% beam params
switch mapcoord
  case 'xpyp'
    rot = 0;
  case 'azel_ap'
    rot = 1;
  otherwise
    error('mapcoord must be xpyp or azel_ap')
end

% Make .mat file
if makemat

  % Load up each beam map.  Only use windowed map!
  for ii = 1:length(bm.number)
    
    filename = bm.filename{ii};
    disp(['Loading map ' num2str(ii) ': ' filename]);
    wmapname = [bmdir,'/mapwin_',filename,'_',demodtype,...
                '_',mapcoord,'_',sizestr,suffix];

    w = load(wmapname);
    w = w.bm;
    
    % Calculate beam parameters with rotation
    amp(:,ii) = w.A(1,:);
    x_bin_cent(:,ii) = w.A(2,:);
    y_bin_cent(:,ii) = w.A(3,:);
    [ss pp cc] = ffbm_findingparam(w.A,p,w.dk,rot,expt);
    sig(:,ii) = ss;
    e_p(:,ii) = pp;
    e_c(:,ii) = cc;
    % This step is a little redundant for mapcoord = xpyp
    [xx yy] = ffbm_calcdiffpoint(w.A,p,ind,w.dk,expt,rot);
    x(:,ii) = xx;
    y(:,ii) = yy;
    
    % Add up all the cuts
    totcut = zeros(length(p.gcp),1);
    for jj = 1:length(cutlist)
      totcut = totcut + getfield(cut(ii),cutlist{jj});
    end
    badind = find(totcut ~= 0);
    
    % Apply the cuts
    amp(badind,ii) = NaN;
    x_bin_cent(badind,ii) = NaN;
    y_bin_cent(badind,ii) = NaN;
    sig(badind,ii) = NaN;
    e_p(badind,ii) = NaN;
    e_c(badind,ii) = NaN;
    x(badind,ii) = NaN;
    y(badind,ii) = NaN;
  end
  
  % Prep output struct
  % Correct beamwidth for finite chopper aperture
  sig = correct_sigma(sig,expt,year);
  beamparam.amp = amp;
  beamparam.x_bin_cent = x_bin_cent;
  beamparam.y_bin_cent = y_bin_cent;
  beamparam.sig = sig;
  beamparam.e_p = e_p;
  beamparam.e_c = e_c;
  beamparam.x = x;
  beamparam.y = y;
  beamparam.compopt = compopt;
  beamparam.bm = bm;
  savename = [beamparamdir '/' beamparamfile];
  save(savename,'beamparam');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% That made the basic beamparam file, containing all individual
% measurements.  For everything else we just want summary statistics
% provided by ffbm_evalbeamparam, the output of which is stored in
% the beam .csv files, and which are plotted by ffbm_plot_beamparam and
% ffbm_plot_beamparam_tile

% Make csv/fp_data files -- only works if we saved the .mat file above
if(makecsv | makefpdata)
  
  % This loads the beamparam file and calculates statistics thereon
  params = ffbm_evalbeamparam(compopt);
  
  % Move back to r,theta coordinates
  [theta,r] = cart2pol(params.x_med,params.y_med);
  params.r = r;
  params.theta = theta * 180.0/pi;
  
  % Find non-NaN values
  cc = find(~isnan(params.aboffset_err_x));
  pr.gcp = p.gcp;
  pr.r = p.r;
  pr.theta = p.theta;
  pr.r(cc) = params.r(cc);
  pr.theta(cc) = params.theta(cc);
  pr.aboffset_err_x = NaN(length(params.sig_med),1);
  pr.aboffset_err_y = NaN(length(params.sig_med),1);
  pr.aboffset_err_x(cc) = params.aboffset_err_x(cc);
  pr.aboffset_err_y(cc) = params.aboffset_err_y(cc);
  pr.aboffset_err_x(find(pr.aboffset_err_x == 0)) = NaN;
  pr.aboffset_err_y(find(pr.aboffset_err_x == 0)) = NaN;
  
  % Assign median sigma for each band, then replace with measured if it
  % exists
  pr.sigma = NaN(length(params.sig_med),1);
  bands = unique(p.band(find(p.band))); 
  for ii = 1:length(bands)
    bandind = find(p.band == bands(ii));
    pr.sigma(bandind) = nanmedian(params.sig_med(bandind));
  end
  cc = find(~isnan(params.sig_med));
  pr.sigma(cc) = params.sig_med(cc);
  pr.p = zeros(length(params.sig_med),1);
  pr.p(cc) = params.e_p_med(cc);
  pr.c = zeros(length(params.sig_med),1);
  pr.c(cc) = params.e_c_med(cc);
  
  pr.err_sigma = zeros(length(params.sig_med),1);
  pr.err_sigma(cc) = params.sig_wid(cc);
  pr.err_p = zeros(length(params.sig_med),1);
  pr.err_p(cc) = params.e_p_wid(cc);
  pr.err_c = zeros(length(params.sig_med),1);
  pr.err_c(cc) = params.e_c_wid(cc);
  
  % Differential parameters: same for A/B
  % Since we subtract the template scaled by the diff coefficient, make sure
  % we don't have NaNs (which could come from the measurement) on any light
  % channels!  Better to not subtract at all than to kill the channel
  pr.dp = zeros(length(params.sig_med),1);
  tmp = params.de_p_med(ind.a);
  tmp(find(isnan(tmp))) = 0;
  pr.dp(ind.a) = tmp;
  pr.dp(ind.b) = tmp;
  pr.dc = zeros(length(params.sig_med),1);
  tmp = params.de_c_med(ind.a);
  tmp(find(isnan(tmp))) = 0;
  pr.dc(ind.a) = tmp;
  pr.dc(ind.b) = tmp;
  
  pr.err_dp = zeros(length(params.sig_med),1);
  tmp = params.de_p_wid(ind.a);
  tmp(find(isnan(tmp))) = 0;
  pr.err_dp(ind.a) = tmp;
  pr.err_dp(ind.b) = tmp;
  pr.err_dc = zeros(length(params.sig_med),1);
  tmp = params.de_c_wid(ind.a);
  tmp(find(isnan(tmp))) = 0;
  pr.err_dc(ind.a) = tmp;
  pr.err_dc(ind.b) = tmp;

  % NaNs for darks
  pr.sigma(ind.o) = NaN;
  pr.c(ind.o) = NaN;
  pr.p(ind.o) = NaN;
  pr.aboffset_err_x(ind.o) = NaN;
  pr.aboffset_err_y(ind.o) = NaN;
  pr.err_sigma(ind.o) = NaN;
  pr.err_p(ind.o) = NaN;
  pr.err_c(ind.o) = NaN;
  pr.dp(ind.o) = NaN;
  pr.dc(ind.o) = NaN;
  pr.err_dp(ind.o) = NaN;
  pr.err_dc(ind.o) = NaN;

  if makecsv
    switch expt
      case 'keck'
        % Get comments and save per-rx
        rx = unique(p.rx);
        for ii = 1:length(rx)
          k = get_csv_comments(year,rx(ii),expt);
          cutind = (p.rx == rx(ii));
          prx = structcut(pr,cutind);
          day = datestr(now,'yyyymmdd');
          filename = [beamparamdir '/beams_' expt '_obs' num2str(year) '_rx' ...
                      num2str(rx(ii)) '_' day '.csv'];
          ParameterWrite(filename,prx,k);
        end    
    end
    % Save total file
    k = get_csv_comments(year,'all',expt);
    filename = [beamparamdir '/beams_' expt '_obs_' num2str(year) '0101.csv'];
    ParameterWrite(filename,pr,k);
     
  end % makecsv
  
  if makefpdata % taken from ffbm_makefpdata_beam
    % Get FPU names
    switch expt
      case 'keck'
        switch year
          case 2012
            [P,K] = ParameterRead('aux_data/fp_data/fp_data_master_20111226.csv');
          case 2013
            [P,K] = ParameterRead('aux_data/fp_data/fp_data_master_20130201.csv');
          case 2014
            [P,K] = ParameterRead('aux_data/fp_data/fp_data_master_20131201.csv');
          case 2015
            [P,K] = ParameterRead('aux_data/fp_data/fp_data_master_20150101.csv');
          case 2016
            [P,K] = ParameterRead('aux_data/fp_data/fp_data_master_20160101.csv');
          case 2017
            [P,K] = ParameterRead('aux_data/fp_data/fp_data_master_20170101.csv');
        end
      case 'bicep3'
        switch year
          case 2016
            [P,K] = ParameterRead('aux_data/fp_data/fp_data_master_20151201.csv');
          case 2017
            [P,K] = ParameterRead('aux_data/fp_data/fp_data_master_20161201.csv');
        end
    end
    % Already have p/ind from earlier

    for ii = 1:length(P.rx)

      cutind = (p.rx == P.rx(ii));
      pcut = structcut(p,cutind);
      indcut = make_ind(pcut);
      % Cut down parameters (pr) struct
      paramscut = structcut(params,cutind);
      day = datestr(now,'yyyymmdd');

      % Load the real fp_data file for the first few fields
      % Using the stardate for rx4 after 2012 doesn't work...
      switch P.rx(ii)
	case 4
	  switch year
	    case {2013,2014,2015,2016}
	      P.startdate(ii) = 20130201;
	  end
      end
      fpdatain = ['aux_data/fp_data/fp_data_' P.fpu{ii} '_' num2str(P.startdate(ii)) '.csv'];
      [tmp,kk] = ParameterRead(fpdatain);
      
      pn.gcp = pcut.gcp;
      pn.tile = pcut.tile;
      pn.det_col = pcut.det_col;
      pn.det_row = pcut.det_row;
      pn.pol = pcut.pol;
      pn.mce_col = pcut.mce_col;
      pn.mce_row = pcut.mce_row;
      pn.type = pcut.type;
      
      % A/B Offset flags
      filename = ['fp_data_aboffsetflags_' P.fpu{ii} '_' num2str(year) '0101.csv'];
      kn = get_fpdata_comments(year,expt,'aboffsetflags',filename,kk);      
      pn.aboffset(indcut.e) = NaN;
      pn.aboffset(indcut.a) = sqrt(paramscut.aboffset_x(indcut.a).^2 + ...
	  paramscut.aboffset_y(indcut.a).^2);
      pn.aboffset(indcut.b) = pn.aboffset(indcut.a);
      % Median sigma for this band
      sigmarx = nanmedian(paramscut.sig_med);
      kn.comments{length(kn.comments)+1} = ['Median rx sigma: ' num2str(sigmarx) ' deg'];
      pn.aboffset_pct = pn.aboffset./sigmarx;
      filename = [beamparamdir '/' filename];
      ParameterWrite(filename,pn,kn);

      % Beam shape flags
      filename = ['fp_data_beamshapeflags_' P.fpu{ii} '_' num2str(year) '0101.csv'];
      kn = get_fpdata_comments(year,expt,'beamshapeflags',filename,kk);
      kn.comments{length(kn.comments)+1} = ['Median rx sigma: ' num2str(sigmarx) ' deg'];
      % For sigma, make sure we have nominal values for non-measured
      % (i.e. from above)
      prcut = structcut(pr,cutind);
      pn.sigma = prcut.sigma;
      pn.sigma_pct = pn.sigma./sigmarx;
      pn.sigma_diff(indcut.e,1) = NaN;
      pn.sigma_diff(indcut.a,1) = pn.sigma(indcut.a) - pn.sigma(indcut.b);
      pn.sigma_diff(indcut.b,1) = pn.sigma_diff(indcut.a);
      pn.sigma_diff_pct = pn.sigma_diff./sigmarx;
      % For ellipticity, make sure non-measured are zeros (again from above)
      pn.ellip = sqrt(prcut.p.^2 + prcut.c.^2);
      pn.ellip_diff(indcut.e,1) = NaN;
      pn.ellip_diff(indcut.a,1) = pn.ellip(indcut.a) - pn.ellip(indcut.b);
      pn.ellip_diff(indcut.b,1) = pn.ellip_diff(indcut.a);
      pn.pc_diff(indcut.e,1) = NaN;
      pn.pc_diff(indcut.a,1) = sqrt(prcut.dp(indcut.a).^2 + prcut.dc(indcut.a).^2);
      pn.pc_diff(indcut.b,1) = pn.pc_diff(indcut.a);
      filename = [beamparamdir '/' filename];
      ParameterWrite(filename,pn,kn);
    end
  end % makefpdata
  
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function k = get_csv_comments(year,rxNum,expt)

switch expt
  case 'keck'
    
    switch rxNum
      case 'all'
        
        switch year
          case 2014
            k.comments{1} = '# KECK 2014 all receivers';
            k.comments{2} = '# Data from thermal beam map data taken Feb/Mar 2014';
            k.comments{3} = '# This file is a concatenation of individual receiver beams files';
            k.comments{4} = '# Since the beam maps do not contain absolute pointing information,';
            k.comments{5} = '# the pointing information (r,theta) is derived from adding AB pointing offsets to ideal pointing information';
            k.comments{6} = '# Analysis found here:';
            k.comments{7} = '# ~/spuder/keck_analysis_logbook/analysis/20150602_keck2014_beamparam';
            k.comments{8} = '# generated by script:';
            k.comments{9} = '# keck_analysis/beammap/ffbm_makebeamparamfiles.m';
            k.comments{10} = ['# generated ' datestr(clock)];
            k.filename = ['beams_keck_obs_' num2str(year) '0101.csv'];
            k.created = [datestr(clock,'yyyymmdd') ' KSK'];
            k.fields = {'gcp','r','theta','sigma','p','c','dp','dc',...
                        'aboffset_err_x','aboffset_err_y','err_sigma','err_p','err_c','err_dp','err_dc'};
            k.units = {'(null)','deg','deg','deg','unitless','unitless','unitless','unitless',...
                       'deg','deg','deg','unitless','unitless','unitless','unitless'};
            k.formats = {'integer','double','double','double','double','double','double','double',...
                         'double','double','double','double','double','double','double'};
          case 2015
            k.comments{1} = '# KECK 2015 all receivers';
            k.comments{2} = '# Data from thermal beam map data taken Feb/Mar 2015';
            k.comments{3} = '# This file is a concatenation of individual receiver beams files';
            k.comments{4} = '# Since the beam maps do not contain absolute pointing information,';
            k.comments{5} = '# the pointing information (r,theta) is derived from adding AB pointing offsets to ideal pointing information';
            k.comments{6} = '# Analysis found here:';
            k.comments{7} = '# ~/spuder/keck_analysis_logbook/analysis/20150602_keck_2014_beam_maps';
            k.comments{8} = '# generated by script:';
            k.comments{9} = '# keck_analysis/beammap/ffbm_makebeamparamfiles.m';
            k.comments{10} = ['# generated ' datestr(clock)];
            k.filename = ['beams_keck_obs_' num2str(year) '0101.csv'];
            k.created = [datestr(clock,'yyyymmdd') ' KSK'];
            k.fields = {'gcp','r','theta','sigma','p','c','dp','dc',...
                        'aboffset_err_x','aboffset_err_y','err_sigma','err_p','err_c','err_dp','err_dc'};
            k.units = {'(null)','deg','deg','deg','unitless','unitless','unitless','unitless',...
                       'deg','deg','deg','unitless','unitless','unitless','unitless'};
            k.formats = {'integer','double','double','double','double','double','double','double',...
                         'double','double','double','double','double','double','double'};
          case 2016
            k.comments{1} = '# KECK 2016 all receivers';
            k.comments{2} = '# Data from thermal beam map data taken Feb/Mar 2016';
            k.comments{3} = '# This file is a concatenation of individual receiver beams files';
            k.comments{4} = '# Since the beam maps do not contain absolute pointing information,';
            k.comments{5} = '# the pointing information (r,theta) is derived from adding AB pointing offsets to ideal pointing information';
            k.comments{6} = '# Analysis found here:';
            k.comments{7} = '# ~/spuder/keck_analysis_logbook/analysis/20150602_keck_2014_beam_maps';
            k.comments{8} = '# generated by script:';
            k.comments{9} = '# keck_analysis/beammap/ffbm_makebeamparamfiles.m';
            k.comments{10} = ['# generated ' datestr(clock)];
            k.filename = ['beams_keck_obs_' num2str(year) '0101.csv'];
            k.created = [datestr(clock,'yyyymmdd') ' KSK'];
            k.fields = {'gcp','r','theta','sigma','p','c','dp','dc',...
                        'aboffset_err_x','aboffset_err_y','err_sigma','err_p','err_c','err_dp','err_dc'};
            k.units = {'(null)','deg','deg','deg','unitless','unitless','unitless','unitless',...
                       'deg','deg','deg','unitless','unitless','unitless','unitless'};
            k.formats = {'integer','double','double','double','double','double','double','double',...
                         'double','double','double','double','double','double','double'};
          case 2017
            k.comments{1} = '# KECK 2017 all receivers';
            k.comments{2} = '# Data from thermal beam map data taken Feb/Mar 2017';
            k.comments{3} = '# This file is a concatenation of individual receiver beams files';
            k.comments{4} = '# Since the beam maps do not contain absolute pointing information,';
            k.comments{5} = '# the pointing information (r,theta) is derived from adding AB pointing offsets to ideal pointing information';
            k.comments{6} = '# Analysis found here:';
            k.comments{7} = '# bkcmb/analysis_logbook/analysis/20180917_2017_beamparamcomparison';
            k.comments{8} = '# generated by script:';
            k.comments{9} = '# keck_analysis/beammap/ffbm_makebeamparamfiles.m';
            k.comments{10} = ['# generated ' datestr(clock)];
            k.filename = ['beams_keck_obs_' num2str(year) '0101.csv'];
            k.created = [datestr(clock,'yyyymmdd') ' TSG'];
            k.fields = {'gcp','r','theta','sigma','p','c','dp','dc',...
                        'aboffset_err_x','aboffset_err_y','err_sigma','err_p','err_c','err_dp','err_dc'};
            k.units = {'(null)','deg','deg','deg','unitless','unitless','unitless','unitless',...
                       'deg','deg','deg','unitless','unitless','unitless','unitless'};
            k.formats = {'integer','double','double','double','double','double','double','double',...
                         'double','double','double','double','double','double','double'};
        end % year switch
        
      otherwise % per-rx
        
        switch year
          case 2012
            k.comments{1} = ['# KECK RX' num2str(rxNum) 'FPU'];
            k.comments{2} = '# Data from thermal beam map data taken Feb 2012';
            k.comments{3} = '# Since the beam maps do not contain absolute pointing information,';
            k.comments{4} = '# the pointing information (r,theta) is derived from adding AB pointing offsets to ideal pointing information';
            k.comments{5} = '# Analysis found here:';
            k.comments{6} = '# ~/spuder/keck_analysis_logbook/analysis/20130716_beamparam2012';
            k.comments{7} = '# generated by script:';
            k.comments{8} = '# ~clwong/20130703_beamparam/param_manyfiles.m';
            k.comments{9} = ['# generated ' datestr(clock)];
            k.filename = ['beams_keck_obs' num2str(year) '_rx' num2str(rxNum) '_20130703.csv'];
            k.created = [datestr(clock,'yyyymmdd') ' CLW'];
            k.fields = {'gcp' 'r' 'theta' 'sigma' 'p' 'c' 'aboffset_err_x' 'aboffset_err_y' 'err_sigma','err_p','err_c'};
            k.units = {'(null)' 'deg' 'deg' 'deg' 'unitless' 'unitless' 'deg' 'deg' 'deg' 'unitless' 'unitless'};
            k.formats = {'integer' 'double' 'double' 'double' 'double' 'double' ...
                         'double' 'double' 'double' 'double' 'double'};	
          case 2013
            k.comments{1} = ['# KECK RX' num2str(rxNum) 'FPU'];
            k.comments{2} = '# Data from thermal beam map data taken Feb 2013';
            k.comments{3} = '# Since the beam maps do not contain absolute pointing information,';
            k.comments{4} = '# the pointing information (r,theta) is derived from adding AB pointing offsets to ideal pointing information';
            k.comments{5} = '# Analysis found here:';
            k.comments{6} = '# ~/spuder/keck_analysis_logbook/analysis/20131114_beamparam2013';
            k.comments{7} = '# generated by script:';
            k.comments{8} = '# ~clwong/20131112_channelcuts2013/param_manyfiles.m';
            k.comments{9} = ['# generated ' datestr(clock)];
            k.filename = ['beams_keck_obs' num2str(year) '_rx' num2str(rxNum) '_20131113.csv'];
            k.created = [datestr(clock,'yyyymmdd') ' CLW'];
            k.fields = {'gcp' 'r' 'theta' 'sigma' 'p' 'c' 'aboffset_err_x' 'aboffset_err_y' 'err_sigma','err_p','err_c'};
            k.units = {'(null)' 'deg' 'deg' 'deg' 'unitless' 'unitless' 'deg' 'deg' 'deg' 'unitless' 'unitless'};
            k.formats = {'integer' 'double' 'double' 'double' 'double' 'double' ...
                         'double' 'double' 'double' 'double' 'double'};	
          case 2014
            day = datestr(now,'yyyymmdd');
            k.comments{1} = ['# KECK RX' num2str(rxNum) 'FPU'];
            k.comments{2} = '# Data from thermal beam map data taken Feb/Mar 2014';
            k.comments{3} = '# Since the beam maps do not contain absolute pointing information,';
            k.comments{4} = '# the pointing information (r,theta) is derived from adding AB pointing offsets to ideal pointing information';
            k.comments{5} = '# Analysis found here:';
            k.comments{6} = '# ~/spuder/keck_analysis_logbook/analysis/20150602_keck2014_beamparam';
            k.comments{7} = '# generated by script:';
            k.comments{8} = '# keck_analysis/beammap/ffbm_makebeamparamfiles.m';
            k.comments{9} = ['# generated ' datestr(clock)];
            k.filename = ['beams_keck_obs' num2str(year) '_rx' num2str(rxNum) ...
                          '_' day '.csv'];
            k.created = [datestr(clock,'yyyymmdd') ' KSK'];
            k.fields = {'gcp','r','theta','sigma','p','c','dp','dc',...
                        'aboffset_err_x','aboffset_err_y','err_sigma','err_p','err_c','err_dp','err_dc'};
            k.units = {'(null)','deg','deg','deg','unitless','unitless','unitless','unitless',...
                       'deg','deg','deg','unitless','unitless','unitless','unitless'};
            k.formats = {'integer','double','double','double','double','double','double','double',...
                         'double','double','double','double','double','double','double'};
          case 2015
            day = datestr(now,'yyyymmdd');
            k.comments{1} = ['# KECK RX' num2str(rxNum) 'FPU'];
            k.comments{2} = '# Data from thermal beam map data taken Feb/Mar 2015';
            k.comments{3} = '# Since the beam maps do not contain absolute pointing information,';
            k.comments{4} = '# the pointing information (r,theta) is derived from adding AB pointing offsets to ideal pointing information';
            k.comments{5} = '# Analysis found here:';
            k.comments{6} = '# ~/spuder/keck_analysis_logbook/analysis/20150602_keck_2014_beam_maps';
            k.comments{7} = '# generated by script:';
            k.comments{8} = '# keck_analysis/beammap/ffbm_makebeamparamfiles.m';
            k.comments{9} = ['# generated ' datestr(clock)];
            k.filename = ['beams_keck_obs' num2str(year) '_rx' num2str(rxNum) ...
                          '_' day '.csv'];
            k.created = [datestr(clock,'yyyymmdd') ' KSK'];
            k.fields = {'gcp','r','theta','sigma','p','c','dp','dc',...
                        'aboffset_err_x','aboffset_err_y','err_sigma','err_p','err_c','err_dp','err_dc'};
            k.units = {'(null)','deg','deg','deg','unitless','unitless','unitless','unitless',...
                       'deg','deg','deg','unitless','unitless','unitless','unitless'};
            k.formats = {'integer','double','double','double','double','double','double','double',...
                         'double','double','double','double','double','double','double'};
          case 2016
            day = datestr(now,'yyyymmdd');
            k.comments{1} = ['# KECK RX' num2str(rxNum) 'FPU'];
            k.comments{2} = '# Data from thermal beam map data taken Feb/Mar 2016';
            k.comments{3} = '# Since the beam maps do not contain absolute pointing information,';
            k.comments{4} = '# the pointing information (r,theta) is derived from adding AB pointing offsets to ideal pointing information';
            k.comments{5} = '# Analysis found here:';
            k.comments{6} = '# ~/spuder/keck_analysis_logbook/analysis/20150602_keck_2014_beam_maps';
            k.comments{7} = '# generated by script:';
            k.comments{8} = '# keck_analysis/beammap/ffbm_makebeamparamfiles.m';
            k.comments{9} = ['# generated ' datestr(clock)];
            k.filename = ['beams_keck_obs' num2str(year) '_rx' num2str(rxNum) ...
                          '_' day '.csv'];
            k.created = [datestr(clock,'yyyymmdd') ' KSK'];
            k.fields = {'gcp','r','theta','sigma','p','c','dp','dc',...
                        'aboffset_err_x','aboffset_err_y','err_sigma','err_p','err_c','err_dp','err_dc'};
            k.units = {'(null)','deg','deg','deg','unitless','unitless','unitless','unitless',...
                       'deg','deg','deg','unitless','unitless','unitless','unitless'};
            k.formats = {'integer','double','double','double','double','double','double','double',...
                         'double','double','double','double','double','double','double'};
        case 2017
            day = datestr(now,'yyyymmdd');
            k.comments{1} = ['# KECK RX' num2str(rxNum) 'FPU'];
            k.comments{2} = '# Data from thermal beam map data taken Feb/Mar 2017';
            k.comments{3} = '# Since the beam maps do not contain absolute pointing information,';
            k.comments{4} = '# the pointing information (r,theta) is derived from adding AB pointing offsets to ideal pointing information';
            k.comments{5} = '# Analysis found here:';
            k.comments{6} = '# bkcmb/analysis_logbook/analysis/20180917_2017_beamparamcomparison';
            k.comments{7} = '# generated by script:';
            k.comments{8} = '# keck_analysis/beammap/ffbm_makebeamparamfiles.m';
            k.comments{9} = ['# generated ' datestr(clock)];
            k.filename = ['beams_keck_obs' num2str(year) '_rx' num2str(rxNum) ...
                          '_' day '.csv'];
            k.created = [datestr(clock,'yyyymmdd') ' TSG'];
            k.fields = {'gcp','r','theta','sigma','p','c','dp','dc',...
                        'aboffset_err_x','aboffset_err_y','err_sigma','err_p','err_c','err_dp','err_dc'};
            k.units = {'(null)','deg','deg','deg','unitless','unitless','unitless','unitless',...
                       'deg','deg','deg','unitless','unitless','unitless','unitless'};
            k.formats = {'integer','double','double','double','double','double','double','double',...
                         'double','double','double','double','double','double','double'};
        end % year switch
    end % rxNum switch

  case 'bicep3'
    switch year
      case 2016
        k.comments{1} = '# BICEP3 2016';
        k.comments{2} = '# Data from thermal beam map data taken Feb/Mar 2016';
        k.comments{3} = '# Since the beam maps do not contain absolute pointing information,';
        k.comments{4} = '# the pointing information (r,theta) is derived from adding AB pointing offsets to ideal pointing information';
        k.comments{5} = '# Analysis found here:';
        k.comments{6} = '# bicep3/analysis_logbook/20160611_b3_2016_beamparams';
        k.comments{7} = '# generated by script:';
        k.comments{8} = '# keck_analysis/beammap/ffbm_makebeamparamfiles.m';
        k.comments{9} = ['# generated ' datestr(clock)];
        k.filename = ['beams_keck_obs_' num2str(year) '0101.csv'];
        k.created = [datestr(clock,'yyyymmdd') ' KSK'];
        k.fields = {'gcp','r','theta','sigma','p','c','dp','dc',...
                    'aboffset_err_x','aboffset_err_y','err_sigma','err_p','err_c','err_dp','err_dc'};
        k.units = {'(null)','deg','deg','deg','unitless','unitless','unitless','unitless',...
                   'deg','deg','deg','unitless','unitless','unitless','unitless'};
        k.formats = {'integer','double','double','double','double','double','double','double',...
                     'double','double','double','double','double', ...
                     'double','double'};
      case 2017
        k.comments{1} = '# BICEP3 2017';
        k.comments{2} = '# Data from thermal beam map data taken Feb/Mar 2017';
        k.comments{3} = '# Since the beam maps do not contain absolute pointing information,';
        k.comments{4} = '# the pointing information (r,theta) is derived from adding AB pointing offsets to ideal pointing information';
        k.comments{5} = '# Analysis found here:';
        k.comments{6} = '# bkcmb/analysis_logbook/analysis/20180917_2017_beamparamcomparison';
        k.comments{7} = '# generated by script:';
        k.comments{8} = '# keck_analysis/beammap/ffbm_makebeamparamfiles.m';
        k.comments{9} = ['# generated ' datestr(clock)];
        k.filename = ['beams_keck_obs_' num2str(year) '0101.csv'];
        k.created = [datestr(clock,'yyyymmdd') ' TSG'];
        k.fields = {'gcp','r','theta','sigma','p','c','dp','dc',...
                    'aboffset_err_x','aboffset_err_y','err_sigma','err_p','err_c','err_dp','err_dc'};
        k.units = {'(null)','deg','deg','deg','unitless','unitless','unitless','unitless',...
                   'deg','deg','deg','unitless','unitless','unitless','unitless'};
        k.formats = {'integer','double','double','double','double','double','double','double',...
                     'double','double','double','double','double', ...
                     'double','double'};
    end % year switch
end % expt switch

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function kn = get_fpdata_comments(year,expt,type,filename,kk)

day = datestr(now,'yyyymmdd');

switch type
  case 'aboffsetflags'

    switch expt
      case 'keck'
        switch year
          case 2012
            kn.comments{1} = ['# Keck 2012'];
            kn.comments{2} = '# Following bicep2_aux_data/fp_data/fp_data_aboffsetflags_bicep2_20100105.csv';
            kn.comments{3} = '# Generated by ~clwong/20130703_beamparam/makefpdatabeamshapeflags.m';
            kn.comments{4} = '# Channel type code: L=light detector, D=dark detector, O=open input,,,,,,,,,,';
            kn.comments{5} = '# A/B offset magnitude is calculated from far field beam maps taken in Feb 2012'; 
            kn.comments{6} = '# detailed in http://bmode.caltech.edu/~spuder/keck_analysis_logbook/analysis/20130716_beamparam2012/';
            kn.comments{7} = '# A/B offset as a ratio of the keck beam sigma (0.2131 deg) is expressed in the column aboffset_pct; it has the same information but is merely scaled';
            kn.filename = ['aux_data/fp_data/' filename]; 
            kn.created = [day ' CLW'];
            kn.fields = [kk.fields(1:8) {'aboffset' 'aboffset_pct'}];
            kn.units = [kk.units(1:8) {'degrees' 'unitless'}];
            kn.formats = [kk.formats(1:8) {'double' 'double'}];
          case 2013
            kn.comments{1} = ['# Keck 2013'];
            kn.comments{2} = '# Following bicep2_aux_data/fp_data/fp_data_aboffsetflags_bicep2_20100105.csv';
            kn.comments{3} = '# Generated by ~clwong/20131112_channelcuts2013/makefpdatabeamshapeflags_2013.m';
            kn.comments{4} = '# Channel type code: L=light detector, D=dark detector, O=open input,,,,,,,,,,';
            kn.comments{5} = '# A/B offset magnitude is calculated from far field beam maps taken in Feb 2013'; 
            kn.comments{6} = '# detailed in http://bmode.caltech.edu/~spuder/keck_analysis_logbook/analysis/20131114_beamparam2013/';
            kn.comments{7} = '# A/B offset as a ratio of the keck beam sigma (0.2131 deg) is expressed in the column aboffset_pct; it has the same information but is merely scaled';
            kn.filename = ['aux_data/fp_data/' filename]; 
            kn.created = [day ' CLW'];
            kn.fields = [kk.fields(1:8) {'aboffset' 'aboffset_pct'}];
            kn.units = [kk.units(1:8) {'degrees' 'unitless'}];
            kn.formats = [kk.formats(1:8) {'double' 'double'}];
          case 2014
            kn.comments{1} = ['# Keck 2014'];
            kn.comments{2} = '# Following bicep2_aux_data/fp_data/fp_data_aboffsetflags_bicep2_20100105.csv';
            kn.comments{3} = '# Generated by keck_analysis/beammap/ffbm_makefpdata_beam.m';
            kn.comments{4} = '# Channel type code: L=light detector, D=dark detector, O=open input,,,,,,,,,,';
            kn.comments{5} = '# A/B offset magnitude is calculated from far field beam maps taken in Feb/Mar 2014'; 
            kn.comments{6} = '# detailed in http://bmode.caltech.edu/~spuder/keck_analysis_logbook/analysis/20140708_keck2014_beamparam/';
            kn.comments{7} = '# A/B offset as a ratio of the keck beam sigma (0.3042/0.2131 deg) is expressed in the column aboffset_pct; it has the same information but is merely scaled';
            kn.filename = ['aux_data/fp_data/' filename]; 
            kn.created = [day ' KSK'];
            kn.fields = [kk.fields(1:8) {'aboffset' 'aboffset_pct'}];
            kn.units = [kk.units(1:8) {'degrees' 'unitless'}];
            kn.formats = [kk.formats(1:8) {'double' 'double'}];
          case 2015
            kn.comments{1} = ['# Keck 2015'];
            kn.comments{2} = '# Following bicep2_aux_data/fp_data/fp_data_aboffsetflags_bicep2_20100105.csv';
            kn.comments{3} = '# Generated by keck_analysis/beammap/ffbm_makebeamparamfiles.m';
            kn.comments{4} = '# Channel type code: L=light detector, D=dark detector, O=open input,,,,,,,,,,';
            kn.comments{5} = '# A/B offset magnitude is calculated from far field beam maps taken in Feb/Mar 2015'; 
            kn.comments{6} = '# detailed in http://bmode.caltech.edu/~spuder/keck_analysis_logbook/analysis/20140708_keck2014_beamparam/';
            kn.comments{7} = '# A/B offset as a ratio of the keck beam sigma is expressed in the column aboffset_pct; it has the same information but is merely scaled';
            kn.filename = ['aux_data/fp_data/' filename]; 
            kn.created = [day ' KSK'];
            kn.fields = [kk.fields(1:8) {'aboffset' 'aboffset_pct'}];
            kn.units = [kk.units(1:8) {'degrees' 'unitless'}];
            kn.formats = [kk.formats(1:8) {'double' 'double'}];
          case 2016
            kn.comments{1} = ['# Keck 2016'];
            kn.comments{2} = '# Following bicep2_aux_data/fp_data/fp_data_aboffsetflags_bicep2_20100105.csv';
            kn.comments{3} = '# Generated by keck_analysis/beammap/ffbm_makebeamparamfiles.m';
            kn.comments{4} = '# Channel type code: L=light detector, D=dark detector, O=open input,,,,,,,,,,';
            kn.comments{5} = '# A/B offset magnitude is calculated from far field beam maps taken in Feb/Mar 2016'; 
            kn.comments{6} = '# detailed in http://bmode.caltech.edu/~spuder/keck_analysis_logbook/analysis/20140708_keck2014_beamparam/';
            kn.comments{7} = '# A/B offset as a ratio of the keck beam sigma is expressed in the column aboffset_pct; it has the same information but is merely scaled';
            kn.filename = ['aux_data/fp_data/' filename]; 
            kn.created = [day ' TSG'];
            kn.fields = [kk.fields(1:8) {'aboffset' 'aboffset_pct'}];
            kn.units = [kk.units(1:8) {'degrees' 'unitless'}];
            kn.formats = [kk.formats(1:8) {'double' 'double'}];	
          case 2017
            kn.comments{1} = ['# Keck 2017'];
            kn.comments{2} = '# Following bicep2_aux_data/fp_data/fp_data_aboffsetflags_bicep2_20100105.csv';
            kn.comments{3} = '# Generated by keck_analysis/beammap/ffbm_makebeamparamfiles.m';
            kn.comments{4} = '# Channel type code: L=light detector, D=dark detector, O=open input,,,,,,,,,,';
            kn.comments{5} = '# A/B offset magnitude is calculated from far field beam maps taken in Feb/Mar 2017'; 
            kn.comments{6} = '# detailed in bkcmb/analysis_logbook/analysis/20180917_2017_beamparamcomparison';
            kn.comments{7} = '# A/B offset as a ratio of the keck beam sigma is expressed in the column aboffset_pct; it has the same information but is merely scaled';
            kn.filename = ['aux_data/fp_data/' filename]; 
            kn.created = [day ' TSG'];
            kn.fields = [kk.fields(1:8) {'aboffset' 'aboffset_pct'}];
            kn.units = [kk.units(1:8) {'degrees' 'unitless'}];
            kn.formats = [kk.formats(1:8) {'double' 'double'}];
        end
        
      case 'bicep3'
        switch year
          case 2016
            kn.comments{1} = ['# BICEP3 2016'];
            kn.comments{2} = '# Following bicep2_aux_data/fp_data/fp_data_aboffsetflags_bicep2_20100105.csv';
            kn.comments{3} = '# Generated by keck_analysis/beammap/ffbm_makebeamparamfiles.m';
            kn.comments{4} = '# Channel type code: L=light detector, D=dark detector, O=open input,,,,,,,,,,';
            kn.comments{5} = '# A/B offset magnitude is calculated from far field beam maps taken in Feb/Mar 2016'; 
            kn.comments{6} = '# '; %'# detailed in http://bmode.caltech.edu/~spuder/keck_analysis_logbook/analysis/20140708_keck2014_beamparam/';
            kn.comments{7} = '# A/B offset as a ratio of the keck beam sigma is expressed in the column aboffset_pct; it has the same information but is merely scaled';
            kn.filename = ['aux_data/fp_data/' filename]; 
            kn.created = [day ' TSG'];
            kn.fields = [kk.fields([2:3 5:10]) {'aboffset' 'aboffset_pct'}];
            kn.units = [kk.units([2:3 5:10]) {'degrees' 'unitless'}];
            kn.formats = [kk.formats([2:3 5:10]) {'double' 'double'}];	
          case 2017
            kn.comments{1} = ['# BICEP3 2017'];
            kn.comments{2} = '# Following bicep2_aux_data/fp_data/fp_data_aboffsetflags_bicep2_20100105.csv';
            kn.comments{3} = '# Generated by keck_analysis/beammap/ffbm_makebeamparamfiles.m';
            kn.comments{4} = '# Channel type code: L=light detector, D=dark detector, O=open input,,,,,,,,,,';
            kn.comments{5} = '# A/B offset magnitude is calculated from far field beam maps taken in Feb/Mar 2017'; 
            kn.comments{6} = '# detailed in bkcmb/analysis_logbook/analysis/20180917_2017_beamparamcomparison';
            kn.comments{7} = '# A/B offset as a ratio of the keck beam sigma is expressed in the column aboffset_pct; it has the same information but is merely scaled';
            kn.filename = ['aux_data/fp_data/' filename]; 
            kn.created = [day ' TSG'];
            kn.fields = [kk.fields([2:3 5:10]) {'aboffset' 'aboffset_pct'}];
            kn.units = [kk.units([2:3 5:10]) {'degrees' 'unitless'}];
            kn.formats = [kk.formats([2:3 5:10]) {'double' 'double'}];
        end
    end
   
  case 'beamshapeflags'
    
    switch expt
      case 'keck'
        switch year
          case 2012
            kn.comments{1} = ['# Keck 2012'];
            kn.comments{2} = '# Following bicep2_aux_data/fp_data/fp_data_beamshapeflags_bicep2_20100105.csv';
            kn.comments{3} = '# Generated by ~clwong/20130703_beamparam/makefpdatabeamshapeflags.m';
            kn.comments{4} = '# Channel type code: L=light detector, D=dark detector, O=open input,,,,,,,,,,';
            kn.comments{5} = ['# sigma, p, and c are taken from the file keck_aux_data/beams/beams_keck_obs2012_rx*_' bfiledate '.csv'];
            kn.comments{6} = '# Beam parameters detailed in http://bmode.caltech.edu/~spuder/keck_analysis_logbook/analysis/20130716_beamparam2012/';
            kn.comments{7} = '# sigma - per detector beamwidth';
            kn.comments{8} = '# sigma_pct - per detector beamwidth divided by the nominal beamwidth (median sigma over rgls = .2131 degrees)';
            kn.comments{9} = '# sigma_diff - differential beamwidth';
            kn.comments{10} = '# sigma_diff_pct - differential beamwidth divided by nominal beamwidth';
            kn.comments{11} = '# ellip - per detector (p^2+c^2)^(1/2)';
            kn.comments{12} = '# ellip_diff - differential (p^2+c^2)^(1/2)';
            kn.comments{13} = '# pc_diff - quadrature sum of differential p and differential c ((p_a-p_b)^2 + (c_a-c_b)^2)^(1/2)';
            kn.filename = ['aux_data/fp_data/' filename]; 
            kn.created = [day ' CLW'];
            kn.fields = [kk.fields(1:8) {'sigma' 'sigma_pct' 'sigma_diff' 'sigma_diff_pct' 'ellip' 'ellip_diff' 'pc_diff'}];
            kn.units = [kk.units(1:8) {'degrees' 'unitless' 'unitless' 'unitless' 'unitless' 'unitless' 'unitless'}];
            kn.formats = [kk.formats(1:8) {'double' 'double' 'double'  'double' ...
                                'double' 'double' 'double'}];
          case 2013
            kn.comments{1} = ['# Keck 2013'];
            kn.comments{2} = '# Following bicep2_aux_data/fp_data/fp_data_beamshapeflags_bicep2_20100105.csv';
            kn.comments{3} = '# Generated by ~clwong/20131113_channelcuts2013/makefpdatabeamshapeflags_2013.m';
            kn.comments{4} = '# Channel type code: L=light detector, D=dark detector, O=open input,,,,,,,,,,';
            kn.comments{5} = ['# sigma, p, and c are taken from the file keck_aux_data/beams/beams_keck_obs2012_rx*_' bfiledate '.csv'];
            kn.comments{6} = '# Beam parameters detailed in http://bmode.caltech.edu/~spuder/keck_analysis_logbook/analysis/20131114_beamparam2013/';
            kn.comments{7} = '# sigma - per detector beamwidth';
            kn.comments{8} = '# sigma_pct - per detector beamwidth divided by the nominal beamwidth (median sigma over rgls = .2131 degrees)';
            kn.comments{9} = '# sigma_diff - differential beamwidth';
            kn.comments{10} = '# sigma_diff_pct - differential beamwidth divided by nominal beamwidth';
            kn.comments{11} = '# ellip - per detector (p^2+c^2)^(1/2)';
            kn.comments{12} = '# ellip_diff - differential (p^2+c^2)^(1/2)';
            kn.comments{13} = '# pc_diff - quadrature sum of differential p and differential c ((p_a-p_b)^2 + (c_a-c_b)^2)^(1/2)';
            kn.filename = ['aux_data/fp_data/' filename]; 
            kn.created = [day ' CLW'];
            kn.fields = [kk.fields(1:8) {'sigma' 'sigma_pct' 'sigma_diff' 'sigma_diff_pct' 'ellip' 'ellip_diff' 'pc_diff'}];
            kn.units = [kk.units(1:8) {'degrees' 'unitless' 'unitless' 'unitless' 'unitless' 'unitless' 'unitless'}];
            kn.formats = [kk.formats(1:8) {'double' 'double' 'double'  'double' ...
                                'double' 'double' 'double'}];
          case 2014
            kn.comments{1} = ['# Keck 2014'];
            kn.comments{2} = '# Following bicep2_aux_data/fp_data/fp_data_beamshapeflags_bicep2_20100105.csv';
            kn.comments{3} = '# Generated by ~keck_analysis/beammap/ffbm_makefpdata_beam.m';
            kn.comments{4} = '# Channel type code: L=light detector, D=dark detector, O=open input,,,,,,,,,,';
            kn.comments{5} = ['# sigma, p, and c are taken from the file keck_aux_data/beams/beams_keck_obs2012_rx*_20140705'  '.csv'];
            kn.comments{6} = '# Beam parameters detailed in http://bmode.caltech.edu/~spuder/keck_analysis_logbook/analysis/20140708_keck2014_beamparam/';
            kn.comments{7} = '# sigma - per detector beamwidth';
            kn.comments{8} = '# sigma_pct - per detector beamwidth divided by the nominal beamwidth (median sigma over rgls = .2131 degrees)';
            kn.comments{9} = '# sigma_diff - differential beamwidth';
            kn.comments{10} = '# sigma_diff_pct - differential beamwidth divided by nominal beamwidth';
            kn.comments{11} = '# ellip - per detector (p^2+c^2)^(1/2)';
            kn.comments{12} = '# ellip_diff - differential (p^2+c^2)^(1/2)';
            kn.comments{13} = '# pc_diff - quadrature sum of differential p and differential c ((p_a-p_b)^2 + (c_a-c_b)^2)^(1/2)';
            kn.filename = ['aux_data/fp_data/' filename]; 
            kn.created = [day ' KSK'];
            kn.fields = [kk.fields(1:8) {'sigma' 'sigma_pct' 'sigma_diff' 'sigma_diff_pct' 'ellip' 'ellip_diff' 'pc_diff'}];
            kn.units = [kk.units(1:8) {'degrees' 'unitless' 'unitless' 'unitless' 'unitless' 'unitless' 'unitless'}];
            kn.formats = [kk.formats(1:8) {'double' 'double' 'double'  'double' ...
                                'double' 'double' 'double'}];
          case 2015
            kn.comments{1} = ['# Keck 2015'];
            kn.comments{2} = '# Following bicep2_aux_data/fp_data/fp_data_beamshapeflags_bicep2_20100105.csv';
            kn.comments{3} = '# Generated by keck_analysis/beammap/ffbm_makebeamparamfiles.m';
            kn.comments{4} = '# Channel type code: L=light detector, D=dark detector, O=open input,,,,,,,,,,';
            kn.comments{5} = '# sig, p, c, dp, dc are taken from the file beammaps/beamparams/beamparams_2015.mat';
            kn.comments{6} = '# Beam parameters detailed in keck_analysis_logbook/analysis/20150602_keck_2015_beam_maps';
            kn.comments{7} = '# sigma: per-detector beamwidth';
            kn.comments{8} = '# sigma_pct: per detector beamwidth divided by the nominal beamwidth (median sigma over rgls)';
            kn.comments{9} = '# sigma_diff: diff beamwidth, sig(A) - sig(B)';
            kn.comments{10} = '# sigma_diff_pct: diff beamwidth divided by nominal beamwidth, (sig(A)-sig(B))/sig_med';
            kn.comments{11} = '# ellip: per-detector total ellipticity, (p^2 + c^2)^(1/2)';
            kn.comments{12} = '# ellip_diff: diff total ellipticity, (p(A)^2 + c(A)^2)^(1/2) - (p(B)^2 + c(B)^2)^(1/2)';
            kn.comments{13} = '# pc_diff: quadrature sum of diff p and diff c, ((p(A) - p(B))^2 + (c(A) - c(B))^2)^(1/2)';
            kn.filename = ['aux_data/fp_data/' filename]; 
            kn.created = [day ' KSK'];
            kn.fields = [kk.fields(1:8) {'sigma' 'sigma_pct' 'sigma_diff' ...
                                'sigma_diff_pct' 'ellip' 'ellip_diff' 'pc_diff'}];
            kn.units = [kk.units(1:8) {'degrees' 'unitless' 'unitless' ...
                                'unitless' 'unitless' 'unitless' 'unitless'}];
            kn.formats = [kk.formats(1:8) {'double' 'double' 'double' ...
                                'double' 'double' 'double' 'double'}];
          case 2016
            kn.comments{1} = ['# Keck 2016'];
            kn.comments{2} = '# Following bicep2_aux_data/fp_data/fp_data_beamshapeflags_bicep2_20100105.csv';
            kn.comments{3} = '# Generated by keck_analysis/beammap/ffbm_makebeamparamfiles.m';
            kn.comments{4} = '# Channel type code: L=light detector, D=dark detector, O=open input,,,,,,,,,,';
            kn.comments{5} = '# sig, p, c, dp, dc are taken from the file beammaps/beamparams/beamparams_2016.mat';
            kn.comments{6} = '# Beam parameters detailed in keck_analysis_logbook/analysis/';
            kn.comments{7} = '# sigma: per-detector beamwidth';
            kn.comments{8} = '# sigma_pct: per detector beamwidth divided by the nominal beamwidth (median sigma over rgls)';
            kn.comments{9} = '# sigma_diff: diff beamwidth, sig(A) - sig(B)';
            kn.comments{10} = '# sigma_diff_pct: diff beamwidth divided by nominal beamwidth, (sig(A)-sig(B))/sig_med';
            kn.comments{11} = '# ellip: per-detector total ellipticity, (p^2 + c^2)^(1/2)';
            kn.comments{12} = '# ellip_diff: diff total ellipticity, (p(A)^2 + c(A)^2)^(1/2) - (p(B)^2 + c(B)^2)^(1/2)';
            kn.comments{13} = '# pc_diff: quadrature sum of diff p and diff c, ((p(A) - p(B))^2 + (c(A) - c(B))^2)^(1/2)';
            kn.filename = ['aux_data/fp_data/' filename]; 
            kn.created = [day ' TSG'];
            kn.fields = [kk.fields(1:8) {'sigma' 'sigma_pct' 'sigma_diff' ...
                                'sigma_diff_pct' 'ellip' 'ellip_diff' 'pc_diff'}];
            kn.units = [kk.units(1:8) {'degrees' 'unitless' 'unitless' ...
                                'unitless' 'unitless' 'unitless' 'unitless'}];
            kn.formats = [kk.formats(1:8) {'double' 'double' 'double' ...
                                'double' 'double' 'double' 'double'}];
          case 2017
            kn.comments{1} = ['# Keck 2017'];
            kn.comments{2} = '# Following bicep2_aux_data/fp_data/fp_data_beamshapeflags_bicep2_20100105.csv';
            kn.comments{3} = '# Generated by keck_analysis/beammap/ffbm_makebeamparamfiles.m';
            kn.comments{4} = '# Channel type code: L=light detector, D=dark detector, O=open input,,,,,,,,,,';
            kn.comments{5} = '# sig, p, c, dp, dc are taken from the file beammaps/beamparams/beamparams_2017.mat';
            kn.comments{6} = '# Beam parameters detailed in bkcmb/analysis_logbook/analysis/20180917_2017_beamparamcomparison';
            kn.comments{7} = '# sigma: per-detector beamwidth';
            kn.comments{8} = '# sigma_pct: per detector beamwidth divided by the nominal beamwidth (median sigma over rgls)';
            kn.comments{9} = '# sigma_diff: diff beamwidth, sig(A) - sig(B)';
            kn.comments{10} = '# sigma_diff_pct: diff beamwidth divided by nominal beamwidth, (sig(A)-sig(B))/sig_med';
            kn.comments{11} = '# ellip: per-detector total ellipticity, (p^2 + c^2)^(1/2)';
            kn.comments{12} = '# ellip_diff: diff total ellipticity, (p(A)^2 + c(A)^2)^(1/2) - (p(B)^2 + c(B)^2)^(1/2)';
            kn.comments{13} = '# pc_diff: quadrature sum of diff p and diff c, ((p(A) - p(B))^2 + (c(A) - c(B))^2)^(1/2)';
            kn.filename = ['aux_data/fp_data/' filename]; 
            kn.created = [day ' TSG'];
            kn.fields = [kk.fields(1:8) {'sigma' 'sigma_pct' 'sigma_diff' ...
                                'sigma_diff_pct' 'ellip' 'ellip_diff' 'pc_diff'}];
            kn.units = [kk.units(1:8) {'degrees' 'unitless' 'unitless' ...
                                'unitless' 'unitless' 'unitless' 'unitless'}];
            kn.formats = [kk.formats(1:8) {'double' 'double' 'double' ...
                                'double' 'double' 'double' 'double'}];
        end
        
      case 'bicep3'
        switch year
          case 2016
            kn.comments{1} = ['# BICEP3 2016'];
            kn.comments{2} = '# Following bicep2_aux_data/fp_data/fp_data_beamshapeflags_bicep2_20100105.csv';
            kn.comments{3} = '# Generated by keck_analysis/beammap/ffbm_makebeamparamfiles.m';
            kn.comments{4} = '# Channel type code: L=light detector, D=dark detector, O=open input,,,,,,,,,,';
            kn.comments{5} = '# sig, p, c, dp, dc are taken from the file beammaps/beamparams/beamparams_2016.mat';
            kn.comments{6} = '# '; %'# Beam parameters detailed in keck_analysis_logbook/analysis/';
            kn.comments{7} = '# sigma: per-detector beamwidth';
            kn.comments{8} = '# sigma_pct: per detector beamwidth divided by the nominal beamwidth (median sigma over rgls)';
            kn.comments{9} = '# sigma_diff: diff beamwidth, sig(A) - sig(B)';
            kn.comments{10} = '# sigma_diff_pct: diff beamwidth divided by nominal beamwidth, (sig(A)-sig(B))/sig_med';
            kn.comments{11} = '# ellip: per-detector total ellipticity, (p^2 + c^2)^(1/2)';
            kn.comments{12} = '# ellip_diff: diff total ellipticity, (p(A)^2 + c(A)^2)^(1/2) - (p(B)^2 + c(B)^2)^(1/2)';
            kn.comments{13} = '# pc_diff: quadrature sum of diff p and diff c, ((p(A) - p(B))^2 + (c(A) - c(B))^2)^(1/2)';
            kn.filename = ['aux_data/fp_data/' filename]; 
            kn.created = [day ' TSG'];
            kn.fields = [kk.fields([2:3 5:10]) {'sigma' 'sigma_pct' 'sigma_diff' ...
                                'sigma_diff_pct' 'ellip' 'ellip_diff' 'pc_diff'}];
            kn.units = [kk.units([2:3 5:10]) {'degrees' 'unitless' 'unitless' ...
                                'unitless' 'unitless' 'unitless' 'unitless'}];
            kn.formats = [kk.formats([2:3 5:10]) {'double' 'double' 'double' ...
                                'double' 'double' 'double' 'double'}];
          case 2017
            kn.comments{1} = ['# BICEP3 2017'];
            kn.comments{2} = '# Following bicep2_aux_data/fp_data/fp_data_beamshapeflags_bicep2_20100105.csv';
            kn.comments{3} = '# Generated by keck_analysis/beammap/ffbm_makebeamparamfiles.m';
            kn.comments{4} = '# Channel type code: L=light detector, D=dark detector, O=open input,,,,,,,,,,';
            kn.comments{5} = '# sig, p, c, dp, dc are taken from the file beammaps/beamparams/beamparams_2017.mat';
            kn.comments{6} = '# Beam parameters detailed in bkcmb/analysis_logbook/analysis/20180917_2017_beamparamcomparison';
            kn.comments{7} = '# sigma: per-detector beamwidth';
            kn.comments{8} = '# sigma_pct: per detector beamwidth divided by the nominal beamwidth (median sigma over rgls)';
            kn.comments{9} = '# sigma_diff: diff beamwidth, sig(A) - sig(B)';
            kn.comments{10} = '# sigma_diff_pct: diff beamwidth divided by nominal beamwidth, (sig(A)-sig(B))/sig_med';
            kn.comments{11} = '# ellip: per-detector total ellipticity, (p^2 + c^2)^(1/2)';
            kn.comments{12} = '# ellip_diff: diff total ellipticity, (p(A)^2 + c(A)^2)^(1/2) - (p(B)^2 + c(B)^2)^(1/2)';
            kn.comments{13} = '# pc_diff: quadrature sum of diff p and diff c, ((p(A) - p(B))^2 + (c(A) - c(B))^2)^(1/2)';
            kn.filename = ['aux_data/fp_data/' filename]; 
            kn.created = [day ' TSG'];
            kn.fields = [kk.fields([2:3 5:10]) {'sigma' 'sigma_pct' 'sigma_diff' ...
                                'sigma_diff_pct' 'ellip' 'ellip_diff' 'pc_diff'}];
            kn.units = [kk.units([2:3 5:10]) {'degrees' 'unitless' 'unitless' ...
                                'unitless' 'unitless' 'unitless' 'unitless'}];
            kn.formats = [kk.formats([2:3 5:10]) {'double' 'double' 'double' ...
                                'double' 'double' 'double' 'double'}];
        end
    end
end

return