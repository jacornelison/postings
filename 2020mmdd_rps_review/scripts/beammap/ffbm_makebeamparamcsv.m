function ffbm_makebeamparamcsv(experiment,year,rxNum)
% function ffbm_makebeamparamcsv(experiment,year,rxNum)
%
% Make the beam parameter .csv files for BICEP2 and Keck
%
% Takes the beam parameter summary files output from ffbm_param_manyfiles.
% The .csv should then be copied and committed to aux_data/beams.
%
% Note that if the beamwidth isn't measured, we set it to a fiducial
% beamwidth which is defined at the top
%
% INPUTS
%        experiment:    'bicep2','keck'
%        year:           2012/2013 for keck, rewrite for keck 2014
%                        2011/2012/'both' for bicep2
%        rxNum:          0/1/2/3/4 only relevant for keck
%                        'all' will make it for all receivers.
%
% OUTPUTS
%         .csv is saved in ./beamfiles (make sure this exists)
%
% CLW 2014-06-16

medsig_100 = 0.3042;
medsig_150 = 0.2131;

% Get array info
switch experiment
  case 'bicep2'
    switch year
      case 2011
	[p ind] = get_array_info('20110201');
      case 2012
	[p ind] = get_array_info('20120201');
      otherwise
	[p ind] = get_array_info('20120201');
    end
  case 'keck'
    switch year
      case 2012
	[p ind] = get_array_info('20120201');
      case 2013
	[p ind] = get_array_info('20130201');
      case 2014
	[p ind] = get_array_info('20140201');
    end
    if ~strcmp(rxNum,'all')
      cutind = (p.rx == rxNum);
      p = structcut(p,cutind);
      ind = make_ind(p);
    end
    p.theta2 = p.theta;
    p.theta = p.theta2 + p.drumangle;
end

% Choose input beam parameter files
switch experiment
  case 'bicep2'
    prestr = '~clwong/20140310_B2beamparam/'
    load([prestr 'beamparam/beamparam_bicep2_' num2str(year)]);
  case 'keck'
    switch year
      case 2012
	prestr = '~clwong/20130703_beamparam/';
      case 2013
	prestr = '~clwong/20131112_channelcuts2013/';
      case 2014
	prestr = './'; % Assuming you're in keckpipe
    end
    if strcmp(rxNum,'all')
      load([prestr 'beamparam_' num2str(year) '/beamparam_all']);
    else
      load([prestr 'beamparam_' num2str(year) '/beamparam_rx' num2str(rxNum)]);
    end
end

% Load up the mask
switch experiment
  case 'bicep2'
    switch year
      case 'both'
	numFiles = 24;
      otherwise
	M = load([prestr 'bicep2_mask_' num2str(year) '/mask_b2.mat'])
	mask = M.mask;
	numFiles = length(mask);
    end
  case 'keck'
    if ~strcmp(rxNum,'all')
      M = load([prestr 'mask_ffbm_' num2str(year) '/mask_rx' num2str(rxNum) '.mat']);
      mask = M.mask;
      numFiles = length(mask);
    else
      switch year
	case 2012
	  numFiles = 5;
	case 2013
	  numFiles = 14;
	case 2014
	  numFiles = '';
      end
    end
end    

% Process the beam parameters
params.sig_mean = nanmedian(sig,2);
params.sig_std = nanstd(sig,[],2);
params.e_p_mean = nanmedian(e_p,2);
params.e_p_std = nanstd(e_p,[],2);
params.e_c_mean = nanmedian(e_c,2);
params.e_c_std = nanstd(e_c,[],2);
params.x = nanmedian(x,2);
params.y = nanmedian(y,2);

dx(ind.e,numFiles) = NaN;
dy(ind.e,numFiles) = NaN;
dx(ind.a,:) = x(ind.a,:) - x(ind.b,:);
dy(ind.a,:) = y(ind.a,:) - y(ind.b,:);

% Move back to r, theta coordinates
[theta,r] = cart2pol(params.x,params.y);
params.r = r;
params.theta = theta * 180.0/pi;
params.aboffset_err_x(ind.e,1) = NaN;
params.aboffset_err_y(ind.e,1) = NaN;
params.aboffset_err_x(ind.a,1) = nanstd(dx(ind.a,:),[],2);
params.aboffset_err_x(ind.b,1) = params.aboffset_err_x(ind.a);
params.aboffset_err_y(ind.a,1) = nanstd(dy(ind.a,:),[],2);
params.aboffset_err_y(ind.b,1) = params.aboffset_err_y(ind.a)

% Find non-NaN values
cc = find(~isnan(params.aboffset_err_x));
pr.gcp = p.gcp;
pr.r = p.r;
pr.theta = p.theta;
pr.r(cc) = params.r(cc);
pr.theta(cc) = params.theta(cc);
pr.aboffset_err_x = NaN(528,1);
pr.aboffset_err_y = NaN(528,1);
pr.aboffset_err_x(cc) = params.aboffset_err_x(cc);
pr.aboffset_err_y(cc) = params.aboffset_err_y(cc);
thing = find(pr.aboffset_err_x == 0);
pr.aboffset_err_x(thing) = NaN;
pr.aboffset_err_y(thing) = NaN;

cc = find(~isnan(params.sig_mean));
switch p.band(1)
  case 100
    pr.sigma = medsig_100 * ones(528,1);
  case 150
    pr.sigma = medsig_150 * ones(528,1); 
end
pr.sigma(cc) = params.sig_mean(cc);
pr.c = zeros(528,1);
pr.c(cc) = params.e_c_mean(cc);
pr.p = zeros(528,1);
pr.p(cc) = params.e_p_mean(cc);

pr.err_sigma = zeros(528,1);
pr.err_sigma(cc) = params.sig_std(cc);
pr.err_c = zeros(528,1);
pr.err_c(cc) = params.e_c_std(cc);
pr.err_p = zeros(528,1);
pr.err_p(cc) = params.e_p_std(cc);

% NaNs for dark SQUIDs
pr.sigma(ind.o) = NaN;
pr.c(ind.o) = NaN;
pr.p(ind.o) = NaN;
pr.aboffset_err_x(ind.o) = NaN;
pr.aboffset_err_y(ind.o) = NaN;
pr.err_sigma(ind.o) = NaN;
pr.err_c(ind.o) = NaN;
pr.err_p(ind.o) = NaN;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Comments for the .csv
switch experiment
  case 'bicep2'
    switch year
      case 2011
	k.comments{1}=['# Bicep2'];
	k.comments{2}='# Data from thermal beam map data taken Nov/Dec 2011';
	k.comments{3}='# Since the beam maps do not contain absolute pointing information,';
	k.comments{4}='# the pointing information (r,theta) is derived from adding AB pointing offsets to ideal pointing information';
	k.comments{5}='# Analysis found here:';
	k.comments{6}='# ';
	k.comments{7}='# generated by script:';
	k.comments{8}='# ~clwong/20140310_B2beamparam/makebeamparamcsv.m';
	k.comments{9}=['# generated ' datestr(clock)];
	k.filename=['beams_bicep2_' year '.csv'];
	k.created=[datestr(clock,'yyyymmdd') ' CLW'];
	k.fields={'gcp' 'r' 'theta' 'sigma' 'p' 'c' 'aboffset_err_x' 'aboffset_err_y' 'err_sigma','err_p','err_c'};
	k.units={'(null)' 'deg' 'deg' 'deg' 'unitless' 'unitless' 'deg' 'deg' 'deg' 'unitless' 'unitless'};
	k.formats={'integer' 'double' 'double' 'double' 'double' 'double' ...
	    'double' 'double' 'double' 'double' 'double'};
      case 2012
	k.comments{1}=['# Bicep2'];
	k.comments{2}='# Data from thermal beam map data taken Nov 2012';
	k.comments{3}='# Since the beam maps do not contain absolute pointing information,';
	k.comments{4}='# the pointing information (r,theta) is derived from adding AB pointing offsets to ideal pointing information';
	k.comments{5}='# Analysis found here:';
	k.comments{6}='# ';
	k.comments{7}='# generated by script:';
	k.comments{8}='# ~clwong/20140310_B2beamparam/makebeamparamcsv.m';
	k.comments{9}=['# generated ' datestr(clock)];
	k.filename=['beams_bicep2_' num2str(year) '.csv'];
	k.created=[datestr(clock,'yyyymmdd') ' CLW'];
	k.fields={'gcp' 'r' 'theta' 'sigma' 'p' 'c' 'aboffset_err_x' 'aboffset_err_y' 'err_sigma','err_p','err_c'};
	k.units={'(null)' 'deg' 'deg' 'deg' 'unitless' 'unitless' 'deg' 'deg' 'deg' 'unitless' 'unitless'};
	k.formats={'integer' 'double' 'double' 'double' 'double' 'double' ...
	    'double' 'double' 'double' 'double' 'double'};
      otherwise
	k.comments{1}=['# Bicep2'];
	k.comments{2}='# Data from thermal beam map data taken 2011/2012';
	k.comments{3}='# Since the beam maps do not contain absolute pointing information,';
	k.comments{4}='# the pointing information (r,theta) is derived from adding AB pointing offsets to ideal pointing information';
	k.comments{5}='# Analysis found here:';
	k.comments{6}='# ';
	k.comments{7}='# generated by script:';
	k.comments{8}='# ~clwong/20140310_B2beamparam/makebeamparamcsv.m';
	k.comments{9}=['# generated ' datestr(clock)];
	k.filename=['beams_bicep2_' num2str(year) '.csv'];
	k.created=[datestr(clock,'yyyymmdd') ' CLW'];
	k.fields={'gcp' 'r' 'theta' 'sigma' 'p' 'c' 'aboffset_err_x' 'aboffset_err_y' 'err_sigma','err_p','err_c'};
	k.units={'(null)' 'deg' 'deg' 'deg' 'unitless' 'unitless' 'deg' 'deg' 'deg' 'unitless' 'unitless'};
	k.formats={'integer' 'double' 'double' 'double' 'double' 'double' ...
	    'double' 'double' 'double' 'double' 'double'};
    end % year switch

  case 'keck'
    switch year
      case 2012
	k.comments{1}=['# KECK RX' num2str(rxNum) 'FPU'];
	k.comments{2}='# Data from thermal beam map data taken Feb 2012';
	k.comments{3}='# Since the beam maps do not contain absolute pointing information,';
	k.comments{4}='# the pointing information (r,theta) is derived from adding AB pointing offsets to ideal pointing information';
	k.comments{5}='# Analysis found here:';
	k.comments{6}='# ~/spuder/keck_analysis_logbook/analysis/20130716_beamparam2012';
	k.comments{7}='# generated by script:';
	k.comments{8}='# ~clwong/20130703_beamparam/param_manyfiles.m';
	k.comments{9}=['# generated ' datestr(clock)];
	k.filename=['beams_keck_obs' num2str(year) '_rx' num2str(rxNum) '_20130703.csv'];
	k.created=[datestr(clock,'yyyymmdd') ' CLW'];
	k.fields={'gcp' 'r' 'theta' 'sigma' 'p' 'c' 'aboffset_err_x' 'aboffset_err_y' 'err_sigma','err_p','err_c'};
	k.units={'(null)' 'deg' 'deg' 'deg' 'unitless' 'unitless' 'deg' 'deg' 'deg' 'unitless' 'unitless'};
	k.formats={'integer' 'double' 'double' 'double' 'double' 'double' ...
	    'double' 'double' 'double' 'double' 'double'};	
      case 2013
	k.comments{1}=['# KECK RX' num2str(rxNum) 'FPU'];
	k.comments{2}='# Data from thermal beam map data taken Feb 2013';
	k.comments{3}='# Since the beam maps do not contain absolute pointing information,';
	k.comments{4}='# the pointing information (r,theta) is derived from adding AB pointing offsets to ideal pointing information';
	k.comments{5}='# Analysis found here:';
	k.comments{6}='# ~/spuder/keck_analysis_logbook/analysis/20131114_beamparam2013';
	k.comments{7}='# generated by script:';
	k.comments{8}='# ~clwong/20131112_channelcuts2013/param_manyfiles.m';
	k.comments{9}=['# generated ' datestr(clock)];
	k.filename=['beams_keck_obs' num2str(year) '_rx' num2str(rxNum) '_20131113.csv'];
	k.created=[datestr(clock,'yyyymmdd') ' CLW'];
	k.fields={'gcp' 'r' 'theta' 'sigma' 'p' 'c' 'aboffset_err_x' 'aboffset_err_y' 'err_sigma','err_p','err_c'};
	k.units={'(null)' 'deg' 'deg' 'deg' 'unitless' 'unitless' 'deg' 'deg' 'deg' 'unitless' 'unitless'};
	k.formats={'integer' 'double' 'double' 'double' 'double' 'double' ...
	    'double' 'double' 'double' 'double' 'double'};	
      case 2014
	day=datestr(now,'yyyymmdd');
	k.comments{1}=['# KECK RX' num2str(rxNum) 'FPU'];
	k.comments{2}='# Data from thermal beam map data taken Feb/Mar 2014';
	k.comments{3}='# Since the beam maps do not contain absolute pointing information,';
	k.comments{4}='# the pointing information (r,theta) is derived from adding AB pointing offsets to ideal pointing information';
	k.comments{5}='# Analysis found here:';
	k.comments{6}='# ~/spuder/keck_analysis_logbook/analysis/20140708_keck2014_beamparam';
	k.comments{7}='# generated by script:';
	k.comments{8}='# keck_analysis/beammap/ffbm_makebeamparamcsv.m';
	k.comments{9}=['# generated ' datestr(clock)];
	k.filename=['beams_keck_obs' num2str(year) '_rx' num2str(rxNum) ...
	      '_' day '.csv'];
	k.created=[datestr(clock,'yyyymmdd') ' KSK'];
	k.fields={'gcp' 'r' 'theta' 'sigma' 'p' 'c' 'aboffset_err_x' 'aboffset_err_y' 'err_sigma','err_p','err_c'};
	k.units={'(null)' 'deg' 'deg' 'deg' 'unitless' 'unitless' 'deg' 'deg' 'deg' 'unitless' 'unitless'};
	k.formats={'integer' 'double' 'double' 'double' 'double' 'double' ...
	    'double' 'double' 'double' 'double' 'double'};
    end % year switch
end % experiment switch

% Write out the .csv
switch experiment
  case 'bicep2'
    switch year
      case 2011
	ParameterWrite(['beamfiles/beams_b2_obs' num2str(year) ...
	      '.csv'],pr,k);
      case 2012
	ParameterWrite(['beamfiles/beams_b2_obs' num2str(year) ...
	      '.csv'],pr,k);
      otherwise
	ParameterWrite(['beamfiles/beams_b2_obs.csv'],pr,k);
    end
  case 'keck'
    switch year
      case 2012
	ParameterWrite(['beamfiles/beams_keck_obs' num2str(year) '_rx' ...
	      num2str(rxNum) '_20130703.csv'],pr,k);	
      case 2013
	ParameterWrite(['beamfiles/beams_keck_obs' num2str(year) '_rx' ...
	      num2str(rxNum) '_20131113.csv'],pr,k); 	
      case 2014
	ParameterWrite(['beamfiles/beams_keck_obs' num2str(year) '_rx' ...
	      num2str(rxNum) '_' day '.csv'],pr,k);
    end
end

return


