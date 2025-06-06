function ffbm_makesinglebeamfile(year,bfiledate)
% ffbm_makesinglebeamfile(year,bfiledate)
% 
% Take the individual rx beam files for Keck and concatenate them
% 
% INPUTS
%         year:      2014
%         bfiledate: '20140710'
%
% KSK 2014-07-10

% File to save as
filename = ['beams_keck_obs_' num2str(year) '0101.csv']

% Date of individual beams files
% bfiledate = '20140710';

for ii = 1:5
  [p{ii},k{ii}] = ParameterRead(['beamfiles/beams_keck_obs' num2str(year) '_rx' num2str(ii-1) '_' bfiledate '.csv'])
end

pn = structcat(1,[p{1},p{2},p{3},p{4},p{5}]);

switch year
  case 2014
    kn.comments{1}=['# KECK 2014 all receivers'];
    kn.comments{2}='# Data from thermal beam map data taken Feb/Mar 2014';
    kn.comments{3}='# This file is a concatenation of individual receiver beams file';
    kn.comments{4}='# Since the beam maps do not contain absolute pointing information,';
    kn.comments{5}='# the pointing information (r,theta) is derived from adding AB pointing offsets to ideal pointing information';
    kn.comments{6}='# Analysis found here:';
    kn.comments{7}='# ~/spuder/keck_analysis_logbook/analysis/20140708_keck2014_beamparam';
    kn.comments{8}='# generated by script:';
    kn.comments{9}='# keck_analysis/beammap/ffbm_makesinglebeamfile.m';
    kn.comments{10}=['# generated ' datestr(clock)];
    kn.filename=filename;
    kn.created=[datestr(clock,'yyyymmdd') ' KSK'];
    kn.fields={'gcp' 'r' 'theta' 'sigma' 'p' 'c' 'aboffset_err_x' 'aboffset_err_y' 'err_sigma','err_p','err_c'};
    kn.units={'(null)' 'deg' 'deg' 'deg' 'unitless' 'unitless' 'deg' 'deg' 'deg' 'unitless' 'unitless'};
    kn.formats={'integer' 'double' 'double' 'double' 'double' 'double' ...
	'double' 'double' 'double' 'double' 'double'};
    
    ParameterWrite(['beamfiles/' filename],pn,kn)
end

return