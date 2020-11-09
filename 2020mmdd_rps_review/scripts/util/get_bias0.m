function bias0=get_bias0(yr,expt)
%%%
% bias0=get_bias0(yr,expt) is a helper function for reduc_applycal where yr is the year
% of the tag as a string. expt is an optional argument.  If not present, it's set with
% get_experiment_name. Acceptable inputs for expt are 'keck' and 'bicep3'
%
% With the changes to the way tod's are calibrated in 2014, we multiply the elnod_gain
% of each detector by its tes bias.  We then rescale by the so-called bias0 factor to
% minimize the effect of this multiplication
%
% In this function, bias0 should have as many values as there are unique frequencies in
% the tod.  They should be listed in frequency-ascending order

if ~exist('yr','var')
  error('A year must be used as an input argument to this function')
end

% make sure yr is a string
if ~ischar(yr)
  error('Please input the year as a string')
end


% get experiment name, if necessary
if ~exist('expt','var') || isempty(expt)
  expt=get_experiment_name;
end


switch expt
 case 'keck'
  
  switch yr
   case '2014'
    bias0=[1000,2600];
   case '2015'
    bias0=[1000,2600,3300];
   case '2016'
    % 150, 210, 220 GHz
    bias0=[2600,2500,3300];
   case '2017'
    % 210, 220, 270 GHz
    % temp values
    bias0=[2500,3300,2500];
   case '2018' 
    bias0=[2500,3300,2217];
   case '2019' 
    %  150 (copied from 2015), 210, 220, 270 GHz
    bias0=[2600, 2500,3300,2217];
  end
  
 case 'bicep3'
  
  switch yr
   case '2015'
    bias0=240;
   case '2016' 
    bias0=270;
   case '2017'
    % temp value
    bias0=270;
   case '2018'
    bias0=270;
  end
  
end

  
