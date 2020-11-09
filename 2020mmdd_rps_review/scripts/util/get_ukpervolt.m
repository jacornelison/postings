function ukpervolt=get_ukpervolt(tag,expt)
% ukpervolt=get_ukpervolt(tag,expt)
%
% keep the ukpervolt cal factor in one place to make it easy to change
% - the value comes from reduc_abscal
%
% tag = tag to determine season
% expt = bicep2, keck or bicep3
%        for keck seasons 2014 onward per rx values are returned to
%        deal with multi frequency

if(~exist('tag','var'))
  tag=[];
end

if(isempty(tag))
  tag='2012';
end

if ~exist('expt','var')||isempty(expt)
  expt=get_experiment_name;
end

switch expt
  case 'bicep2'
   
   switch tag(1:4)
     
    case '2010'
     % order of mag value
     % ukpervolt=[5000,5000];
     % preliminary abscal
     ukpervolt=3150;
     
    case '2011'
     % order of mag value
     % ukpervolt=[5000,5000];
     % preliminary abscal
     ukpervolt=3150;
     
    case '2012'
     % order of mag value
     % ukpervolt=[5000,5000];
     % preliminary abscal
     ukpervolt=3150;
     
    otherwise
     error(['Don''t have an ukpervolt value for year ' tag(1:4)]);
     
   end
   
 case 'keck'
  
  switch tag(1:4)

   case '2012'
    % http://bicep.caltech.edu/~spuder/keck_analysis_logbook/analysis/20130719_abscal/
    ukpervolt=3400;

   case '2013'
    % http://bicep.caltech.edu/~spuder/keck_analysis_logbook/analysis/20130124_abscal/
    ukpervolt=2900;

   case '2014'
    % first half season:
    % http://bicep.caltech.edu/~spuder/keck_analysis_logbook/analysis/20140718_keck2014_final-rx-abscals/
    % on per rx basis to cope with the different frequencies.
    %ukpervolt=[2021,2991,2026,3108,3018];

    % Full season analysis:
    % http://bicep.rc.fas.harvard.edu/keck/analysis_logbook/analysis/20141215_keck2014_rx-abscals/
    ukpervolt=[1942,2966,1926,3033,3065];

   case '2015'
    % All values updated according to analysis at
    % http://bicep.rc.fas.harvard.edu/bkcmb/analysis_logbook/analysis/20160404_K2015_abscal/
    ukpervolt=[1921,10947,1925,11269,2871];

   case '2016'
    % These values are based on the analysis at
    % http://bicep.rc.fas.harvard.edu/bkcmb/analysis_logbook/analysis/20180331_K2016_abscal_update/
    ukpervolt=[5826,10872,5925,11307,2724];
    
   case '2017'
     % Rx0,1,2,3 abscals are from the posting
     % http://bicep.rc.fas.harvard.edu/bkcmb/analysis_logbook/analysis/20180923_K2017_abscal/
     % The abscal of Rx4 (this Rx is excluded from BK17) is from this posting
     % http://bicep.rc.fas.harvard.edu/bkcmb/analysis_logbook/analysis/20180413_MSIP_proposal_K2017maps/
     ukpervolt=[5682,10518,5762,10963,16376];

   case '2018'
     % Preliminary
     % Posting: http://bicep.rc.fas.harvard.edu/bkcmb/analysis_logbook/analysis/20180527_K2018prelim_Performance/
     ukpervolt=[6354,11145,6249,11847,12382];

   case '2019'
     % Preliminary
     % Posting: http://bicep.rc.fas.harvard.edu/bkcmb/analysis_logbook/analysis/20180527_K2018prelim_Performance/
     % add a placeholder for rx1(smurf receiver)
     ukpervolt=[6354,10000,6249,11847,12382];

   otherwise
    error(['Don''t have an ukpervolt value for year ' tag(1:4)]);
    
  end
  
 case 'bicep3'
  switch tag(1:4)
   case {'2014'}
    % For test runs
    ukpervolt=3150; % Value copied from previous bicep2
   case {'2015'}
    % http://bicep.rc.fas.harvard.edu/bicep3/analysis_logbook/20150921_abscal/
    ukpervolt=3827;
   case {'2016'}
    % http://bicep.rc.fas.harvard.edu/bkcmb/analysis_logbook/analysis/20180322_B2016_abscal_poly/
    ukpervolt=4243;
   case {'2017'}
    % http://bicep.rc.fas.harvard.edu/bkcmb/analysis_logbook/analysis/20180921_B2017_abscal/
    ukpervolt=3246;
   case '2018'
    % From B2018 preliminary analysis
    % http://bicep.rc.fas.harvard.edu/bkcmb/analysis_logbook/analysis/20180608_B2018_prelim/
    ukpervolt=3249;
   case '2019'
    % Placeholder from B2018 values (need to update later)
    ukpervolt=3249;

   otherwise
    error(['Don''t have an ukpervolt value for year ' tag(1:4)]);

  end
end

return
