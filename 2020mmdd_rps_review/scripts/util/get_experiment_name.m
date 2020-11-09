function expt=get_experiment_name(tag)
% expt=get_experiment_name(tag)
%
% Returns the experiment name, e.g. bicep2, keck, bicep3

if(exist('tag','var')&&~isempty(tag)&&isletter(tag(1)))
  % unified style - if tag arg provided and first char is letter
  % use this to determine which expt we are
  switch tag(1:2)
   case 'B2'
    expt='bicep2';
   case 'KA'
    expt='keck';
   case 'B3'
    expt='bicep3';
  end
  
else
  % traditional style - get expt name from file in aux_data
  
  % Read in from a file in aux_data
  tmp=textread('aux_data/experiment_name.txt','%s');
  
  % Check for empty result
  if length(tmp)==0
    expt='';
  else
    
    % Skip any comment lines
    for i=1:length(tmp)
      expt=strtrim(tmp{i});
      if length(expt)>0 && expt(1)~='%'
        break
      end
    end
  end
end
  
return
