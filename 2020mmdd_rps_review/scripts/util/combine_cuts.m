function [cm,msgout]=combine_cuts(c,ind)
% cm=combine_cuts(c,ind)
%
% Combine together the cut masks in structure c to make over all cut
% mask cm
%
% ind is used to show statistics over "good" channels when printing stats

if(exist('ind','var'))
  disp('combining cuts...');
end

fields=fieldnames(c);
fields=fields(~ismember(fields,{'nhs','nch','nrx','nmce'}));

% combine together reporting stats at each step
i=1;
for f=fields'
  f=f{1};
  
  x=getfield(c,f);

  if(~exist('cm','var'))
    cm=true(size(x));
  end
  
  % logical and into overall mask
  cm=cm&x;

  % print stats
  if(exist('ind','var'))
    msg=sprintf('cut %s pass fraction %4.4f - aggregate so far %4.4f',f,...
      sum(rvec(x(:,ind.rgl)))/numel(x(:,ind.rgl)),...
      sum(rvec(cm(:,ind.rgl)))/numel(cm(:,ind.rgl)));
    disp(msg);
    msgout{i}=msg;
    i=i+1;
  end
end

return
