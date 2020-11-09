function r=add_xfac(r)
% r=add_xfac(r)
%
% Add the "x factors" aka the noise offsets
% These are used in the offset-lognormal transformation to generate
% the right shape bandpower likelihood
%
% For auto spectra they are just the total that has been subtracted
% (debiased) out of the real data
% For cross spectra (or combined bandpower containing cross spectra)
% it's not clear what is the right thing to do
% see: http://find.spa.umn.edu/~quad/logbook/20080530/bp_uncer.pdf

for j=1:length(r)
  for i=1:size(r(j).sim,2)
    % recipe used for QUaD - supposed to work when combined including cross spectra
    %r(j).xfac(:,i)=mean(r(j).sig(:,i,:),3).*std(r(j).noi(:,i,:),[],3)./std(r(j).sig(:,i,:),[],3);

    % most basic recipe
    r(j).xfac(:,i)=mean(r(j).noi(:,i,:),3);
    
    % if r.db exists is record of all debiasing - noise and signal - so
    % use this
    if(isfield(r,'db'))
      r(j).xfac(:,i)=r(j).db(:,i);
    end
    
  end
end

return
