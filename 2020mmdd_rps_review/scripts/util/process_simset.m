function apsset=process_simset(apsset)
% function apsset=process_simset(apsset)
%
% Calc mean, std, eom, and inp for a set of aps from simulation
% provided in apsset

% calc mean, std and eom
for j=1:numel(apsset)
  for i=1:size(apsset(j).Cs_l,2)
    apsset(j).mean(:,i)=nanmean(apsset(j).Cs_l(:,i,:),3);
    apsset(j).std(:,i)=nanstd(apsset(j).Cs_l(:,i,:),0,3);
  end
  apsset(j).eom=apsset(j).std/sqrt(size(apsset(j).Cs_l,3));
  apsset(j).eos=apsset(j).std/sqrt(2*size(apsset(j).Cs_l,3));
end

return
