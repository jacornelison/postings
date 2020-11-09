function aps=aps_bin_merge(aps,n)
% function aps=aps_bin_merge(aps,n)
%
% Merge aps bins to downgrade resolution
%  
%080513: dB adds modification so 1st fine bin is not counted in 1st coarse bin ( but not fully implemented)

for j=1:length(aps)
  
  m=size(aps(j).Cs_l);
  %deal with first l bin 
 % aps(j).Cs_l(1,:,:)=squeeze(mean(reshape(aps(j).Cs_l(2:n,:),[n-1,1,m(2:end)]),1));
 % aps(j).l(1)=mean(reshape(aps(j).l(2:n),[n-1,1]),1)';

 %deal with rest of bins
  aps(j).Cs_l=squeeze(mean(reshape(aps(j).Cs_l,[n,m(1)/n,m(2:end)]),1));
  aps(j).l=mean(reshape(aps(j).l,[n,m(1)/n]),1)';
end

% if expv (expectation value) value present take mean also
if(isfield(aps,'expv'))
  for j=1:length(aps)
    m=size(aps(j).expv);
    
    %deal with first bin
  %  aps(j).expv(1,:,:)=squeeze(mean(reshape(aps(j).expv(2:n,:),[n-1,1,m(2:end)]),1));
    
    %rest of bins
    aps(j).expv=squeeze(mean(reshape(aps(j).expv,[n,m(1)/n,m(2:end)]),1));
  end
end

return
