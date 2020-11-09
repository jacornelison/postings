function [aps,mod]=strip_aps_lowl(aps,mod,type)
% function [aps,mod]=strip_aps_lowl(aps,mod,type)
%
% Remove high l points from an aps

% determine max l
switch type
  case 'quad05'
    %max_l=aps(1).l(end)+1; % already limited in makeaps
    % I messed up in makecorrmap - sims only contain power up to
    % ell=2000 so must limit here for now
    max_l=1936; % leaves 96 bins
  case 'spud'
    max_l=150;
end

% strip aps
ind=aps(1).l<max_l;
for j=1:length(aps)
  aps(j).l=aps(j).l(ind);
  aps(j).Cs_l=aps(j).Cs_l(ind,:,:);
  
  if(isfield(aps(j),'mean'))
    aps(j).mean=aps(j).mean(ind,:);
    aps(j).std=aps(j).std(ind,:);
    aps(j).eom=aps(j).eom(ind,:);
  end
end

% strip model
ind=mod.l<(max_l+100);
mod.l=mod.l(ind);
mod.Cs_l=mod.Cs_l(ind,:);

return
