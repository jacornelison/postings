function temporaljack_splitdate=get_evenweight_splitdate(realmap,tag)
% temporaljack_splitdate = get_evenweight_splitdate(realmap, tag)
% returns the date in realmap where the pairdiff weights accumulated to half
% of the total pairdiff weight.  The format of temporaljack_splitdate is YYYYMMDD.
% This function does not distinguish between frequencies.
% 
% tag is passed to get_array_info.  If not specified, the first tag in the coaddopt
% of realmap is used.
%
% Example: 
% temporaljack_splitdate = get_evenweight_splitdate('maps/1351/real_d_filtp3_weight3_gs_dp0000_jack01.mat','20140404')
% was used for full season 2014 Keck data,
% and it returns 20140717
% To use this date, set:
% coaddopt.temporaljack_splitdate = '20140717';

x=load(sprintf('%s',realmap),'coaddopt');

if ~exist('tag','var') || isempty(tag)
  tag=x.coaddopt.tags{1};
end

% get array info
[p ind]=get_array_info(tag);

if numel(x.coaddopt.hsmax)>1
  hsmax=structcat(1,x.coaddopt.hsmax);
end

% This is the total pairdiff weight in realmap
wtot=nansum(nansum(hsmax.w(:,ind.lb)));

runsum=0;

% Silly loop to figure out when the running sum of weights exceeds half the total
for k=1:size(hsmax.w,1)
  runsum=runsum+nansum(hsmax.w(k,ind.lb));
  if runsum >= 0.5*wtot
    break
  end
end

splittag=x.coaddopt.tags{k};
temporaljack_splitdate=splittag(1:8);

return

end
