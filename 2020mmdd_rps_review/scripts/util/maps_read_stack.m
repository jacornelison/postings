function [mapu,mapd,tag]=maps_read_stack(dates,series)
% [mapu,mapd,tag]=maps_read_stack(ltags,stag)
%
% Read in a series of 1 day files

% read in all the maps
ctag=[];
for z=1:length(dates)
  dates{z}
  
  load(['data/',dates{z},series]);
  
  amapu(z,:,:)=mapu;
  amapd(z,:,:)=mapd;
  ctag=[ctag,'+',dates{z}];
end

mapu=amapu;
mapd=amapd;

% massage tag string
if(length(dates)>8)
  tag=sprintf('%s through %s',dates{1},dates{end});
else
  tag=ctag(2:end);
end
  
return
