function hnull(lims,opt);
% h=hnull(lims);
%
% Make an empty plot frame
%
% If opt='S' superimpose invisible axes on
% current figure without rescaling.

if(~exist('lims','var'))
  lims=[];
end
if(~exist('opt','var'))
  opt=' ';
end

if(isempty(lims))
  lims=[0,10,0,10];
end

if(any(opt=='S'))
  pos=get(gca,'Position');
  axes('Position',pos,'Xlim',lims(1:2),'Ylim',lims(3:4));
  set(gca,'Visible','off');
  return
end

axis(lims);
set(gca,'Box','on');
