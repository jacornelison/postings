function [h,s]=stairhist(data,edges,color,fcolor,donorm,ishist)
% [h,s]=stairhist(data,edges,color,fcolor,donorm,ishist)
% makes a stair histogram using histc() and stairs() internally
%  
% data  :   to  be histogrammed
% edges :   the bin edges, also see help of histc 
% color :   ('k') the edge color of the hist
% fcolor:   ([]) the face color, no
% donorm:   (0) convenience - (1): norm to integral unity, (2): norm to max
% ishist:   (0) data is already a histogram, just plot it
%
% Example:  stairhist(randn(1000,1,1),linspace(-4,4,50))

if(~exist('color','var')) color = 'k'; end
if(~exist('fcolor','var'))|| isempty(fcolor) fcolor= [] ; end
if(~exist('donorm','var'))|| isempty(donorm) donorm= 0  ; end
if(~exist('ishist','var'))|| isempty(ishist) ishist= 0  ; end

if ~ishist
  h = histc(data,edges);
else
  h = data;
end

if donorm==1
  h=h/length(data);
elseif donorm==2
  h=h/max(h);
end

ih = ishold;

s = [];
if ~isempty(fcolor)
  % patches don't clear the axis even if hold is off
  % do it explicitly:
  if ~ih cla; end
  
  % draw filled patches:
  pb = [];
  ph = [];
  for ii=1:length(edges)-1
    pb=[pb,edges(ii),edges(ii),edges(ii+1),edges(ii+1)];
    ph=[ph,0,h(ii),h(ii),0];
  end
  s = patch(pb,ph,fcolor,'EdgeColor','none');
  
  % if edges are to be drawn on top, switch on hold
  hold on;
end

% draw regular stairs
if ~isempty(color)
  s2 = stairs(edges,h);
  set(s2,'color',color);
  s = [s,s2];
end

% bring the axis back to old state:
if ih
  hold on;
else
  hold off;
end
  
return
