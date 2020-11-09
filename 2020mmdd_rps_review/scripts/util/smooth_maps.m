function map=smooth_maps(m,map,width)
% map=smooth_maps(m,map)
%
%same as smooth_varmaps but smooths the normal map (T,Q,U)
 
if(~exist('width','var'))
  width=[];
end
if isempty(width)
  width=1.0; % in deg
end
 
% convert width to pixels
width=width/m.pixsize;

% make the grid
tic=-5*width:+5*width;
[xg,yg]=meshgrid(tic,tic);

% generate the smoothing kernel
g=egauss2([1,0,0,width,width,0],xg,yg);
g=g./sum(g(:));

% for each map
for i=1:numel(map)
   
  % smooth the maps
  map(i).T=smooth_map(map(i).T,g,width);
  if isfield(map(i),'Q')
    map(i).Q=smooth_map(map(i).Q,g,width);
    map(i).U=smooth_map(map(i).U,g,width);
  end
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%
function m=smooth_map(m,g,width)

% find the nan region
n=isnan(m);

% do the smoothing in mask (1/var) space
%m=1./m;

% set nan's to zero
m(n)=0;

% whittle back from the edge by 3 smoothing widths
% - this is so the taper down to edge becomes smooth
for i=1:3*width
  m(bwperim(m))=0;
end

% convolve
m=conv2(m,g,'same');

% set nan region back to nan
m(n)=NaN;

% go back to var
%m=1./m;

return
