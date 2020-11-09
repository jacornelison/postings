function map=smooth_varmaps(m,map,width)
% map=smooth_varmaps(m,map,width)
%
% We mask the maps using 1/var map
% If var map contains high freq power this mixes modes
% By smoothing the var map we can reduce this effect at the cost of a
% slightly sub-optimal weighting

% sigma in pixels of Gaussian smoothing kernel
% - note that the roll down to zero for an edge step occurs over
% approx 6 times this range rather than 3
if ~exist('width','var')
  width=[];
end
if isempty(width)
  width=0.5; % in deg
end

% convert width to pixels
width=width/m.pixsize;

% make the grid
x=-5*width:+5*width;
[xg,yg]=meshgrid(x,x);
      
% generate the smoothing kernel
g=egauss2([1,0,0,width,width,0],xg,yg);
g=g./sum(g(:));

% for each map
for i=1:numel(map)
  % smooth the maps
  map(i).Tvar=smooth_map(map(i).Tvar,g,width);
  
  if isfield(map,'Tdifvar')
    map(i).Tdifvar=smooth_map(map(i).Tdifvar,g,width);
  end

  if isfield(map,'Dvar')
    map(i).Dvar=smooth_map(map(i).Dvar,g,width);
  end

  if isfield(map(i),'Q')
    map(i).Qvar=smooth_map(map(i).Qvar,g,width);
    map(i).Uvar=smooth_map(map(i).Uvar,g,width);
  end

end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%
function m=smooth_map(m,g,width)

% find the nan region
n=isnan(m);

% do the smoothing in mask (1/var) space
m=1./m;

% set nan's to zero
m(n)=0;

% whittle back from the edge by 5 smoothing widths
% - this is so the taper down to edge becomes smooth
for i=1:5*width
  m(bwperim(m))=0;
end

% convolve
m=conv2(m,g,'same');

% set nan region back to nan
m(n)=NaN;

% go back to var
m=1./m;

return
