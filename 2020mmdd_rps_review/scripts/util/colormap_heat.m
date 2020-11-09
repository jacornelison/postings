function heat = colormap_heat(n,rev)
% white-red colormap.  use like this:
% heat = colormap_heat();
% colormap(heat)

% Improvement on matlab's heat colormap. This is 1/2 of the lint colormap
% that you can grab using colormap_lint

% Optional arguments:
% n: resolution (default 128)
% rev: reverse colorscale

% v1. 05.20.2013. rwa


if ~exist('n','var') || isempty(n)
  n=128;
end

d = [247 247 247;...
  253 219 199;...
  244 165 130;... 
  214 96 77;...
  178 24 43;...
  103 0 31;...
  1 1 1];
  
usa = d/255;

iv = ceil(128/7);

heat = interp1(1:iv:n,usa,1:n,'spline');
heat(heat<0)=0;
heat(heat>1)=1;

heat = flipud(heat);

if exist('rev','var')
  for ii=1:3
    heat(:,ii) = wrev(heat(:,ii));
  end
end