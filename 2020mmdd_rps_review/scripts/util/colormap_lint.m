function lint = colormap_lint(n,rev)
% red-white-blue colormap.  use like this:
% lint = colormap_lint();
% colormap(lint)

% Improvement on matlab's jet colormap.  ripped off of colorscale 
% that I grabbed from colorbrewer2.org

% Optional arguments:
% n: resolution (default 128)
% rev: reverse colorscale

% v1. 05.20.2013. rwa


if ~exist('n','var') || isempty(n)
  n=128;
end

d = [103 0 31 ;...
  178 24 43;...
  214 96 77;...
  244 165 130;...
  253 219 199;...
  247 247 247;...
  209 229 240;...
  146 197 222;...
  67 147 195;...
  33 102 172;...
  5 48 97];

usa = d/255;

iv = ceil(n/11);

lint = interp1(1:iv:n,usa,1:n,'spline');
lint(lint<0)=0;

lint = flipud(lint);

if exist('rev')
  for ii=1:3
    lint(:,ii) = wrev(lint(:,ii));
  end
end
