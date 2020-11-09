function [m,map]=upsample_map(m,map,upsample_factor)

% simple interpolation routine for plotting purposes only
% 'map' is a single panel for plot_map, not the full map structure

if ~exist('upsample_factor','var')
  upsample_factor=2;
end

upsample_factor=max(1,ceil(upsample_factor));

% interpolate the tic
dx=mean(diff(m.x_tic));
dy=mean(diff(m.y_tic));
x=linspace(...
  m.x_tic(1)-(upsample_factor-1)*dx,...
  m.x_tic(end)+(upsample_factor-1)*dx,...
  numel(m.x_tic)*(2*upsample_factor-1));
y=linspace(...
  m.y_tic(1)-(upsample_factor-1)*dy,...
  m.y_tic(end)+(upsample_factor-1)*dy,...
  numel(m.y_tic)*(2*upsample_factor-1));

% interpolate the map
[rin,cin]=meshgrid(m.x_tic,m.y_tic);
[rout,cout]=meshgrid(x,y);
map(isnan(map))=0;
map=interp2(rin,cin,map,rout,cout);
map(map==0)=nan;

m.x_tic=x;
m.y_tic=y;