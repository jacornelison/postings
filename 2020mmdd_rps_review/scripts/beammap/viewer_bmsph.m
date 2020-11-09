function viewer_bmsph(map, channel, varargin)
% viewer_bmsph(map, channel, ...)
%
%   Draw beams from maps generated in beam map coordinates.
%     Defined in http://bmode.caltech.edu/~spuder/analysis_logbook/analysis/20131005_sidelobe_coordinates/sidelobe_coordinates.pdf
%
%   The maps are binned in x and y for a spherical coordinate
%   system with the boresight of the pixel at the north pole. 
%
%   Inputs:
%     map       Map data structure 
%     channel   Channel number to plot.
%   Keyword arguments:
%     'Log'     Specify this keyword to plot the map on a logarithmic
%               scale (base 10). The abs function is applied to the
%               map before taking the logarithm.
%     'CLim', [low high]
%               Specify this keyword to manually set the range of the
%               color scale. By default, the color scale will be set
%               to the full range of the map. Note that range
%               specified with the 'CLim' keyword is not transformed
%               based on the presence of the 'Log' keyword.
%     'NPix', value
%     'Graticule', [r_step, theta_step]
%               Specify this keyword to add a graticule to the
%               plot. By default, the graticule has lines spaced by 30
%               degrees in r and 60 degrees in theta. The line
%               spacing can be specified by following the keyword
%               with a vector argument containing the spacing in r
%               and theta.

% 2011-12-06 CAB
% 2014-04-28 IDB

% Parse optional arguments.
plotlog = 0;
clim = [];
npix = 201;
thetamax = 90;
graticule = 0;
graticule_r_step = 30;
graticule_theta_step = 60;
for i=1:length(varargin)
  % Logarithmic map.
  if strcmp(varargin{i}, 'Log')
    plotlog = 1;
  end
  
  % Select color scale.
  if strcmp(varargin{i}, 'CLim')
    value = double(varargin{i + 1});
    if (length(value) == 2) & (value(1) < value(2))
      clim = value;
      i = i + 1;
    end
  end

  % Select plot size.
  if strcmp(varargin{i}, 'NPix')
    value = double(varargin{i + 1});
    if (value > 0)
      npix = value;
      i = i + 1;
    end
  end
  
  % Select plot range.
  if strcmp(varargin{i}, 'ThetaMax')
    value = double(varargin{i + 1});
    if (value > 0) & (value < 90)
      thetamax = value;
      i = i + 1;
    end
  end

  % Plot graticule.
  if strcmp(varargin{i}, 'Graticule')
    graticule = 1;
    % Check if graticule spacing is specified.
    if (length(varargin) > i)
      if isvector(varargin{i + 1}) & (length(varargin{i + 1}) == 2)
	value = double(varargin{i + 1});
	graticule_r_step = value(1);
	graticule_theta_step = value(2);
      end
    end
  end
end

% Get input map, with logarithmic scaling if the log keyword is set.
zin = map.map(:, :, channel);
if plotlog
  zin = log10(abs(map.map(:, :, channel)));
  zin(isinf(zin)) = NaN;
end
  
% Default color scale is the full map range.
if isempty(clim)
  clim = [min(zin(:)) max(zin(:))];
end


% Draw projected image.
imagesc(map.x_bin, map.y_bin, zin, clim);
set(gca,'YDir','normal');
colorbar();
axis image;
axis square;
xlabel('x_P');
ylabel('y_P');

% Draw graticule.
%   30 degree steps in r.
%   60 degree steps in theta.
if graticule
  hold on;
  % Lines of constant r are circles.
  theta = [-pi:0.05:pi] *180/pi;
  for r=[graticule_r_step:graticule_r_step:180]
    [x, y] = rtheta_to_xbyb(r, theta);
    plot(x, y, 'k');
  end
  % Lines of constant theta are straight lines.
  r = [0 180];
  for theta=[0:graticule_theta_step:360]
    [x, y] = rtheta_to_xbyb(r, theta);
    plot(x, y, 'k');
  end
end
