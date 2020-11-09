function [ch, min_dist] = rps_find_channels(scan, p, rpsopt, coverage_plot)
% [ch, min_dist] = rps_find_channels(scan, p, rpsopt, coverage_plot)
%
% Analyzes the mount trajectory during an rps scan and returns a
% list of channels that are on source.
%
% NOTE: Current version of this function doesn't know anything about which 
% receivers are off the ffflat, assumes all beams are on the mirror.
%
% [Arguments]
%   scan    Data structure describing a single rps scan. See rps_log_read.m
%   p       Data structure with detector properties, returned by
%           get_array_info.m
%   rpsopt  Data structure controlling options for RPS analysis.
%   coverage_plot  Set this optional argument to produce a plot showing 
%                  the scan coverage in (r,theta) coordinates and which 
%                  detectors were on source. If you set this argument 
%                  to a string like 'myplot.eps' or 'myplot.png', then
%                  the figure will be saved to the desired file name and 
%                  format (only eps and png supported).
%
%   Fields for rpsopt data structure:
%     rpsopt.mount     Mount model data structure. Default value supplied by
%                      keck_beam_map_pointing.m
%     rpsopt.mirror    Mirror model data structure. Default value supplied
%                      by keck_beam_map_pointing.m
%     rpsopt.source    Source position data structure. Default value supplied
%                      by keck_beam_map_pointing.m
%     rpsopt.min_dist  Minimum angular distance that the detector must 
%                      achieve relative to the source to select that
%                      channel. Defaults to 0.5 degrees.
%     rpsopt.flip_dk   Flips the sign of dk encoder offset for BICEP2 RPS 
%                      schedules that use mirror online pointing model.
%
% [Returns]
%   ch        List of detectors that are on source for this scan. Indices 
%             run across multiple receivers, so channel 529 indicates the 
%             first channel of rx1, etc.
%   min_dist  For all detectors, the minimum angular offset from the source 
%             achieved during the scan.

% Last update: 2014-05-07 CAB

% Check optional arguments.
if nargin < 4
  coverage_plot = [];
end

% Fill out rpsopt structure with default parameters, if necessary.
if ~exist('rpsopt', 'var')
  rpsopt = [];
end
% Use default mount model from keck_beam_map_pointing.m
if ~isfield(rpsopt, 'mount')
  rpsopt.mount = [];
end
% Use default mirror model from keck_beam_map_pointing.m
if ~isfield(rpsopt, 'mirror')
  rpsopt.mirror = [];
end
% Use default source location from keck_beam_map_pointing.m
if ~isfield(rpsopt, 'source')
  rpsopt.source = [];
end
% Minimum distance for selecting detectors defaults to 0.4 degrees.
if ~isfield(rpsopt, 'min_dist')
  rpsopt.min_dist = 0.4;
end
% For BICEP2, set the flip_dk keyword to 1 so that it correctly
% handles the online pointing model that was used for RPS observations.
if ~isfield(rpsopt, 'flip_dk')
  rpsopt.flip_dk = 0;
end

% Read in the scan pattern for one RPS angle.
data = load_arc('arc', mjd2datestr(scan.t1), ...
                mjd2datestr(scan.t2), {'antenna0'});
% If flip_dk keyword is set, flip the sign for the dk axis parts of
% encoder_off.
if rpsopt.flip_dk
  data.antenna0.tracker.encoder_off(:,3) = ...
      -1 * data.antenna0.tracker.encoder_off(:,3);
end

expt = get_experiment_name();



% Schedule encoder_zeros -52.621, -7.0, 2.626       
% Normal encoder_zeros 127.37899, -82.41876, -2.626 # 2016 Feb 2 - BICEP3 from 20160122
% Schdule encoder_cals: 2304000, -2304000, -574400, ant=all # az,el,dk counts per turn, reverse el&dk direction for ffflat 
% Normal encoder_cals: 2304000, 2304000, 574400

data.antenna0.tracker.encoder_mul(1,:) = [2304000, 2304000, 574400];
data.antenna0.tracker.encoder_off = repmat(([127.37899, -82.41876, -2.626])*3.6e6, ...
    size(data.antenna0.tracker.encoder_off,1),1);

% Apply offline pointing model.
if strcmp(expt, 'bicep3')
    pm = get_pointing_model(scan.t1,0, data);
    pm.az_tilt_ha = 0;
    pm.az_tilt_lat = 0;
    pm.el_tilt = 0;
	data = invpointing_model(data, pm);
	az = data.pointing.hor.az;
	el = data.pointing.hor.el;
	dk = data.pointing.hor.dk;
else
	pm = get_pointing_model(scan.t1, 0, data);
	data = invpointing_model(data, pm, 0);
	az = data.pointing.hor.az;
	el = data.pointing.hor.el;
	dk = data.pointing.hor.dk;
end


% Calculate boresight trajectory relative to source.
% This is a shortcut because we don't really need to calculate all
% of the detector trajectories (which is slow!).
bs.r = 0;
bs.theta = 0;
bs.chi = 0;
bs.chi_thetaref = 0;
bs.drumangle = 0;
mount.aperture_offr = 0; % Locate boresight 'aperture' in center of drum.
[r, theta, psi] = keck_beam_map_pointing(az, el, dk, rpsopt.mount, ...
                                         rpsopt.mirror, rpsopt.source, bs);

% Calculate rx, ry coordinates for boresight trajectory.
% Projection used here doesn't matter much, because we are looking
% for points where the detector pointing approaches very close to
% the source.
rx = 2 * sind(r / 2) .* cosd(theta) * 180 / pi;
ry = 2 * sind(r / 2) .* sind(theta) * 180 / pi;

% Calculate rx, ry coordinates for detectors.
prx = 2 * sind(p.r / 2) .* cosd(p.theta) * 180 / pi;
pry = 2 * sind(p.r / 2) .* sind(p.theta) * 180 / pi;

%%%%%%%%%%%%%%%%
% This was the old way. I didn't like it because it'd often pick
% detectors whose beams were just barely off coverage 
%which makes for crappy fits.
%
% Loop over detectors to find which ones hit the source in this scan.
%for i=1:length(p.r)
%  distance = sqrt((rx - prx(i)).^2 + (ry - pry(i)).^2);
%  min_dist(i) = min(distance);
%end
%ch = find(min_dist < rpsopt.min_dist);
%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%
% The new way creates a bounding polygon around the covered area
% and finds beam centers that fall within that area. 
% To ensure that most of the beam is within the area, I determine
% first which channels are within the boundary and then cut the channels
% that are too close determined by rpsopt.min_dist.

md=rpsopt.min_dist;
vi = convhull(rx,ry);

%Find indices that are inside and outside bounding polygon
[in on]=inpolygon(prx,pry,rx(vi),ry(vi));

ch = find(in|on)';
min_dist = zeros(size(ch));
for i = 1:length(ch)
  distance = sqrt((rx(vi) - prx(ch(i))).^2 + (ry(vi) - pry(ch(i))).^2);
  min_dist(i) = min(distance);
end

ch = ch(min_dist > md);

% Optional coverage plot.
if ~isempty(coverage_plot)
  % Default file name.
  plot_prefix = ['rps_coverage_' mjd2datestr(scan.t1, 'yyyymmdd_HHMMSS')];
  
  % Plot formats: 
  %   0: eps (default)
  %   1: png
  plot_format = 0;

  % Can specify prefix and/or format by passing a string to the
  % coverage plot argument.
  if isstr(coverage_plot)
    % Is .eps specified?
    tokens = regexp(coverage_plot, '(.*).eps', 'tokens');
    if ~isempty(tokens)
      plot_format = 0; % eps
      plot_prefix = tokens{1}{1};
    end

    % Is .png specified?
    tokens = regexp(coverage_plot, '(.*).png', 'tokens');
    if ~isempty(tokens)
      plot_format = 1; % png
      plot_prefix = tokens{1}{1};
    end
  end

  % Make a map of scan coverage.
  range = 15; % Covers [-range, range] for rx, ry.
  %binsize = 0.2; % Size of square pixels.
  %x = [-range:binsize:range];
  %y = x;
  %coverage = grid_map(rx, ry, ones(size(rx)), x, y);
  %coverage(isnan(coverage)) = 0;

  % Set up figure.
  fig = figure('Visible', 'off');
  set(fig, 'Position', [0 0 600 600], 'PaperPosition', [0 0 6 6]);

  % Draw the coverage region.
  %contour(x, y, coverage, 1, 'LineWidth', 2, 'Color', [0.5 0.5 0.5]);
  plot(rx(vi),ry(vi),'r')
  legend('coverage')
  hold all;

  % Plot detector positions. Highlight the selected detectors in red.
  plot(prx, pry, '.');
  plot(prx(ch), pry(ch), 'r.');

  % Tidy up the plot.
  xlim([-range, range]);
  ylim([-range, range]);
  xlabel('r_x');
  ylabel('r_y');
  axis square;

  % Print figure.
  if plot_format == 0
    figname = [plot_prefix '.eps'];
    print(figname, fig, '-depsc');
  elseif plot_format == 1
    figname = [plot_prefix '.png'];
    print(figname, fig, '-dpng');
  end
  close(fig);
end
