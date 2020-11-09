% load_arc_daily.m
%
% load_arc wrapper to load specified registers one day at a time,
% and save to .mat files.
%
% SYNTAX:
% load_arc_daily(ti,tf,rs,ds,fd)
%
% INPUTS:
% ti --- Date string of starting day (INCLUSIVE), 'yyyy-mmm-dd'
% tf --- Date string of ending day (INCLUSIVE), 'yyyy-mmm-dd'
% rl --- Cell array of GCP register string names to load 
%        (time is included by default)
% rs --- Cell array of GCP register string names to save
% ns --- Cell array of string names to assign to registers
% ds --- Downsample factor to downsample timestream data
% fd --- Directory to save the .mat files to
% nf --- No fast time register saved if == 1 
%
% OUTPUTS:
% Saves .mat file in format 'yyyymmdd'. 
%
% EXAMPLE:
% rl = { 'antenna0.compressor.h_alp[0]', ...};
% rs = { 'antenna0.compressor.h_alp', ...};
% ns = { 'alp', ...};
% fd = '/n/home09/jgrayson/analysis/jag_analysis/compressor/data/';
% load_arc_daily( '2015-Feb-01', '2015-Feb-01', rl, rs, ns, 10, fd, 0)
%
% 20150330 JAG: created. 

function load_arc_daily( ti, tf, rl, rs, ns, ds, fd, nf)

if ~exist('nf','var')
  nf=0
end

% Print to screen:
disp( 'load_arc_daily:')
disp(['  Loading data from ' ti ' to ' tf])
disp( '  For variables:')
for i = 1:length(rl)
  disp(['    ' char(rl(i))])
end
disp(['  with downsampling factor of n = ' num2str(ds)])
disp(['  to directory: ' fd])

% Arc file path
path = 'arc';
% Load start and end days as date numbers
firstday = datenum( ti, 'yyyy-mmm-dd');
lastday = datenum( tf, 'yyyy-mmm-dd');
if lastday<firstday
  error('Requested last day occurs before requested first day.')
end

% Add time to register list
slow = 'antenna0.time.utcslow';
fast = 'antenna0.time.utcfast';
if nf == 1
  rl = [rl slow];
  rs = [rs slow];
  ns = [ns 'mjd_slow'];
elseif nf == 0
  rl = [rl slow fast];
  rs = [rs slow fast];
  ns = [ns 'mjd_slow' 'mjd_fast'];
end
if length(rs)~=length(ns)
  error('Size of GCP register cell arrays do not match.')
end

% Load and save data one day at a time
for i = firstday:1:lastday

  % Get start and end time in load_arc format
  si = datestr( i, 'yyyy-mmm-dd:HH:MM:SS');
  sf = datestr( i+1, 'yyyy-mmm-dd:HH:MM:SS'); % ends 1 day later

  % Load the arc data for this day
  disp(['Loading arc day: ' si])
  d = load_arc( path, si, sf, rl);

  % Check if data is empty from load_arc
  if isempty(fieldnames(d))
    disp('No data returned from load_arc ...')
    continue
  end

  % Save variables in structure 'D' and downsample
  D = struct;
  for j = 1:length(rs)
    sr = char(rs(j));
    sn = char(ns(j));
    eval(['D.' sn '= d.' sr ';']) % save registers to variable names 
    eval(['D.' sn '= downsample( D.' sn ',' num2str(ds) ');']) % downsample array
  end
  % Convert mjd days and seconds to mjd days
  if ~isempty(D.mjd_slow)
    D.mjd_slow = D.mjd_slow(:,1) + D.mjd_slow(:,2)./86400;
    if nf ~= 1
      D.mjd_fast = D.mjd_fast(:,1) + D.mjd_fast(:,2)./86400;
    end
  end
  clear d

  % Save to file
  fn = datestr( i, 'yyyymmdd');
  save( [fd fn], '-struct', 'D');
  disp(['  Saving arc day: ' si])

end

return
