function str = mjd2datestr(mjd, format)
% str = mjd2datestr(mjd)
%
% Convert decimal MJD values to strings.
%
% Example: mjd2datestr(55000.2678) returns '2009-Jun-18:06:25:37'.
% 
% [Arguments]
%   mjd     Modified Julian Date.
%   format  Optionally specify an alternate format for the date string. 
%           Default format is 'yyyy-mmm-dd:HH:MM:SS'.
%
% [Returns]
%   str     Formatted date string.

% 2014-02-06 CAB

% Default format.
if nargin < 2
  format = 'yyyy-mmm-dd:HH:MM:SS';
end

% MJD=0 corresponds to 1858-Nov-17:00:00:00.
str = datestr(mjd + datenum(1858, 11, 17), format);
