% UTC2DATENUM  Convert GCP time registers, such as
%   antenna0.time.utcslow, to Matlab datenum format
%   (time since January 1, year zero, in days and
%   fractions of days).
%
% Can also be used to convert UTC time strings,
%   such as '2009-aug-13' or '2009-aug-13:01:14:05',
%   to Matlab datenums.
%
% RWO 090518
function m = utc2datenum (u)

if (nargin < 1) || isempty (u)
	error (['Not enough input arguments.']);
end;

% String like '2009-aug-13:01:14:05'
if ischar (u)
	FMT_LIST = {'dd-mmm-yyyy:HH:MM:SS','yyyy-mmm-dd:HH:MM:SS', 'yyyy-mmm-dd:HH:MM', 'yyyy-mmm-dd:HH', ...
		    'yyyy-mmm-dd', ...
		    'yyyy-mm-dd:HH:MM:SS', 'yyyy-mm-dd:HH:MM', 'yyyy-mm-dd:HH', ...
		    'yyyy-mm-dd'};
	for (ii = 1:length(FMT_LIST))
           disp(FMT_LIST{ii});
		try,
            % Specify pivot year=1 to avoid inconsistent attempts in various Matlab
            % versions to coax out-of-range values into plausible dates.
            m = datenum (u, FMT_LIST{ii}, 1);
            if m<7E5, error(' '); end
            return;
		catch,
            % Wrong format, try another
		end;
	end;
	error (['Unrecognized UTC format ' u '.']);
end;

if (nargin < 1) || isempty (u)
        error (['Not enough input arguments.']);
end;

if ~isnumeric (u) || ~(isa (u, 'double') || isa (u, 'uint64'))
        error (['UTC times should be doubles or 64-bit unsigned integers.']);
end;

clear dd tt

% 64-bit uint, not unpacked at all
if isa (u, 'uint64')
        % First cast to uint32, to separate date and time
        tt = typecast (u, 'uint32');
        dd = double (tt (1:2:end));
        tt = double (tt (2:2:end)) / 1000;
end;

% Doubles, unpacked into MJD + time in seconds
if isa (u, 'double')
        if size(u,2)==2
                dd = u(:,1);
		tt = u(:,2);
        elseif size(u,1)==2
                dd = u(1,:);
		tt = u(2,:);
        else
                error (['UTC times should have two columns, MJD and time of day.']);
        end;
end;

if ~exist('dd','var') || ~exist ('tt','var')
	error (['Unrecognized UTC time format.']);
end;

% convert to days + fractions of days
m = dd + tt / (24 * 3600);

% Change offset from MJD zero to year-zero
% http://tycho.usno.navy.mil/mjd.html
% Note that in the Navy document on MJD,
% the Day of Year (DOY) is from 1, not 0
m = m - (48987+1) + datenum (1993, 1, 1);
