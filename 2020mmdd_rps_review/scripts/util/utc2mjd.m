% UTC2MJD  Extract the Modified Julian Day number from
%   GCP time registers, such as antenna0.time.utcslow.
%
% RWO 090518
function m = utc2mjd (u)

if (nargin < 1) || isempty (u)
	error (['Not enough input arguments.']);
end;

if ~isnumeric (u) || ~(isa (u, 'double') || isa (u, 'uint64'))
        error (['UTC times should be doubles or 64-bit unsigned integers.']);
end;

% 64-bit uint, not unpacked at all
if isa (u, 'uint64')
        % First cast to uint32, to separate date and time
        tt = typecast (u, 'uint32');
        dd = double (tt (1:2:end));
        tt = double (tt (2:2:end)) / 1000;

        % return date part
        m = dd;
        return;
end;

% Doubles, unpacked into MJD + time in seconds
if isa (u, 'double')
        if size(u,2)==2
                m = u(:,1);
        elseif size(u,1)==2
                m = u(1,:);
        else
                error (['UTC times should have two columns, MJD and time of day.']);
        end;
end;

