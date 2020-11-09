% LIST_ARC_FILES  Find arc files in a directory
%   matching a specified range of UTC times.
%
%   FL=LIST_ARC_FILES(D,T1,T2,FEXT,ONEBEFORE)
%      D    = directory containing data files
%      T1   = starting UTC time
%      T2   = ending UTC time
%      FEXT = file extension pattern
%             (optional; default is .dat*)
%      ONEBEFORE = toggle last-before-t1 option
%             (optional; default is true)
%
%   Note that LIST_ARC_FILES returns not only
%   files with file names falling between T1
%   and T2, but also the last file before T1.
%   This is because a file starting before T1
%   may contain frames extending past T1. Set
%   ONEBEFORE=false to turn this behavior off.

function fl = list_arc_files (arc_dir, t1, t2, fext, onebefore)

	% default extension is .dat*, for arc files
	% (compressed or otherwise)
	if nargin<4 || isempty(fext)
		fext = '.dat*';
	end
	if fext(1)~='.'
		fext=['.' fext];
	end

	% get starting and ending UTC times
        if (nargin < 2) || isempty(t1)
                t1 = '1900-jan-1';
        end;
        if (nargin < 3) || isempty(t2)
                t2 = '3000-jan-1';
        end;
        utc1 = txt2utc (t1);
        utc2 = txt2utc (t2);

	% by default, onebefore=true
	if (nargin < 5) || isempty(onebefore)
		onebefore = true;
	end

	% Two ways to get the directory listings.
	% When requesting 75 days or fewer, list the files
	% day by day.  When requesting more, get a listing
	% for the entire directory at once.  This generally
	% gives good performance.

        % tic
        if utc2-utc1 > 75 
	    % disp('Listing all arc files');
            d = dir (fullfile (arc_dir, ['*' fext]));
            fl = {d(:).name};
        else
	    % disp('Listing selected days');
            fl = {};
            for utc_day = floor(utc1-1):floor(utc2)
		    d = dir (fullfile (arc_dir, [datestr(utc_day - 48987 + datenum(1993,1,1) - 1,'yyyymmdd') '_*' fext]));
		    fl = {fl{:}, d(:).name};
            end
        end
        % toc

	% Now select the files that are actually within the
	% time range of interest.
	cc = true (size (fl));
	uu = -1 * ones(size(fl));
	for (ifl = 1:length(d))
		if (fl{ifl}(1) ~= '2')
			cc(ifl) = false;
		end;
	end;
	if isempty({fl{cc}})
		fl = {};
		return
	end
	dl = datenum ({fl{cc}}, 'yyyymmdd_HHMMSS');
	uu(cc) = dl + 48987 - datenum (1993,1,1) + 1;
        [uu ifl] = sort (uu);
        % uu = uu(ifl);
        fl = {fl{ifl}};

	cc = (uu >= utc1) & (uu <= utc2);
	cc = cc(:);

	% if requested, do the ONEBEFORE thing
	if onebefore
		if sum(cc)==0
			% If no files were between t1 and t2,
			% Return the last file before t1.
			% It presumably contains data running
			% through t1 and t2.
			tdiff = utc1 - uu;
			tdiff (tdiff < 0) = Inf;
			[tdiff idx] = min (tdiff);
			cc(idx) = 1;
		else
			% Keep one file before t1, since it may contain
			% data running through t1.
			cc = (cc | [cc(2:end); false]);
		end
	end

	fl = {fl{cc}};
