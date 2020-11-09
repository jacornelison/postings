function [C F SCH] = parse_run_log (logdir,t1,t2)
%[C F S]=parse_run_log(logdir,T1,T2)
%
% Read the GCP log files found in logdir, and
% parse out the observing phase start and end
% times.  The schedule names, phase codes and
% start and end times are output in structure
% C.  Fridge cycle information is returned in
% the structure F, and schedule start and end
% times in S, gleaned from the same log files.
%
% The optional input parameters T1 and T2 are
% the start and stop times to set a period of
% interest.  They can be specified as strings
% in UTC, like '2011-jan-1' or '2010-dec-1:10'.
%
% Output structures:
%
% C.t1,t2 = tag start and end time (UTC, text)
%  .sch   = schedule file name, if any
%  .t0    = date of schedule start
%  .scset = scanset / field name, if any
%  .ph    = phase letter
%  .ch    = number of tag within phase
%  .sc    = scan name, if any
%  .nsc   = number of scans within scanset
%  .dk    = deck angle
%  .elofs = elevation offset
%  .livetime = live time in minutes spent in:
%           1. Total
%           2. On source
%           3. El nods
%           4. Partial load curves
%           5. Fridge cycles
%           6. Sky dips
%           7. Full load curves
%  .nmark = number of completed operations
%           in same catgories as above
%
% F.t     = cycle start and stop (as datenums)
%  .ph    = phase letter (default is A)
%
% S.schedule = schedule name
%  .t        = start and stop times (datenums)
%

if (nargin < 1) || isempty(logdir)
	logdir = '/data/bicep2daq/log/';
end;
if (nargin < 2) || isempty(t1)
	t1 = '2010-jan-1';
end;
if (nargin < 3) || isempty(t2)
	t2 = '3000-dec-31';
end;

% Select arc files up to three days before start
% of time range.  This ensures we parse the schedule
% start, which is needed to build the tag names.
% This number will need to be tweaked if we ever
% go to a schedule lasting > 3 days.
d = list_arc_files (logdir, ...
	datestr(utc2datenum(t1)-3,'yyyy-mmm-dd:HH:MM:SS'), ...
	t2, '.log');

IN_FNAME = '';
IN_SCHED = '';
IN_PHASE = '';
IS_CRYO_PHASE = false;
IN_SCSET = '';
SCAN_NAME = '';
SCAN_NREPS = [];
ON_DK = [];
ON_EL = [];
EL_AT_START = [];

IN_CYCLE = 0;
IN_CAL = 0;

SCSET_IDX = 0;
PHCAL_IDX = 0;
CHUNK_IDX = 0;
CALDY_IDX = 0;

MARKBITS = clear_markbits();

ON_TIME = [];
ON_TIME_STR = '';
CURR_START_TIME = [];
SCHED_START_DATE = [];
LAST_CAL_START_DATE = [];
LAST_PHCAL_DATE_PH = [];

C = [];
F = [];
SCH = [];

for (ii = 1:length(d))
	if isempty (regexp (d{ii}, '^(\d){8}_(\d){6}.log'))
		continue;
	end;
	fname = fullfile (logdir, d{ii});
	f = fopen (fname, 'rt');
	IN_FNAME = f;
	disp (fname);

	while ~feof (f)
		ll = fgetl (f);
		if ~ischar(ll) || isempty(ll)
			continue;
		end;
		llcell = textscan (ll, '%s %s %[^\n]');
		if (length(llcell) < 3) || isempty(llcell{3})
			continue;
		end;
		tmp_date = llcell{1}{1};
		tmp_time = llcell{2}{1};
		ll = strtrim(llcell{3}{1});

		% Remember the string from the current line's time.
		% But don't convert to a datenum yet, because datestr
		% is slow, and most lines are irrelevant.
		ON_TIME_STR = [tmp_date ' ' tmp_time];
		ON_TIME = [];
		% disp ([tmp_date ' ' tmp_time ' ' datestr(ON_TIME) ' ' ll]);
		s = check_for_sched_start (ll);
		if ~isempty (s)
			% disp (['In schedule ' s]);
			if isempty(ON_TIME)
				ON_TIME = datenum (ON_TIME_STR, 'yymmdd HH:MM:SS');
			end;
                        if ~isempty(IN_SCHED) && ~isempty(IN_PHASE) && ~isempty(IN_SCSET)
                                if strcmp(IN_SCSET,'perphasecals')
                                        CHUNK_IDX = -1 * PHCAL_IDX;
					PHCAL_IDX = PHCAL_IDX + 1;
                                else
                                        CHUNK_IDX = SCSET_IDX;
                                end;
                                % print_chunk_line(CURR_START_TIME, ON_TIME, IN_SCHED, SCHED_START_DATE, IN_PHASE, CHUNK_IDX);
				C = add_chunk(C, CURR_START_TIME, ON_TIME, IN_SCHED, SCHED_START_DATE, IN_SCSET, IN_PHASE, CHUNK_IDX, SCAN_NAME, SCAN_NREPS, ON_DK, EL_AT_START, MARKBITS);
                                % disp ([datestr(CURR_START_TIME) ', ' datestr(ON_TIME) ', ' IN_SCHED ', ' SCHED_START_DATE IN_PHASE num2str(CHUNK_IDX)]);
                        end;
			if ~isempty(IN_SCHED) && length(SCH)>0
				SCH(end).t(2) = ON_TIME;
			end;
			SCH(end+1).schedule = s;
			SCH(end).t(1) = ON_TIME;
			IN_SCHED = s;
			tmp = datestr (ON_TIME, 'yyyymmdd');
			if length(LAST_PHCAL_DATE_PH)<9 || ~strcmp (tmp, LAST_PHCAL_DATE_PH(1:8))
				PHCAL_IDX = 0;
			end;
			SCHED_START_DATE = tmp;
			[tmp1 tmp2] = fileparts(s);
			IN_CAL = strncmp('7_',tmp2,2);
			if IN_CAL
				if isempty(LAST_CAL_START_DATE) || ~strcmp(LAST_CAL_START_DATE, SCHED_START_DATE)
					CALDY_IDX = 0;
				end;
				% CALDY_IDX = CALDY_IDX + 1;
				SCSET_IDX = CALDY_IDX;
				LAST_CAL_START_DATE = SCHED_START_DATE;
				IN_PHASE = '';
				IS_CRYO_PHASE = false;
				% disp (['In cals!']);
			else
				SCSET_IDX = 0;
				IN_PHASE = '';
				IS_CRYO_PHASE = false;
			end;
			CHUNK_IDX = 0;
			SCAN_NAME = '';
			SCAN_NREPS = [];
			ON_DK = [];
			ON_EL = [];
			EL_AT_START = [];
			IN_SCSET = '';
			if (IN_CYCLE == 1)
				F(end).t(2) = ON_TIME;
				IN_CYCLE = 0;
			end;
			MARKBITS = clear_markbits();
		end;
		[s,is_cryo] = check_for_phase_start (ll);
		if ~isempty (s)
			% disp (['In phase ' s]);
                        if isempty(ON_TIME)
                                ON_TIME = datenum (ON_TIME_STR, 'yymmdd HH:MM:SS');
                        end;
			% If we see a phase start in a 7_ schedule, it's a CMB-like
			% observation, not a mast raster or whatever.
			if IN_CAL
				IN_CAL = 0;
				SCSET_IDX = 0;
			end
			if ~isempty(IN_SCHED) && ~isempty(IN_PHASE) && ~isempty(IN_SCSET)
				if strcmp(IN_SCSET,'perphasecals')
					CHUNK_IDX = -1 * PHCAL_IDX;
					PHCAL_IDX = PHCAL_IDX + 1;
				else
					CHUNK_IDX = SCSET_IDX;
				end;
                                C = add_chunk(C, CURR_START_TIME, ON_TIME, IN_SCHED, SCHED_START_DATE, IN_SCSET, IN_PHASE, CHUNK_IDX, SCAN_NAME, SCAN_NREPS, ON_DK, EL_AT_START, MARKBITS);
				% disp ([datestr(CURR_START_TIME) ', ' datestr(ON_TIME) ', ' IN_SCHED ', ' SCHED_START_DATE IN_PHASE num2str(CHUNK_IDX)]);
			end;
			% Make sure not to get bad phase name listing for 2015-04-21 / 2015-04-22 BICEP3 schedules
			if isempty(findstr(s,'phaseletter'))
				IN_PHASE = s;
			end
			IS_CRYO_PHASE = is_cryo;
			if length(LAST_PHCAL_DATE_PH)<9 || ~strcmp(LAST_PHCAL_DATE_PH(1:8),SCHED_START_DATE) || ~strcmp(LAST_PHCAL_DATE_PH(9), s)
				PHCAL_IDX = 1;
			end;
			LAST_PHCAL_DATE_PH = [SCHED_START_DATE IN_PHASE];
			SCSET_IDX = 0;
			CHUNK_IDX = 0;
                        SCAN_NAME = '';
                        SCAN_NREPS = [];
			IN_SCSET = 'perphasecals';
			CURR_START_TIME = ON_TIME;
			if (IN_CYCLE == 1)
				F(end).t(2) = ON_TIME;
				IN_CYCLE = 0;
			end;
			MARKBITS = clear_markbits();
			EL_AT_START = ON_EL;
		end;
		s = check_for_scanset_start (ll, IN_CAL);
		if ~isempty(s) && strcmp(s,'perphasecals') && strcmp(IN_SCSET,'perphasecals')
			s = '';
		end;
		if ~isempty(s)
			% disp (['In scanset ' s]);
                        if isempty(ON_TIME)
                                ON_TIME = datenum (ON_TIME_STR, 'yymmdd HH:MM:SS');
                        end;
			if ~isempty (IN_SCHED) && ~isempty(IN_PHASE) && ~isempty(IN_SCSET)
				if strcmp(IN_SCSET,'perphasecals')
					CHUNK_IDX = -1 * PHCAL_IDX;
					PHCAL_IDX = PHCAL_IDX + 1;
				else
					CHUNK_IDX = SCSET_IDX;
				end;
				C = add_chunk(C, CURR_START_TIME, ON_TIME, IN_SCHED, SCHED_START_DATE, IN_SCSET, IN_PHASE, CHUNK_IDX, SCAN_NAME, SCAN_NREPS, ON_DK, EL_AT_START, MARKBITS);
				% disp ([datestr(CURR_START_TIME) ', ' datestr(ON_TIME) ', ' IN_SCHED ', ' SCHED_START_DATE IN_PHASE num2str(CHUNK_IDX)]);
			end;

			CURR_START_TIME = ON_TIME;
			if strcmp(s,'SKIP')
				IN_SCSET = '';
			else
				IN_SCSET = s;
			end
			if (IN_CAL == 1)
				IN_PHASE = '';
				IS_CRYO_PHASE = false;
			elseif ~strcmp (s, 'perphasecals')
				SCSET_IDX = SCSET_IDX+1;
			end;
			if (IN_CYCLE == 1)
				F(end).t(2) = ON_TIME;
				IN_CYCLE = 0;
			end;
			MARKBITS = clear_markbits();
			EL_AT_START = ON_EL;
		end;
		[tmp_sc tmp_scn] = check_for_scan (ll, IN_CAL);
		if ~isempty(tmp_sc)
			SCAN_NAME = tmp_sc;
			if IN_CAL & isempty(IN_PHASE)
				IN_PHASE = 'X';
				IS_CRYO_PHASE = false;
				CALDY_IDX = CALDY_IDX+1;
				SCSET_IDX = CALDY_IDX;
				% CALDY_IDX = SCSET_IDX;
				% SCSET_IDX = SCSET_IDX+1;
			end;
		end;
		if ~isempty(tmp_scn)
			SCAN_NREPS = tmp_scn;
		end;

		tmp = check_for_cycle_boundary (ll);
		if (tmp ~= 0)
			% disp (['Entering / leaving cycle, tmp=' num2str(tmp)]);
			if isempty(ON_TIME)
				ON_TIME = datenum (ON_TIME_STR, 'yymmdd HH:MM:SS');
			end;
                        if ~isempty(IN_SCHED) && ~isempty(IN_PHASE) && ~isempty(IN_SCSET) && ((IN_PHASE ~= 'A') && ~IS_CRYO_PHASE)
				% Handle weird cases where Robert has started a scanset while one RX is still
				% cycling.  Don't end the scanset just because last RX finished cycle.
				if (tmp < 0)
					continue;
				end
                                if strcmp(IN_SCSET,'perphasecals')
                                        CHUNK_IDX = -1 * PHCAL_IDX;
                                        PHCAL_IDX = PHCAL_IDX + 1;
                                else
                                        CHUNK_IDX = SCSET_IDX;
                                end;
                                C = add_chunk(C, CURR_START_TIME, ON_TIME, IN_SCHED, SCHED_START_DATE, IN_SCSET, IN_PHASE, CHUNK_IDX, SCAN_NAME, SCAN_NREPS, ON_DK, EL_AT_START, MARKBITS);
                                % disp ([datestr(CURR_START_TIME) ', ' datestr(ON_TIME) ', ' IN_SCHED ', ' SCHED_START_DATE IN_PHASE num2str(CHUNK_IDX)]);
				IN_PHASE = '';
                        end;
			IN_SCSET = '';
			SCSET_IDX = 0;
			CHUNK_IDX = 0;
			SCAN_NAME = '';
			SCAN_NREPS = [];
			if (tmp > 0)
				if (IN_CYCLE == 0)
					F(end+1).t = [ON_TIME NaN];
					% Fill in phase if available; otherwise default to 'A'
					if ~isempty(IN_PHASE) && IS_CRYO_PHASE
						F(end).ph = IN_PHASE;
					else
						F(end).ph = 'A';
					end
					IN_CYCLE = 1;
				end;
			else
				if length(F)>0 && (isnan(F(end).t(2)) || (1 > F(end).t(2)-F(end).t(1)))
					F(end).t(2) = ON_TIME;
				end;
				IN_CYCLE = 0;
			end;
			MARKBITS = clear_markbits();
			EL_AT_START = [];
		end;

		tmp = check_for_dk_offset (ll);
		if ~isempty (tmp)
			% disp (['Changed deck offset, dk=' num2str(tmp)]);
			ON_DK = tmp;
		end;

		tmp = check_for_el_offset (ll, ON_EL);
		if ~isempty (tmp)
			ON_EL = tmp;
		end

		[fplus,fminus] = check_for_mark(ll);
		if ~isempty(fplus) || ~isempty(fminus)
			if isempty(ON_TIME)
				ON_TIME = datenum (ON_TIME_STR, 'yymmdd HH:MM:SS');
			end

			for b=fplus
				if (b+1)>length(MARKBITS)
					continue
				end
				if MARKBITS(b+1).val==0
					MARKBITS(b+1).val=1;
					MARKBITS(b+1).ton=[MARKBITS(b+1).ton ON_TIME];
					% disp (['Turned on feature bit ' num2str(b)]);
				end
			end

			for b=fminus
				if (b+1)>length(MARKBITS)
					continue
				end
				if MARKBITS(b+1).val==1
					MARKBITS(b+1).val=0;
					MARKBITS(b+1).toff=[MARKBITS(b+1).toff ON_TIME];
					% disp (['Turned off feature bit ' num2str(b)]);
				end
			end
		end



	end;
	fclose (f);
end;

% Get rid of events that are outside requested time range
dn1 = utc2datenum(t1);
dn2 = utc2datenum(t2);

cut = true(length(C),1);
for ii=1:length(C)
	Cdn1 = datenum(C(ii).t1,'dd-mmm-yyyy:HH:MM:SS');
	Cdn2 = datenum(C(ii).t2,'dd-mmm-yyyy:HH:MM:SS');
	if (Cdn2<dn1) || (Cdn1>dn2)
		cut(ii) = 0;
	end;
end;
C = C(cut);
cut = true(length(F),1);
for ii=1:length(F)
	if (F(ii).t(1)>dn2)
		cut(ii) = 0;
	end;
	if length(F(ii).t)>1 && (F(ii).t(2)<dn1)
		cut(ii) = 0;
	end;
end;
F = F(cut);
cut = true(length(SCH),1);
for ii=1:length(SCH)
        if (SCH(ii).t(1)>dn2)
		cut(ii) = 0;
	end; 
        if length(SCH(ii).t)>1 && (SCH(ii).t(2)<dn1)
                cut(ii) = 0;
        end;
end;
SCH = SCH(cut);

return;

function s = check_for_sched_start (ll)
	s = '';
	sss = 'Starting schedule: ';
	if strncmp (ll, sss, length(sss))
		s = ll((1+length(sss)):end);
		return;
	end;
	sss = 'abort_schedule';
	if strncmp (ll, sss, length(sss))
		s = sss;
		return;
	end;
	return;

function [s,is_cryo] = check_for_phase_start (ll)
        s = '';
	is_cryo = false;
        sss = 'Now starting Phase ';
        if strncmp (lower(ll), lower(sss), length(sss))
                s = ll((1+length(sss)):end);
                s = strtok(s);
		is_cryo = ~isempty(strfind(ll,'cryo'));
		return;
        end;
	sss = 'Now starting 360 degree';
	if strncmp (ll, sss, length(sss))
		s = 'X';
		return;
	end;
        % Deal with BICEP3 early CMB schedules with broken log messages
        C = textscan (ll, 'do_obs_phase "%[^"]%*["]');
        if (length(C) > 0) && ~isempty(C{1})
                s = C{1}{1};
		return;
        end;
        C = textscan (ll, 'do_cryo_phase "%[^"]%*["]');
        if (length(C) > 0) && ~isempty(C{1})
                s = C{1}{1};
		is_cryo = true;
                return;
        end;   

	return;

function s = check_for_scanset_start (ll, iscal)
	if (nargin < 2) || isempty (iscal)
		iscal = 0;
	end;
        s = '';
	if (iscal == 1)
		sss = 'track ';
        	if strncmp (ll, sss, length(sss))
                	s = ll((1+length(sss)):end);
                	return;
        	end;
        	% In cals, "scanset" identification is
        	% harder.  Scans at different dk angles
        	% should certainly be separated.
		if strncmp (ll, 'offset', length('offset')) && ~isempty(strfind(ll, 'dk='))
			s = 'SKIP';
			return;
		end;
	else
        	sss = 'Now starting scanset ';
                if strncmp (ll, sss, length(sss))
                        s = ll((1+length(sss)):end);
                        return;
                end;

		sss = 'Too late for scanset ';
		if strncmp (ll, sss, length(sss))
 			s = 'SKIP';
			return;
		end;

		% Check for per-phase cals
		sss = 'do_skydip';
                if strncmp (ll, sss, length(sss))
			s = 'perphasecals';
			return;
		end;
		sss = 'do_loadcurve';
		if strncmp (ll, sss, length(sss))
			s = 'perphasecals';
			return;
		end;
	end;

	return;

function [s n] = check_for_scan (ll, iscal)
	if (nargin < 2) || isempty (iscal)
		iscal = 0;
	end;
	s = '';
	n = [];
	if (iscal == 1)
                C = textscan (ll, 'do_scans_%*s %d, %[^,]');
	else
		C = textscan (ll, 'do_scans 1, %d, %[^,]');
	end;
	if (length(C) == 2) && ~isempty(C{1}) && ~isempty(C{2})
		s = C{2}{1};
		n = C{1}(1);
	end;
	return;

function s = check_for_cycle_boundary (ll)
	s = 0;
	sss = 'cycleFridge';
	if strncmp (ll, sss, length(sss))
		tmp = ll((1+length(sss)):end);
		tmp = strtrim (tmp);
		if isempty(tmp)
			s = 1;
		elseif strncmp (tmp, 'state=0', 7) || strncmp (tmp, 'state=-1', 8)
			s = -1;
		end;
	end;
	return;

function s = check_for_dk_offset (ll)
	s = [];
	sss = 'offset';
	if strncmp (ll, sss, length(sss))
		tmp = regexp (ll, 'dk=([-+.:0-9]*[0-9])', 'tokens');
                if length(tmp)>=1 && length(tmp{1})>=1
			s = parse_angle (tmp{1}{1});
			s = mod(s, 360);
		end
	end;
	return;

function s = check_for_el_offset (ll, elofs)
	s = [];
	sss = 'offset';
	if strncmp (ll, sss, length(sss))
		sss = 'offset/add';
		if isempty (elofs) || ~strncmp (ll, sss, length(sss))
			elofs = 0;
		end
                tmp = regexp (ll, 'el=([-+.:0-9]*[0-9])', 'tokens');
                if length(tmp)>=1 && length(tmp{1})>=1
                        s = elofs + parse_angle (tmp{1}{1});
                end
	end;
	return;

function fbits = clear_markbits()
	for k=1:30
        	fbits(k).val=0;
        	fbits(k).ton=[];
        	fbits(k).toff=[];
	end
	return;

function [fplus,fminus]=check_for_mark(ll)
	fplus=[];
	fminus=[];
	time=NaN;
	if ~strncmp(ll, 'mark', 4)
		return
	end
	sss = 'mark add';
	if strncmp (ll, sss, length(sss))
		[tmp,tpl]=strtok(ll,',');
		tpl=strtrim(tpl);
		if tpl(1)==','
			tpl(1)='';
		end
		tpl=strtok(tpl,',');
		tpl(tpl=='f' | tpl=='+') = ' ';
		fplus=str2num(tpl);
		% fplus(fplus==2) = [];
		return
	end
	sss = 'mark remove';
	if strncmp (ll, sss, length(sss))
		[tmp,tpl]=strtok(ll,',');
		tpl=strtrim(tpl);
		if tpl(1)==','
			tpl(1)='';
		end
		tpl=strtrim(tpl);
		if strcmp(tpl,'all')
			fminus=[0:30];
			% fminus(fminus==2) = [];
			return
		end
		tpl=strtok(tpl,',');
		tpl(tpl=='f' | tpl=='+') = ' ';
		fminus=str2num(tpl);
		% fminus(fminus==2) = [];
		return
	end
	% s=check_for_cycle_boundary(ll);
	% if s==1
	% 	fplus=2;
	% 	return
	% elseif s==-1
	% 	fminus=2;
	% 	return
	% end
	return;

function C = add_chunk (C, t1, t2, sch, t0, set, ph, ch, sc, nsc, dk, elofs, fbits)
	[tmppath tmpfname tmpext] = fileparts (sch);
	sch = [tmpfname tmpext];
	tmp = regexp (sch, '.*dk(\d{1,3})_','tokens');
	if ~isempty(tmp)
		dk = str2num(tmp{1}{1});
	end;
	% print_chunk_line (t1, t2, sch, t0, ph, ch, sc, nsc, dk);

	Ctmp = [];

	Ctmp.t1 = datestr(t1, 'dd-mmm-yyyy:HH:MM:SS');
	Ctmp.t2 = datestr(t2, 'dd-mmm-yyyy:HH:MM:SS');
	Ctmp.sch = sch;
	Ctmp.t0 = t0;
        Ctmp.scset = set;
	Ctmp.ph = ph;
	Ctmp.ch = ch;
	Ctmp.sc = sc;
	Ctmp.nsc = nsc;
	Ctmp.dk = dk;
	Ctmp.elofs = elofs;

	% Number of starts of various operations
	Ctmp.nmark(1,1) = 1;						% 1. Scansets (=1)
	Ctmp.nmark(1,2) = length(fbits(1+1).ton);			% 2. On source periods
	Ctmp.nmark(1,3) = length(fbits(3+1).ton);			% 3. El nods
	Ctmp.nmark(1,4) = length(fbits(14+1).ton);			% 4. Partial load curves
	Ctmp.nmark(1,5) = length(fbits(2+1).ton);			% 5. Fridge cycles
	Ctmp.nmark(1,6) = length(fbits(7+1).ton);			% 6. Sky dips
	Ctmp.nmark(1,7) = length(fbits(6+1).ton);			% 7. Full load curves

	% Number of stops of various operations
        Ctmp.nmark(2,1) = 1;                                            % 1. Scansets (=1)
        Ctmp.nmark(2,2) = length(fbits(1+1).toff);                      % 2. On source periods
        Ctmp.nmark(2,3) = length(fbits(3+1).toff);                      % 3. El nods
        Ctmp.nmark(2,4) = length(fbits(14+1).toff);                     % 4. Partial load curves
        Ctmp.nmark(2,5) = length(fbits(2+1).toff);                      % 5. Fridge cycles
        Ctmp.nmark(2,6) = length(fbits(7+1).toff);                      % 6. Sky dips
        Ctmp.nmark(2,7) = length(fbits(6+1).toff);                      % 7. Full load curves

        % times spent in various operations:
	% first fill in missing end times
	for k=1:length(fbits)
		if length(fbits(k).toff)<length(fbits(k).ton)
			fbits(k).toff(end+1:length(fbits(k).ton))=t2;
		end
	end
        Ctmp.livetime(1) = t2-t1;					% 1. Total time
	Ctmp.livetime(2) = sum(fbits(1+1).toff-fbits(1+1).ton);		% 2. On source
	Ctmp.livetime(3) = sum(fbits(3+1).toff-fbits(3+1).ton);		% 3. El nod
        Ctmp.livetime(4) = sum(fbits(14+1).toff-fbits(14+1).ton);	% 4. Partial load curve
        Ctmp.livetime(5) = sum(fbits(2+1).toff-fbits(2+1).ton);		% 5. Fridge cycle
        Ctmp.livetime(6) = sum(fbits(7+1).toff-fbits(7+1).ton);		% 6. Sky dip
        Ctmp.livetime(7) = sum(fbits(6+1).toff-fbits(6+1).ton);		% 7. Full load curve
	Ctmp.livetime = Ctmp.livetime * 24 * 60;			% convert to minutes

	% Don't add "per-phase cals" tags that don't actually contain any action
	if all(Ctmp.livetime(2:end)==0) && all(Ctmp.nmark(1,2:end)==0) && all(Ctmp.nmark(2,2:end)==0)
		return
	end

	if isempty(C)
		C = Ctmp;
	else
		C(end+1) = Ctmp;
	end

	return;

function print_chunk_line (t1, t2, sch, t0, ph, ch, sc, nsc, dk)
	disp (['ae(r,''' datestr(t1,'dd-mmm-yyyy:HH:MM:SS') ''', ''' ...
                       datestr(t2,'dd-mmm-yyyy:HH:MM:SS') ''', ''' ...
                       sch ''', ''' t0 lower(ph) num2str(ch) '''); % ' ...
                       t0 ' Phase ' ph ' part ' num2str(ch) ...
                       ' scan=' sc ', nreps=' num2str(nsc) ...
		       ' dk=' num2str(dk)]);
	return;

function print_all_chunks (f, C)
	for (ii = 1:length(C))
		fprintf (f, ['r=ae(r,''' C(ii).t1 ''', ''' C(ii).t2 ''', ''' ...
                       C(ii).sch ''', ''' C(ii).t0 lower(C(ii).ph) num2str(C(ii).ch) '''); %% ' ...
                       C(ii).t0 ' Phase ' C(ii).ph ' part ' num2str(C(ii).ch) '\n']);
	end;
	return;

