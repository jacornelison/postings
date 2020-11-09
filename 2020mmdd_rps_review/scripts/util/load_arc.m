% LOAD_ARC  reads arcfile data as written by GCP.
%
% D = LOAD_ARC(F) reads all registers from the file F,
% creating a structure D.
%
% D = LOAD_ARC(F,T1,T2) selects only data taken between
% times T1 and T2.  These are UTC times specified in
% one of the following formats:
% 2009-Apr-29:14:32:01 (Walt format)
% 29-Apr-2009:14:32:01 (gcp "show" command format)
% 110116 14:20:01      (log file format)
%
% D = LOAD_ARC(F,T1,T2,C), where C is a string or cell
% array of strings, reads only the fields specified
% in C.  Wildcards are allowed.  Channel numbers are
% also allowed, e.g. mce0.data.fb(5:10).  Note that
% channel numbers are indexed from 1.  (Should this
% change to zero?)
%
% If F is a directory instead of a single file, all arc
% files within the specified time range will be included.
%
% F can also be a remote host like bicep0, for remote
% data access:
% D = load_arc ('reuben@bicep0:/data/bicep2daq/arc/', ...)
%
% D = LOAD_ARC(F,T1,T2,C,'DO_UNPACK',1) handles older data
% in which only the raw MCE words were archived.  With
% this option, If the registers mce0.data.fb, .err, or .numfj
% are requested, they will automatically be unpacked
% from mce0.data.raw.
%
% D = LOAD_ARC(F,T1,T2,C,'DO_MERGEFIX',1) tries to repair
% failures in merging antenna and MCE data.  This only works
% when there's no GCP filtering or downsampling, i.e. with
% combine off.  Not for use with CMB data.
%
% RWO 090505

function d = load_arc (f, t1, t2, c, varargin)

USE_TIME = 'array.frame.utc';

% Check for any parameter options
if (nargin < 5) || isempty (varargin)
	varargin = {};
end;
[opts.DO_UNPACK opts.DO_MERGEFIX] = paramparse ({'DO_UNPACK','DO_MERGEFIX'}, {0,0}, varargin{:});

% Check for host name
[f1 f2] = strtok (f, ':');
if ~isempty (f2)
	srv = f1;
	f = f2(2:end);
else
	srv = '';
	if (~exist (f, 'dir') && ~exist (f, 'file'))
		error (['Must specify a file name or directory.']);
	end;
end;

if (nargin >= 4) && ~isempty (c)
	if ischar (c)
		c = {c};
	end;
	c = {c{:}, USE_TIME};
	c = parse_fields_chans (c);
else
	c = {};
end;

if (nargin >= 3) && (~isempty(t1) && ~isempty(t2))
	t1 = txt2txt (t1);
	t2 = txt2txt (t2);
	utc1 = txt2utc(t1);
	utc2 = txt2utc(t2);
else
	t1 = '';
	t2 = '';
        utc1 = [];
        utc2 = [];
end;

d = [];

if (opts.DO_UNPACK == 1)
	% If user asked for err signal, fb, etc., figure out
	% what information we need to calculate these.
	% Add the necessary raw registers to the list.
	[c L] = add_mce_stuff(c);
end;

% disp (['Asking for registers:']);
% for (ii = 1:length(c))
% 	disp (['    ' c{ii}]);
% end;


if isempty (srv)
	d = readarc (f, t1, t2, c);
else
	d = readarc_remote (srv, f, t1, t2, c{:});
end;
if isempty (d)
	return;
end;

% Fix merging failures if requested (and if possible)
if (opts.DO_MERGEFIX == 1)
	d = fix_merging (d);
end;

% Cut on requested time ranges
if (nargin >= 3) && (~isempty (utc1) && ~isempty (utc2))
	d = data_filter (d, utc1, utc2, USE_TIME);
end;

if (opts.DO_UNPACK == 1)
	% If user asked for fb, err signal, etc., unpack these
	% from the raw MCE words
	d = unpack_mce_stuff(d,L);
end;

% Unpack UTC times from 64-bit int to 2-column double [MJD, tod in s]
d = unpack_utc_times(d);

%%%%

function uu = txt2txt (tt)
	FMTLIST = {'yyyy-mmm-dd:HH:MM:SS', 'yyyy-mmm-dd:HH:MM', ...
		   'yyyy-mmm-dd:HH', 'yyyy-mmm-dd', ...
		   'dd-mmm-yyyy:HH:MM:SS', 'dd-mmm-yyyy:HH:MM', ...
		   'dd-mmm-yyyy:HH', 'dd-mmm-yyyy', ...
		   'yymmdd HH:MM:SS'};
	tmp = [];
	for (ii = 1:length (FMTLIST))
		try
			tmp = datenum (tt, FMTLIST{ii});
			% NOTE: As of R2013b datenum throws errors
			% much less often.  It often returns valid
			% dates near the present even with the wrong
			% format.  This bit of code may need
			% rethinking to make it more robust.
		catch
			tmp = 0;
		end;
		if (tmp >= 732313 && tmp <= now) % years 2005-date
			break;
		end;
	end;
	uu = datestr(tmp,'yyyy-mmm-dd:HH:MM:SS');

function uu = txt2utc (tt)
        tmp = datenum (tt,'yyyy-mmm-dd:HH:MM:SS');
	uu = tmp + 48987 - datenum (1993,1,1) + 1;

function fl2 = file_filter (fl1, t1, t2)
	for (ii = 1:length (fl1))
		[a b] = fileparts (fl1{ii});
		try
			uu(ii) = datenum (b, 'yyyymmdd_HHMMSS');
		catch
			disp (['Bad file name ' b '.']);
			uu(ii) = NaN;
		end;
		uu(ii) = uu(ii) + 48987 - datenum ('Jan-1-1993') + 1;
	end;
	[uu ii] = sort (uu);
	fl2 = {fl1{ii}};
	cc = (uu >= t1) & (uu <= t2);

	% identify last file within 1 day of t1;
	% it may start before t1, but extend past t1
	npre = find ((uu < t1) & (t1 - uu < 1));
	cpre = zeros (size (cc));
	if ~isempty (npre)
		cpre(npre(end)) = 1;
	end;
	fl2 = {fl2{(cc|cpre)}};

function d2 = data_filter (d1, t1, t2, USE_TIME)
	if (nargin < 4) || isempty (USE_TIME)
		USE_TIME = 'array.frame.utc';
	end;
	tmp = double (typecast (eval (['d1.' USE_TIME]), 'uint32'));
	tt = (tmp(1:2:end)) + (tmp(2:2:end)) / (24*3600*1000);
	cc = (tt >= t1) & (tt <= t2);
	d2 = structcut (d1, cc);

function d2 = structcut (d1, cc)
	if ~isstruct (d1)
		if (size(d1,1) == length(cc))
			d2 = d1 (cc, :);
		else
			spf = size(d1,1) / length(cc);
			% disp (['spf = ' num2str(spf)]);
			d2 = reshape (d1, spf, length(cc), size(d1,2));
			d2 = d2 (:, cc, :);
			d2 = reshape (d2, sum(cc)*spf, size(d1,2));
		end;
	else
		fl = fieldnames (d1);
		for (ii = 1:length(fl))
			d2.(fl{ii}) = structcut (d1.(fl{ii}), cc);
		end;
	end;

function ch_out = str_chan_list (ch)
	ch_out = '';
	ch = sort (ch);
	idx1 = [1; (diff(ch(:))~=1)];
	idx2 = [(diff(ch(:))~=1); 1];
	idx1 = ch(logical(idx1));
	idx2 = ch(logical(idx2));
	for (ii = 1:length(idx1))
		if ~isempty(ch_out)
			ch_out = [ch_out ','];
		end;
		if (idx1(ii) == idx2(ii))
			ch_out = [ch_out num2str(idx1(ii))];
		else
			ch_out = [ch_out num2str(idx1(ii)) ':' num2str(idx2(ii))];
		end;
	end;
	
function c_out = parse_fields_chans (c)
	c_out = c;
	for (ii = 1:length(c))
		[c_out{ii} tmp] = strtok (c{ii}, '(');
                if isempty(tmp)
                	continue
                end

                if length(tmp)>0 && tmp(1)=='('
                	tmp(1)='';
                end
		[tmp tmp2] = strtok ([' ' tmp], ')');
		tmp = strtrim(tmp);
                if ~isempty(tmp)
                        tmp = str2num(tmp);
			tmp = tmp - 1;
			chans{ii} = str_chan_list(tmp);
			chans{ii} = ['[' chans{ii} ']'];
		else
			chans{ii} = '[]';
		end;
		c_out{ii} = [c_out{ii} chans{ii}];
	end;

% If user asks for mce0.data.fb, mce0.data.err, or mce0.data.numfj,
% we can get mce0.data.raw from the file and do unpacking for them.
% Also be sure to get mce0.rc*.data_mode, mce0.rc*.fw_rev.
% Have to keep track of which registers/channels the user has
% actually asked for, and which ones we just need internally.
function [c_out L] = add_mce_stuff (c)
	ALL_MCE_LIST = 0;

	L.out_raw_list = false (1, 528);
	L.out_fb_list  = false (1, 528);
	L.out_err_list = false (1, 528);
	L.out_nfj_list = false (1, 528);
	L.tmp_rc_dm    = false (1, 2);
	L.tmp_raw_list = false (1, 528);

	c_out = {};

	if isempty(c)
		c = {'*'};
	end;

	for (ii = 1:length (c))
		[rmap tmp] = strtok (c{ii}, '.');
		if strcmp (rmap, 'mce*') || strcmp (rmap, '*')
			mce_num = ALL_MCE_LIST;
                % elseif ~regexp (rmap, '^mce[0-9]+$')
		elseif isempty (regexp (rmap, '^mce[0-9]+$'))
			c_out = {c_out{:}, c{ii}};
			continue;
		else
			mce_num = str2num (rmap (4:end));
		end;
		if ~isempty(tmp) && tmp(1)=='.'
			tmp = tmp(2:end);
		end;
		if isempty(tmp)
			c_out = {c_out{:}, c{ii}};
			tmp = 'data';
			% continue;
		end;
		[brd reg] = strtok (tmp, '.');
		if ~isempty(reg) && reg(1)=='.'
			reg = reg(2:end);
		end;
		if ~strcmp (brd, 'data')
			c_out = {c_out{:}, c{ii}};
			continue;
		end;
		if isempty(reg)
			L.tmp_raw_list(mce_num+1, :) = 1;
			L.tmp_rc_dm(mce_num+1, :) = 1;
			L.out_raw_list(mce_num+1, :) = 1;
			L.out_fb_list(mce_num+1, :) = 1;
			L.out_err_list(mce_num+1, :) = 1;
			L.out_nfj_list(mce_num+1, :) = 1;
			continue;
		end;
		[reg ch] = strtok (reg, '[');
		reg (reg == ' ') = '';
		if ~isempty (ch)
			ch(ch=='-') = ':';
			ch_list = 1 + eval (ch);
		else
			ch_list = 1:528;
		end;
		if strcmp (reg, 'raw')
			L.out_raw_list(mce_num+1,ch_list) = 1;
			L.tmp_raw_list(mce_num+1,ch_list) = 1;
		elseif strcmp (reg, 'fb')
			L.out_fb_list(mce_num+1,ch_list) = 1;
			L.tmp_raw_list(mce_num+1,ch_list) = 1;
			L.tmp_rc_dm(mce_num+1,:) = 1;
		elseif strcmp (reg, 'err')
			L.out_err_list(mce_num+1,ch_list) = 1;
			L.tmp_raw_list(mce_num+1,ch_list) = 1;
			L.tmp_rc_dm(mce_num+1,:) = 1;
		elseif strcmp (reg, 'numfj')
			L.out_nfj_list(mce_num+1,ch_list) = 1;
			L.tmp_raw_list(mce_num+1,ch_list) = 1;
			L.tmp_rc_dm(mce_num+1,:) = 1;
		elseif strcmp (reg, '*')
                        L.out_raw_list(mce_num+1,ch_list) = 1;
                        L.out_fb_list(mce_num+1,ch_list) = 1;
                        L.out_err_list(mce_num+1,ch_list) = 1;
			L.out_nfj_list(mce_num+1,ch_list) = 1;
                        L.tmp_raw_list(mce_num+1,ch_list) = 1;
                        L.tmp_rc_dm(mce_num+1,:) = 1;
		else
			c_out = {c_out{:}, c{ii}};
		end;
	end;
	for (ii = 1:size (L.tmp_raw_list, 1))
		if sum (L.tmp_raw_list(ii,:)) > 0
			c_out = {c_out{:}, ['mce' num2str(ii-1) '.data.raw[' str_chan_list(-1+find(L.tmp_raw_list(ii,:))) ']']};
		end;
		if sum (L.tmp_rc_dm(ii,:)) > 0
			c_out = {c_out{:}, ['mce' num2str(ii-1) '.rc*.data_mode'], ...
					   ['mce' num2str(ii-1) '.rc*.fw_rev']};
		end;
	end;

% Gain of MCE butterworth filter depends on firmware revision number
function g = mce_filter_gain (rc_fw_rev)
	if (rc_fw_rev == hex2dec('05000007'))
		g = 2044;
	else
		g = 1218;
	end;

% Now perform the unpacking, and strip out all the
% unneeded stuff
function d_out = unpack_mce_stuff (d_in, L)
	d_out = d_in;
	clear d_in;
	c_rc1 = mod ((1:528) - 1, 16) < 8;
	c_rc2 = ~c_rc1;
	for (ii = 1:size (L.tmp_raw_list,1))
		c_chan = logical (L.tmp_raw_list(ii,:));
		if sum (c_chan) == 0
			continue;
		end;
		disp (['Unpacking mce' num2str(ii-1)]);
		do_r = sum (L.out_raw_list(ii,:)) > 0;
		do_f = sum (L.out_fb_list(ii,:)) > 0;
		do_e = sum (L.out_err_list(ii,:)) > 0;
		do_n = sum (L.out_nfj_list(ii,:)) > 0;
		if ~do_f & ~do_e & ~do_n
			continue;
		end;

		rr = d_out.(['mce' num2str(ii-1)]).data.raw;
		d_out.(['mce' num2str(ii-1)]).data = rmfield (d_out.(['mce' num2str(ii-1)]).data, 'raw');
		dm_rc1 = d_out.(['mce' num2str(ii-1)]).rc1.data_mode;
		dm_rc2 = d_out.(['mce' num2str(ii-1)]).rc2.data_mode;
		fw_rc1 = d_out.(['mce' num2str(ii-1)]).rc1.fw_rev;
		fw_rc2 = d_out.(['mce' num2str(ii-1)]).rc2.fw_rev;

		if (do_f)
			fb = zeros (size (rr));
		end;
		if (do_e)
			err = zeros (size (rr));
		end;
		if (do_n)
			nfj = zeros (size (rr), 'uint16');
		end;

		NSAMP = size (rr, 1) / length (dm_rc1);

		c_chan = logical (L.tmp_raw_list(ii,:));
                dm_rc1_list = unique (dm_rc1);
		fw_rc1_list = unique (fw_rc1);
		for (jj = 1:length(dm_rc1_list))
			if sum (c_rc1(c_chan)) == 0
				continue;
			end;
			c_dm = reshape (repmat (dm_rc1' == dm_rc1_list(jj), NSAMP, 1), [], 1);
			for (kk = 1:length(fw_rc1_list))
				c_fw = reshape (repmat (fw_rc1' == fw_rc1_list(kk), NSAMP, 1), [], 1);
				if sum(c_fw&c_dm) == 0
					continue;
				end;
				[ftmp atmp] = mce_unpack (rr(c_dm&c_fw,c_rc1(c_chan)), dm_rc1_list(jj), fw_rc1_list(kk));
				if (do_f)
					fb (c_dm&c_fw,c_rc1(c_chan)) = ftmp;
				end;
				clear ftmp
				if ismember (dm_rc1_list(jj), [0 4 6 7])
					if (do_e)
						err (c_dm&c_fw,c_rc1(c_chan)) = atmp;
					end;
				else
					if (do_n)
						nfj (c_dm,c_rc1(c_chan)) = atmp;
					end;
				end;
				clear atmp
			end;
		end;
                dm_rc2_list = unique (dm_rc2);
                fw_rc2_list = unique (fw_rc2);
                for (jj = 1:length(dm_rc2_list))
			if sum (c_rc2(c_chan)) == 0
				continue;
			end;
                        c_dm = reshape (repmat (dm_rc2' == dm_rc2_list(jj), NSAMP, 1), [], 1);
			for (kk = 1:length(fw_rc2_list))
				c_fw = reshape (repmat (fw_rc2' == fw_rc2_list(kk), NSAMP, 1), [], 1);
				if sum(c_fw&c_dm) == 0
					continue;
				end;
                        	[ftmp atmp] = mce_unpack (rr(c_dm&c_fw,c_rc2(c_chan)), dm_rc2_list(jj));
				if (do_f)
                        		fb (c_dm&c_fw,c_rc2(c_chan)) = ftmp;
				end;
				clear ftmp
                        	if ismember (dm_rc2_list(jj), [0 4 6 7])
					if (do_e)
                                		err (c_dm&c_fw,c_rc2(c_chan)) = atmp;
					end;
                        	else
					if (do_n)
                                		nfj (c_dm&c_fw,c_rc2(c_chan)) = atmp;
					end;
                        	end;
				clear atmp
                	end;
		end;
		if (do_r)
			d_out.(['mce' num2str(ii-1)]).data.raw = rr(:, L.out_raw_list(ii,c_chan));
			clear rr
		end;
		if (do_f)
                	d_out.(['mce' num2str(ii-1)]).data.fb = fb(:, L.out_fb_list(ii,c_chan));
			clear fb
		end;
		if (do_e)
                	d_out.(['mce' num2str(ii-1)]).data.err = err(:, L.out_err_list(ii,c_chan));
			clear err
		end;
		if (do_n)
                	d_out.(['mce' num2str(ii-1)]).data.numfj = nfj(:, L.out_nfj_list(ii,c_chan));
			clear nfj
		end;
	end;

% Unpack any UTC times from 64-bit unsigned integer
% to 2-column [MJD, time of day in s]
function d_out = unpack_utc_times (d_in)

	d_out = d_in;
	maps = fieldnames (d_out);
	for (im = 1:length(maps))
		brds = fieldnames (d_out.(maps{im}));
		for (ib = 1:length(brds))
			regs = fieldnames (d_out.(maps{im}).(brds{ib}));
			for (ir = 1:length(regs))
				if strcmp (class (d_out.(maps{im}).(brds{ib}).(regs{ir})), 'uint64')
					tmp = d_out.(maps{im}).(brds{ib}).(regs{ir});
					if isempty(tmp)
						tmp = [zeros(size(tmp),'double'), zeros(size(tmp),'double')];
						continue;
					end;
					tmp = [utc2mjd(tmp), utc2time(tmp)];
					d_out.(maps{im}).(brds{ib}).(regs{ir}) = tmp;
				end;
			end;
		end;
	end;

% Try to resolve merging failures.  Only works with combine off.
% Bit of a hack.  Not for CMB data.
function d_out = fix_merging (d_in)

	d_out = d_in;
	regmaps = fieldnames (d_out);
	fsn_ant_unsorted = d_out.antenna0.syncBox.frame_serialNumber;
	fsn_ant = fsn_ant_unsorted;
	fsn_ant(fsn_ant~=0) = sort (fsn_ant(fsn_ant~=0));
	nframes = length (fsn_ant);
	if (nframes < 2)
		return
	end
	sampratio = length(d_out.antenna0.syncBox.sampleNumber) / nframes;
	for (i = 1:length (regmaps))
		if isfield (d_out.(regmaps{i}), 'syncBox')
			fsn = d_out.(regmaps{i}).syncBox.frame_serialNumber;
		elseif strncmp (regmaps{i},'mce',3)
			error(['Must load mce*.syncBox for DO_MERGEFIX to work.  Please include *.syncBox in register list.']);
		else
			fsn = fsn_ant_unsorted;
		end
		jlist = find (fsn>0 & fsn~=fsn_ant);
		disp([regmaps{i} ' : ' num2str(length(jlist)) ' mismatches.']);
		pslow = zeros (size (fsn));
		pslow(fsn==fsn_ant) = find(fsn==fsn_ant);
		for (j = 1:length(jlist))
			k = find (fsn_ant == fsn(jlist(j)));
			if (length (k) ~= 1 || jlist(j)==k)
				continue
			end;
			% Make sure no "combine" downsampling
			if (d_out.array.frame.nsnap > 1)
				continue
			end
			pslow(k) = jlist(j);
		end
		brds = fieldnames (d_out.(regmaps{i}));
		for (ib = 1:length(brds))
			regs = fieldnames (d_out.(regmaps{i}).(brds{ib}));
			for (ir = 1:length(regs))
				tmp = d_out.(regmaps{i}).(brds{ib}).(regs{ir});
				nperframe = size (tmp, 1) / nframes;
				if (nperframe == 0)
					continue
				end;
				for (l = 1:nperframe)
					tmptmp = tmp (l:nperframe:length(tmp), :);
					tmptmp2 = tmptmp;
					tmptmp2(:) = 0;
					tmptmp2 (pslow>0, :) = tmptmp (pslow(pslow>0), :);
					tmp (l:nperframe:length(tmp), :) = tmptmp2;
				end
				d_out.(regmaps{i}).(brds{ib}).(regs{ir}) = tmp;
               		end;
		end;
	end;

