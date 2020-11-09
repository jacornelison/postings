function rps_reduc_driver(sch, p, p_ind, rpsopt, dirname, rn,QUEUE)
% function rps_reduc_driver(sch, p, rpsopt, dirname, rn,QUEUE)
%
% reduces timestreams given by an RPS schedule file using
% sch = rps_log_read(t1,t2);
%
% [ Arguments ]
%	sch		
% 	rpsopt	Optional variable that can be filled with custom reduc options.
%	p		Focal Plane data returned by get_array_info.
%	rn		Single out a specific scan or specify a range of scans to reduce. 
%			default = all
%
% Run order + description
%	rps_find_channels: returns list of channels that are on source
%	rps_read: returns demodulated data
%	rps_fit_beam_2: Fits data to gaussian beam and returns parameters.
%
% [3 Aug 2016] Both find_channels and rps_read run separate load_arc functions.
% The time that it saves selecting channels is more than made up by
% reloading arc files. We need to find a way to mitigate this.
%
% Last Update: 1 Mar 2018 JAC
% function rps_reduc_driver(sch, p, rpsopt, dirname, rn)

if ~exist('p', 'var')
  p = [];
end

if ~exist('rpsopt', 'var')
  rpsopt = [];
end

if isempty(rpsopt)
  load('work/rps/2018/b3_20180104_rpsopt.mat')
end

% Figure out whether to run all scans or a small subset.
if exist('rn', 'var') & length(rn) ~= 0
	rn_hold = 1;
else
	rn_hold = 0;
end

if ~exist('dirname', 'var') | isempty(dirname)
  %Savefile directory
  dirname = 'rps_data/2019_mirror/';
end

if ~exist('QUEUE','var')
    QUEUE = 'kovac';
end

%implement farming of data
farming = 1000; % # of jobs to run at one time.

LICENSE = 'bicepfs1:17'; % 1500 slots -> 1500/20 = 88 jobs
MAXTIME = 300; % in minutes
MEM=5000;

%Loop over found schedules
for i=1:length(sch)
  %Loop over scans
  schname = sch{i}.name(end-15:end-4);
  if rn_hold==0
	rn = 1:length(sch{i}.scans);
  end
  for j=rn
    if farming > 0
		pid=myprocesses('farmit\(');
		while numel(pid)>=farming
			pid=myprocesses('farmit\(');
			pause(1);
		end
		[farmfile,jobname] = farmfilejobname('BICEP', 'rps_reduc', ...
			sprintf('SCH_%s_SCAN_%03i_%d%d',schname,j,floor(now),floor(rem(now,1)*1e7)));
		disp(jobname)
		farmit(farmfile, ...
			'rps_reduc(sch, i, j, p, p_ind, rpsopt, dirname)', ...
			'jobname',jobname, 'queue',QUEUE, 'account','bicepdata_group', ...
			'mem',MEM, 'maxtime',MAXTIME,'overwrite',true, ...
			'var', {'sch', 'i', 'j', 'p', 'p_ind', 'rpsopt', 'dirname'});
	else
		rps_reduc(sch,i,j,p,p_ind,rpsopt,dirname)
	end
  end
end
end


