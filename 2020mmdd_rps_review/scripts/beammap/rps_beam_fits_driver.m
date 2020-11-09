function rps_beam_fits_driver(sch,p,rpsopt,dirname,rn)
% function rps_beam_fits_driver(sch,p,rpsopt,dirname,rn)
% 
% Farms out rps_fit_beam for every elevation offset in a schedule. It is 
% assumed that all timestreams in the schedule have been reduced into rpstod.mat
% formats.
%
%
%
%
%
% function rps_beam_fits_driver(sch,p,rpsopt,dirname,rn)

% Figure out whether to run all rows or only a small subset.
if exist('rn', 'var')
	rn_hold = 1;
else
	rn_hold = 0;
end

%Savefile directory
if ~exist('dirname','var') & ~isempty(dirname)
dirname = 'rps_data/2019_mirror/';
fprintf(['No directory given. Defaulting to ' dirname '\t'])
end

%implement farming of data
farming = 1000; % # of jobs to run at one time.
QUEUE = 'kovac';
%LICENSE = 'bicepfs1:17'; % 1500 slots -> 1500/20 = 88 jobs
MAXTIME = 300; % in minutes
MEM=8000;

%Loop over found schedules
for i=1:length(sch)
  schname = sch{i}.name(end-15:end-4);
  
  if rn_hold==0
	rn = 1:sch{i}.nrows;
  end
  for j=rn
    if farming > 0
		pid=myprocesses('farmit\(');
		while numel(pid)>=farming
			pid=myprocesses('farmit\(');
			pause(1);
		end
		[farmfile,jobname] = farmfilejobname('BICEP', 'rps_beamfit', ...
			sprintf('SCH_%s_row_%02i_%d%d',schname,j,floor(now),floor(rem(now,1)*1e7)));
		
		farmit(farmfile, ...
			'rps_fit_beam(sch, i, j, p, rpsopt, dirname)', ...
			'jobname',jobname, 'queue',QUEUE, 'account','bicepdata_group', ...
			'mem',MEM, 'maxtime',MAXTIME,'overwrite',true, ...
			'var', {'sch', 'i', 'j', 'p', 'rpsopt', 'dirname'});
	else
		rps_fit_beam(sch, i, j, p, rpsopt, dirname)
	end
  end
end

