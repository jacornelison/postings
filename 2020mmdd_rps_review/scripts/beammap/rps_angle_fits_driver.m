function rps_angle_fits_driver(sch,p,rpsopt,dirname,fd,rn)


% Figure out whether to run all rows or only a small subset.
if ~exist('rn', 'var')
	rn = 1:length(p.gcp);
end

if ~isfield(rpsopt,'fitopt')
	rpsopt.fitopt = [];
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
for ch=rn
    if farming > 0
        pid=myprocesses('farmit\(');
        while numel(pid)>=farming
            pid=myprocesses('farmit\(');
            pause(1);
        end
        [farmfile,jobname] = farmfilejobname('BICEP', 'rps_angfit', ...
            sprintf('CH_%04i_%d%d',ch,floor(now),floor(rem(now,1)*1e7)));
        
        farmit(farmfile, ...
            'rps_fit_angle(sch, p, rpsopt, dirname, fd, ch)', ...
            'jobname',jobname, 'queue',QUEUE, 'account','bicepdata_group', ...
            'mem',MEM, 'maxtime',MAXTIME,'overwrite',true, ...
            'var', {'sch', 'p', 'rpsopt', 'dirname','fd', 'ch'});
    else
        rps_fit_angle(sch, p, rpsopt, dirname, fd, ch)
    end
end

