function rps_angle_allfits_driver(sch,p,rpsopt,dirname,fd,chans,nchans)
% function rps_angle_allfits_driver(sch,p,rpsopt,dirname,fd,chans)
% Fitting across all channels takes forever. Farm out subsets of selected
% channels. Emperically, 6 channels takes about an hour. So do that.


% Figure out whether to run all rows or only a small subset.
if ~exist('chans', 'var') | isempty(chans)
	chans = unique(fd.ch)';
end

if ~exist('nchans', 'var')
	nchans = 6;
end


% Try to avoid bias due to instrumental effects. Jumble selected channels.
if nchans>1
    chans = chans(randperm(length(chans)));
    MAXTIME = 12*60; % in minutes
else
    MAXTIME = 260; % in minutes
end

if ~isfield(rpsopt,'fitopt')
	rpsopt.fitopt = [];
end

%Savefile directory
if ~exist('dirname','var') & ~isempty(dirname)
dirname = 'rps_data/2019_mirror/';
fprintf(['No directory given. Defaulting to ' dirname '\t\n'])
end

%implement farming of data
farming = 1000; % # of jobs to run at one time.
QUEUE = 'kovac';
%LICENSE = 'bicepfs1:17'; % 1500 slots -> 1500/20 = 88 jobs
MAXTIME = 12*60; % in minutes
MEM=8000;

rn = ceil(length(chans)/nchans);
%Loop over found schedules
for id=1:rn
    if id == rn
        ch = chans(((id-1)*nchans+1):end);
    else
        ch = chans((id-1)*nchans+(1:nchans));
    end

    if farming > 0
        pid=myprocesses('farmit\(');
        while numel(pid)>=farming
            pid=myprocesses('farmit\(');
            pause(1);
        end
        [farmfile,jobname] = farmfilejobname('BICEP', 'rps_angfit', ...
            sprintf('ID_%04i_%d%d',id,floor(now),floor(rem(now,1)*1e7)));
        disp(jobname)
        farmit(farmfile, ...
            'rps_fit_angle_all(sch, p, rpsopt, dirname, fd, ch, id)', ...
            'jobname',jobname, 'queue',QUEUE, 'account','bicepdata_group', ...
            'mem',MEM, 'maxtime',MAXTIME,'overwrite',true, 'cleanout',true, ...
            'var', {'sch', 'p', 'rpsopt', 'dirname','fd', 'ch', 'id'});
    else
        rps_fit_angle_all(sch, p, rpsopt, dirname, fd, ch, id)
    end
end

