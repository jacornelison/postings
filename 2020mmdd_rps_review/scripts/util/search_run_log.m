function [gcp_messg,t]=search_run_log(string,start_day,end_day,sched_flag,logdir)
% [gcp_messg,t]=search_run_log(string,start_time,end_time,sched_flag,logdir)
%
% search GCP logs for a string and return time and full GCP log message; searches are
% case sensitive (much faster)
% 
% if schedule flag is set (default no) search GCP logs for specified string and return
% schedule name and start/stop times
%
% string is name of string or schedule name to search for, default '*' (can use UNIX
% string is name of string or schedule name to search for, default '*' (can use UNIX
% format wildcards)  
% start_day and end_day take form 'yyyymmdd' or 'yyyymmdd 01:37:15'
% logdir is directory containing GCP logs (default 'log')
%
% ex. [gcp_message]=find_gcp_messages('A control client');
% ex. [gcp_message]=find_gcp_messages('A control client','20110521','20110523');
% ex. [gcp_message]=find_gcp_messages('starpoint',[],[],1);

if(~exist('start_day','var'))
  start_day=[];
end
if(isempty(start_day))
  start_day='20000101';
end

if(~exist('end_day','var'))
  end_day=[];
end
if(isempty(end_day))
  end_day='20200101';
end

if(~exist('sched_flag','var'))
  sched_flag=[];
end
if(isempty(sched_flag))
  sched_flag=0;
end

if(~exist('logdir','var'))
  logdir=[];
end
if(isempty(logdir))
  logdir='log';
end

% put start/end day into correct format
if(numel(start_day)==8)
  start_day=[start_day,' 00:00:00'];
end
if(numel(end_day)==8)
  end_day=[end_day,' 00:00:00'];
end

% list log files
d=dir([logdir,'/20*_*.log']);
for i=1:numel(d)
  fname{i}=d(i).name;
  fdate(i)=datenum(d(i).name(1:8),'yyyymmdd');
end

% what logfiles do we want to search? search starting from requested start day minus
% one to end day
date1=datenum(start_day(1:8),'yyyymmdd');
date2=datenum(end_day(1:8),'yyyymmdd');
grepind=find(fdate>=(date1-1) & fdate<=date2);


% exit if there are no log files
if(isempty(grepind))
  t=[];
  gcp_messg=[];
  disp('no log files found within date range')
  return
end

% construct file list to pass to grep
filelist='';
for i=1:numel(grepind)
  filelist=[filelist,' ',fullfile(logdir,fname{grepind(i)})];
end

% now search the log files using unix grep
searchstring=sprintf('grep -H "%s" %s',string,filelist);
[status,result]=system_safe(searchstring);

% now parse the result
delimiter=fullfile(logdir,'/20\d{6,6}_\d{6,6}.log:');
x=regexp(result,delimiter,'split');
x=x(~strcmp(x,''));
for i=1:numel(x)
  t{i}=x{i}(1:15);
  gcp_messg{i}=x{i}(17:end);
end

% exit if there is nothing returned by grep
if isempty(x)
  t=[];
  gcp_messg=[];
  disp('no log messages found within date range')
  return
end

% only keep the files within the requested date/time range
tdatenum=datenum(t,'yymmdd HH:MM:SS');
keepind=find(tdatenum>=datenum(start_day,'yyyymmdd HH:MM:SS') ...
             & tdatenum<=datenum(end_day,'yyyymmdd HH:MM:SS'));

% exit if there are no log messages
if(isempty(keepind))
  t=[];
  gcp_messg=[];
  disp('no log messages found within date range')
  return
else
  t=t(keepind);
  tdatenum=tdatenum(keepind);
  gcp_messg=gcp_messg(keepind);
end

if(sched_flag);
  % parse output log messages to find schedules
  % keep schedules that have a clear "schedule" ,"Starting schedule", "Exiting
  % schedule" signature in the GCP logs and has a <5 second lag between "schedule" and
  % "Starting schedule:"
  
  sched_ind=find(strncmp(gcp_messg,'schedule ',9));
  starting_ind=find(strncmp(gcp_messg,'Starting schedule:',18));
  exiting_ind=find(strncmp(gcp_messg,'Exiting schedule:',17));
  
  % l counts the number of completed schedule runs.
  l=0;
  for i=1:length(sched_ind)
    j=find(starting_ind>sched_ind(i),1,'first');
    k=find(exiting_ind>sched_ind(i),1,'first');
    if ~isempty(j) & ~isempty(k) & transpose((tdatenum(starting_ind(j))-tdatenum(sched_ind(i)))*24*3600<5)
      l=l+1;
      gcp_messg_out{l}=gcp_messg{sched_ind(i)}(10:end-1); % schedule name (new line charactor at the end is removed)
      t_out{1,l}=t{starting_ind(j)}; % schedule start time
      t_out{2,l}=t{exiting_ind(k)}; % schedule end time
    end
  end

  if l==0
    % no completed schedules are found.
    t=[];
    gcp_messg=[];
    disp('No schedules found within date range.');
    return
  end

  t=t_out;
  gcp_messg=gcp_messg_out;
end

return
