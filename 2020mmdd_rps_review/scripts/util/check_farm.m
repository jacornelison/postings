% [done_list not_done_list failed_list kill_list]=check_farm(fnames)
%
% Checks status of farm jobs using farmit.
% Input fnames is a cell array of farmit temporary file names.
% The same file names are returned as output, but sorted into
% four categories:
%   done_list : those that have completed successfully
%   not_done_list : those that are running or pending
%   failed_list : those that have failed
%   kill_list : suspended, over time, etc.
%
% check_farm(fnames,doresub,dokill)
%
% If doresub=1, any jobs that have failed will be resubmitted.
% These are then listed under not_done_list rather than failed_list.
%
% If dokill=1, any suspended, over time, etc. jobs will be
% killed.  These are then listed under failed_list.
%
% If dokill=1 and doresub=1, failed_list and kill_list will
% both be empty, and any jobs that are killed / failed are
% included in not_done_list.
function [done_list not_done_list failed_list kill_list] = check_farm(fnames,doresub,dokill)

% By default, don't resubmit anything
if nargin<2 || isempty(doresub)
  doresub=false;
end
% By default, don't kill anything
if nargin<3 || isempty(dokill)
  dokill=false;
end

done_list     ={};
not_done_list ={};
failed_list   ={};
kill_list     ={};

% Get status of jobs using farmit
stat=farmit(fnames,'status');
for i=1:length(stat)

  % Take "no file" to mean job completed & temp file was deleted
  if isempty(stat{i})
    done_list=[done_list,{fnames{i}}];
    continue
  end

  % Interpretation of LSF codes and fake codes NOJOB, OVERTIME
  % Note that Matlab exits and returns 0 regardless of whether
  % the function completed successfully.  So our failed jobs
  % get DONE, not EXIT.
  switch(stat{i})
    case {'RUN','PEND'},
      not_done_list=[not_done_list,{fnames{i}}];
    case {'NOJOB','EXIT','DONE'},
      failed_list=[failed_list,{fnames{i}}];
    case {'SUSP','OVERTIME'},
      kill_list=[kill_list,{fnames{i}}];
    otherwise,
      disp(['Unrecognized job status. Treating as failed.']);
      stat{i}
      failed_list=[failed_list,{fnames{i}}];
  end
end

% Kill, if requested.
if dokill && ~isempty(kill_list)
  disp(['Killing ' num2str(length(kill_list)) ' jobs.']);
  farmit(kill_list,'kill');
  failed_list=[failed_list,kill_list];
  % At this point, we could wait and check again that the jobs
  % have really died; if not, kill 'em harder.
  kill_list={};
end

% Resubmit, if requested.
if doresub && ~isempty(failed_list)
  disp(['Resubmitting ' num2str(length(failed_list)) ' jobs.']);
  for i=1:length(failed_list)
    farmit(failed_list{i},'resubmit');
  end
  not_done_list=[not_done_list,failed_list];
  failed_list={};
end

return
