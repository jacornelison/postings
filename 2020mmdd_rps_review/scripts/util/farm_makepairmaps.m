function [missingpairmaps,missingtags]=farm_makepairmaps(tags,mapopt,varargin)
% farm_makepairmaps(tags,mapopt)
%             or
% [missingpairmaps,missingtags]=farm_makepairmaps(tags,mapopt)
%
% For real data, mapopt.sernum should be a string of the form 'NNNNreal'
% For sims, mapopt.sernum should be a number of form NNNNXXXY
%
% farms reduc_makepairmaps on odyssey or spud
%
% tags = cell array of tags
% mapopt = mapping options. mapopt.sernum must be defined
%
% missingpairmaps and missingtags are cell arrays of the tags and filenames of
% pairmaps that will be made by call to farm_makepairmaps. If OnlyMissing=1 or 2, this
% list is a list of those pairmaps that are missing. 
%
% farm_makepairmaps may be called with the following optional arguments:
%  NTagsPerJob (default 1)
%  Group (default '/[experiment]_pipe_[username]')
%  Queue (default 'short_serial')
%  JobLimit (default 20)
%  MemRequire (default 6000)
%  Realizations (default does nothing, modifies input mapopt.sernum to farm sim
%                realizations from saved TODs)
%  OnlyMissing (default 0; if set to 1 only farm pairmaps that do not exist or have
%               size 0; if set to 2 return list of pairmaps that do not exist but do
%               not submit them for processing; overrides NTagsPerJob, sets to 1.)  
%  FarmJobs (default 1; if set to 1 use farmit to farm jobs
%                     ; if set to 0 do not farm jobs, run in interactive session (for
%                     testing) 
%
% ex. mapopt.sernum='1205real';
%     farm_makepairmaps(tags,mapopt,'Queue','normal_serial','MemRequire',15000);
%
%     mapopt.sernum='1205xxx1';
%     farm_makepairmaps(tags,mapopt,'Realizations',1:10);
%     mapopt.sernum='1205xxx2';
%     farm_makepairmaps(tags,mapopt,'Realizations',1:10,'OnlyMissing',1);

% Default arguments.
Ntagsperjob = 1;
QUEUE = 'short_serial';
JOBLIMIT = 20;
MEM=6000;
maxtime=120;
rlz=[];
onlymissing=0;
FarmJobs=true;
submit=true;
nice=[];

% Group name
username=whoami();
experiment=get_experiment_name;
GROUP=sprintf('/%s_pipe_%s',experiment,username);

% Parse varargin.
for i=1:length(varargin)

  % NTagsPerJob
  if strcmpi(varargin{i}, 'NTagsPerJob')
    if isnumeric(varargin{i+1})
      Ntagsperjob = varargin{i+1};
    end
  end

  % Group
  if strcmpi(varargin{i}, 'Group')
    if ischar(varargin{i+1})
      GROUP = varargin{i+1};
    end
  end
  
  % Queue
  if strcmpi(varargin{i}, 'Queue')
    if ischar(varargin{i+1})
      QUEUE = varargin{i+1};
    end
  end

  % JobLimit
  if strcmpi(varargin{i}, 'JobLimit')
    if isnumeric(varargin{i+1})
      JOBLIMIT = varargin{i+1};
    end
  end

  % MemRequire
  if strcmpi(varargin{i}, 'MemRequire')
    if isnumeric(varargin{i+1})
      MEM = varargin{i+1};
    end
  end
  
  % Max Time
  if strcmpi(varargin{i}, 'maxtime')
    if isnumeric(varargin{i+1})
      maxtime = varargin{i+1};
    end
  end

   % Realizations 
  if strcmpi(varargin{i}, 'Realizations')
    if isnumeric(varargin{i+1})
      rlz = varargin{i+1};
    end
  end

  % OnlyMissing
  if strcmpi(varargin{i}, 'OnlyMissing')
    if isnumeric(varargin{i+1}) | islogical(varargin{i+1})
      onlymissing = varargin{i+1};
    end
  end

  % FarmJobs
  if strcmpi(varargin{i}, 'FarmJobs')
    if isnumeric(varargin{i+1}) | islogical(varargin{i+1})
      FarmJobs = varargin{i+1};
    end
  end

  % submit the job or just create .mat file
  if strcmpi(varargin{i}, 'submit')
    if isnumeric(varargin{i+1}) | islogical(varargin{i+1})
      submit = varargin{i+1};
    end
  end

  if strcmpi(varargin{i}, 'nice')
    if isnumeric(varargin{i+1})
      nice = varargin{i+1};
    end
  end
end

% Uncomment for use with LSF queuing system:
%
%% Create group and limit the number of jobs
%if onlymissing~=2 ...
%    && ~strcmpi(QUEUE,'general') ...
%    && ~strcmpi(QUEUE,'serial_requeue')
%  system(['bgadd ' GROUP]);
%  system(sprintf('bgmod -L %d %s',JOBLIMIT,GROUP));
%end

% Break up tags into chunks for farming
nchunks=ceil(numel(tags)/Ntagsperjob);
for i=1:nchunks
  s=1+(i-1)*Ntagsperjob;
  e=min([i*Ntagsperjob,numel(tags)]);
  tagchunk{i}=tags(s:e);
end

% Make multiple mapopts for different realizations
if ~isempty(rlz)
  for i=1:numel(rlz)
    mapopt_new{i}=mapopt;
    mapopt_new{i}.sernum=sprintf('%s%03d%s',mapopt.sernum(1:4),rlz(i),mapopt.sernum(end));
  end
else
  mapopt_new={mapopt};
end

missingpairmaps={};
missingtags={};

for i=1:numel(tagchunk)
  
  for j=1:numel(mapopt_new)
    
    mapopt=mapopt_new{j};

    tags=tagchunk{i};
    
    % - farmfiles will have form, for example, farmfiles/XXXX/pairmapfilename_YYYZ.xxx
    % - jobnames will have form XXXX_pairmapfilename_YYYZ
    % e.g. farmfiles/1205/20100626H01_dk113_filtp3_weight3_gs_0051.mat &
    % 1205_20100626H01_dk113_filtp3_weight3_gs_0051
    % (The pairmap filename is for that of the first tag being farmed.)
    mapopt0=get_default_mapopt(mapopt);
    pairmapname=get_pairmap_filename(tags{1},mapopt0);
    [dum,pairmapbase]=fileparts(pairmapname);
    pairmapnames=get_pairmap_filename(tags,mapopt0);    

    clear d
    clear dofarm
    if onlymissing>0
      % print some info
      if mod(i,100)==1 & j==1
        disp(sprintf('checking pairmap files for tag block %d of %d',...
                     i, numel(tagchunk)));
      end
      % check to see if files exist.

      for k=1:numel(pairmapnames)
        d=dir(pairmapnames{k});
        if isempty(d) || d.bytes==0
          dofarm(k)=true;
        else
          dofarm(k)=false;
        end
      end
    else
       % make pairmap whether file exists or not
      dofarm=true(size(tags));
    end
    
    missingpairmaps={missingpairmaps{:},pairmapnames{dofarm}};
    missingtags={missingtags{:},tags{dofarm}};
    
    if any(dofarm) & onlymissing~=2
      % only farm these tags
      tags=tags(dofarm);
      
      farmfilebase=sprintf('%s_%s',pairmapbase,mapopt.sernum(5:end));
      farmfiledir=sprintf('farmfiles/%s',mapopt.sernum(1:4));
      if ~exist(farmfiledir,'dir')
        system(['mkdir ' farmfiledir]);
      end
      
      farmfile=sprintf('%s/%s.mat',farmfiledir,farmfilebase);
      outfile=sprintf('%s/%s.out',farmfiledir,farmfilebase);
      errfile=sprintf('%s/%s.err',farmfiledir,farmfilebase);
      
      JOBNAME=sprintf('%s_%s',mapopt.sernum(1:4),farmfilebase);
      
      cmd='reduc_makepairmaps(tags,mapopt)';
      
      if FarmJobs
        if ~exist(farmfile,'file')
          farmit(farmfile,cmd,...
               'group',GROUP,'queue',QUEUE,'mem',MEM,...
               'jobname',JOBNAME,'var',{'tags','mapopt'},...
               'errfile',errfile,'outfile',outfile,'maxtime',maxtime,...
               'nice',nice,'submit',submit);
        end
      else
        eval(cmd);
      end
    
    end % any(dofarm)
  end % loop over multiple mapopts / realizations
end % loop over tag chunks

if onlymissing
  disp(sprintf('%d pairmaps are missing',numel(missingpairmaps)));
end

return


