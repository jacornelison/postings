function [missingmaps,allmaps,missingjacks,missingdp]=farm_coaddpairmaps(tags,coaddopt,varargin)
% farm_coaddpairmaps(tags,coaddopt)
%             or
% [missingmaps,allmaps,missingjacks,missingdp]=farm_coaddpairmaps(tags,coaddopt)
%
% For real data, coaddopt.sernum should be a string of the form 'NNNNreal'
% For sims, coaddopt.sernum should be a number of form NNNNXXXY
%
% farms reduc_coaddpairmaps on odyssey or spud
%
% tags = cell array of tags
% coaddopt= coadd options
%
% missingmaps is a cell array of the filenames of coadded
% maps that will be made by call to farm_coaddpairmaps. If OnlyMissing=1 or 2, this
% list is a list of those maps that are missing. 
%
% farm_coaddpairmaps may be called with the following optional arguments:
%  Group (default '/[experiment]_pipe_[username]')
%  Queue (default none; uses farmit() default.)
%  JobLimit (default 20)
%  MemRequire (default 2000)
%  Realizations (default does nothing, modifies input mapopt.sernum to farm map making
%                from simulated pairmaps)
%  FarmJobs (default 1, farm jobs. If 0, run them in the interactive terminal; good for
%            testing) 
%  FarmJacksSeparately (default 1, needed if coadding many many tags, but annoying if
%                       coadding many realizations of a few tags)
%  OnlyMissing (default 0; if set to 1 only farm maps that do not exist or have
%               size 0; if set to 2 return list of maps that do not exist but do
%               not submit them for processing)  
%
%  SplitSubmission (default 0): split the tag list into days (these are equal or longer 
%                               than a raster) and coadd over these splits
%                               when all maps have been created rerun reduc_coaddpairmaps
%                               and the splits will be coadded.
%                             : = 1 or true : split per day, e.g. 362 splits for 3yrs of B2
%                             : = N, i.e. 50 to split into 50 coadds, 
%                                    maximum the number of days
%
%  UseCompiled - (default 0): Use the standard, non-compiled version of reduc_coaddpairmaps
%                           : = 1 or true : use compiled version of reduc_coaddpairmaps
%
% ex. coaddopt.daughter='a';
%     farm_coaddpairmaps(tags,coaddopt,'MemRequire',6000);
%
%     coaddopt.sernum='1205xxx1';
%     farm_coaddpairmaps(tags,coaddopt,'Realizations',1:10);
%     coaddopt.sernum='1205xxx2';
%     farm_coaddpairmaps(tags,coaddopt,'Realizations',1:10,'OnlyMissing',1);

% Default coaddopts
if ~isfield(coaddopt,'tags');
  coaddopt.tags=tags;
end
coaddopt=get_default_coaddopt(coaddopt);

% Default arguments.
QUEUE = {};
JOBLIMIT = 20;
if strcmp(coaddopt.sernum(5:end),'real')
  MEM = 7000;
else
  MEM = 2500;
end
rlz=[];
OnlyMissing=false;
FarmJobs=true;
FarmJacksSeparately=true;
SplitSubmission=false;
UseCompiled=false;
maxtime=[];
nice=[];
submit=true;

% Group name
username=whoami();
experiment=get_experiment_name;
GROUP=sprintf('/%s_pipe_%s',experiment,username);

% Parse varargin.
for i=1:length(varargin)

  % Group
  if strcmpi(varargin{i}, 'Group')
    if ischar(varargin{i+1})
      GROUP = varargin{i+1};
    end
  end
  
  % Queue
  if strcmpi(varargin{i}, 'Queue')
    if ischar(varargin{i+1}) && ~isempty(varargin{i+1})
      QUEUE = {'queue', varargin{i+1}};
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

   % Realizations 
  if strcmpi(varargin{i}, 'Realizations')
    if isnumeric(varargin{i+1})
      rlz = varargin{i+1};
    end
  end

  % OnlyMissing
  if strcmpi(varargin{i}, 'OnlyMissing')
    if isnumeric(varargin{i+1}) | islogical(varargin{i+1})
      OnlyMissing = varargin{i+1};
    end
  end

  % FarmJobs
  if strcmpi(varargin{i}, 'FarmJobs')
    if isnumeric(varargin{i+1}) | islogical(varargin{i+1})
      FarmJobs = varargin{i+1};
    end
  end
  
  % FarmJacksSeparately
  if strcmpi(varargin{i}, 'FarmJacksSeparately')
    if isnumeric(varargin{i+1}) | islogical(varargin{i+1})
      FarmJacksSeparately = varargin{i+1};
    end
  end

  % SplitSubmission
  if strcmpi(varargin{i}, 'SplitSubmission')
    if isnumeric(varargin{i+1}) | islogical(varargin{i+1})
      SplitSubmission = varargin{i+1};
    end
    % when SplitSubmission is used the vargin are going to 
    % be passed forward to another call to farm_coaddpairmaps.
    % get rid of the SplitSubmission argument:
    varargin{i}=[];
    varargin{i+1}=[];
  end
  
  % UseCompiled
  if strcmpi(varargin{i}, 'UseCompiled')
    if isnumeric(varargin{i+1}) | islogical(varargin{i+1})
      UseCompiled = varargin{i+1};
    end
  end
  
  % submit the job or just create .mat file
  if strcmpi(varargin{i}, 'submit')
    if isnumeric(varargin{i+1}) | islogical(varargin{i+1})
      submit = varargin{i+1};
    end
  end
 
  % maxtime
  if strcmpi(varargin{i}, 'maxtime')
    if isnumeric(varargin{i+1})
      maxtime = varargin{i+1};
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
%global bjobgroup;
%if OnlyMissing~=2
%  % just do it when we have to
%  if ~strcmp(bjobgroup,[GROUP,num2str(JOBLIMIT)])
%    bjobgroup=[GROUP,num2str(JOBLIMIT)];
%    system(['bgadd ' GROUP]);
%    system(sprintf('bgmod -L %d %s',JOBLIMIT,GROUP));
%  end
%end

missingmaps={};
missingjacks='';
missingdp=[];
allmaps={};

if SplitSubmission
  
  % when rlz=[], real data, here a kludge to get into the for
  % loop below anyways:
  if isempty(rlz) 
    rlz=-1;
  end
  
  % there is a convenience loop later over realizations
  % it is easier to perform that loop over realizations
  % here, when SplitSubmission is used  
  for crlz=rlz
    % make a copy of the coaddopt to change it 
    % in the loops:
    coaddopt_split = coaddopt;
    if crlz==-1
      varargin = {varargin{:},'Realizations',[]};
    else
      varargin = {varargin{:},'Realizations',crlz};
    end
    
    % check which maps are missing:
    if OnlyMissing>0
      % run farm_coaddpairmaps with out submitting jobs:
      varargintmp = {varargin{:},'OnlyMissing',[2]};
      [mismaps,newmap,misjacks,missingdp]=farm_coaddpairmaps(tags,coaddopt_split,varargintmp{:});
      % ... which yields a list of what is missing.
      if length(mismaps)==0 
        allmaps = [allmaps, newmap];
        continue;
      end
      % if there are missing maps, just
      % do the corresponding jacks:
      coaddopt_split.jacktype=misjacks;
      coaddopt_split.deproj=missingdp;
    end

    % split the tag list into days, correpsonding to at least
    % a raster worth of data, i.e. 3yrs of B2 are ~360 splits,
    % if SplitSubmission > 1 the number of splits is reduced
    split=split_taglist(tags,SplitSubmission);

    % doing this handels signflip sequence split up in the case
    % of signflip noise:
    coaddopt_split.tags=tags;
    
    % record the missing maps, this will allow us to decide
    % when to run reduc_coadd_coadds:
    missingmaps = {};

    % loop over the splits:
    for ii=1:length(split.name)
      % use the name of the current split (e.g. the datestring)
      % to define the daughter coadd name:
      cname = split.name(ii);
      coaddopt_split.daughter=char(strcat(coaddopt.daughter,'s',num2str(SplitSubmission),cname));
      % coadd over the split
      [mm,allmaps] = farm_coaddpairmaps(split.tags{ii},coaddopt_split,varargin{:});
      missingmaps = union(missingmaps,mm);
    end

    % when no more maps are missing, run reduc_coadd_coadds:
    if length(missingmaps)==0
      % the allmaps variable holds the filenames created during the last
      % split (ends on jack0, jack1 etc).
      for ii=1:numel(allmaps)
        % Explicitly enumerate all component maps we expect to coadd over. Do
        % this rather than use a wildcard to avoid race conditions where
        % reduc_coadd_coadds() may have deleted some maps already (either
        % with two nearly-simultaneously running jobs, or a resubmission of
        % a job which was killed during.
        cmap = cellfunc(@(s) strrep(allmaps{ii}, cname{:}, s), split.name);

        % Then run reduc_coadd_coadds with the file list:
        sernum=strrep(cmap{1}(6:14),'/','');
        cmd = 'reduc_coadd_coadds(cmap,sernum,coaddopt.daughter,1)';
        % Remove tag split from file name to generate the farmfile name.
        farmfile = strrep(['farmfiles/',cmap{end}(6:end)], cname{:}, '');
        [dum,mapbase]=fileparts(farmfile);
        JOBNAME=sprintf('%s_%s',coaddopt.sernum(1:4),mapbase);
        JOBNAME=strrep(JOBNAME,'.mat','');
        % TODO: The above few lines should probably be updated to make use of
        %       farmfilejobname(), but during simset production is not
        %       appropriate time to try that out...
        if OnlyMissing~=2
          if FarmJobs
            if ~exist_file(farmfile)
              farmit(farmfile,cmd,'group',GROUP,...
                'jobname',JOBNAME,QUEUE{:},'mem',MEM,...
                'var',{'cmap','sernum','coaddopt.daughter'},...
                'maxtime',60,'nice',nice,'submit',submit);
            else
              display(['Skipping over existing job: ',farmfile])
            end
          else
            eval(cmd);
          end
        end %OnlyMissing
      end
    end% if length(missingmaps)==0
  end% rlz
  
  % SplitSubmission works recursevely. So finish here, when
  % this option is used:
  return
end % SplitSubmission

% Make multiple coaddopts for different realizations
if ~isempty(rlz)
  for i=1:numel(rlz)
    coaddopt_new{i}=coaddopt;
    coaddopt_new{i}.sernum=sprintf('%s%03d%s',coaddopt.sernum(1:4),rlz(i),coaddopt.sernum(end));
  end
else
  coaddopt_new={coaddopt};
end

% loop over realizations
for j=1:numel(coaddopt_new)

  coaddopt=coaddopt_new{j};
  
  if FarmJacksSeparately
    mj=coaddopt.jacktype;
  else
    mj=1;
  end

  % print some info
  if ~isempty(rlz) & OnlyMissing>0
    disp(sprintf('checking map files for realization %d (%d of %d)',...
                 rlz(j), j, numel(rlz)));
  elseif OnlyMissing>0
    disp('checking map files for realization 1 (1 of 1)');
  end
  
  % loop over jacktypes
  for k=mj
  
    if FarmJacksSeparately
      coaddopt.jacktype=k;
    end
    
    % - farmfiles will have form, for example, farmfiles/XXXX/mapfilename.xxx
    % - jobnames will have form XXXX_mapfilename
    % e.g. farmfiles/1205/0051_a_filtp3_weight3_gs_jack0.mat &
    % 1205_0051_a_filtp3_weight3_gs_jack0
    mapname=get_map_filename(coaddopt);
    if ~iscell(mapname)
      mapname={mapname};
    end
    allmaps={allmaps{:},mapname{:}};
    
    clear d
    clear dofarm
    if OnlyMissing>0

      % check to see if files exist.
      for l=1:size(mapname,1)
        for kk=1:size(mapname,2)
          dofarm(l,kk)=~exist_file(mapname{l,kk});
        end
      end

    else
      % make pairmap whether file exists or not
      dofarm=repmat(true,size(mapname));
    end
   
    if any(dofarm(:))
      missingmaps={missingmaps{:},mapname{dofarm}};
      missingjacks=[missingjacks,coaddopt.jacktype(any(dofarm,2))];
      missingdp=[missingdp' coaddopt.deproj(any(dofarm,1),:)']';
    end

    if any(dofarm(:)) & OnlyMissing~=2

      % excise coaddopt.jacktypes that exist. This only does anything of
      % numel(coaddopt.jacktype)>1, i.e. FarmJacksSeparately=0.
      % also excise dp that exist.
      coaddopt.jacktype=coaddopt.jacktype(any(dofarm,2));
      coaddopt.deproj=coaddopt.deproj(any(dofarm,1),:);

      farmfiledir=sprintf('farmfiles/%s',coaddopt.sernum(1:4));
      if ~exist(farmfiledir,'dir')
        system(['mkdir ' farmfiledir]);
      end
      
      [dum,mapbase]=fileparts(mapname{1});
      farmfile=sprintf('%s/%s.mat',farmfiledir,mapbase);
      outfile=sprintf('%s/%s.out',farmfiledir,mapbase);
      errfile=sprintf('%s/%s.err',farmfiledir,mapbase);
      
      JOBNAME=sprintf('%s_%s',coaddopt.sernum(1:4),mapbase);
      
      % if UseCompiled is false, run reduc_coaddpairmaps as normal
      if ~UseCompiled
        cmd='reduc_coaddpairmaps(tags,coaddopt)';
      % otherwise use compiled version of reduc_coaddpairmaps
      else
        cmd='reduc_coaddpairmaps ';
      end
      if FarmJobs
        % skip over farmfiles that exist and
        % assume the job is still running...
        if ~exist_file(farmfile)
          farmit(farmfile,cmd,...
                'group',GROUP,QUEUE{:},'mem',MEM,...
                'jobname',JOBNAME,'var',{'tags','coaddopt'},...
                'errfile',errfile,'outfile',outfile,'usecompiled',UseCompiled,...
                'maxtime',maxtime,'nice',nice,'submit',submit);
        else
          display(['Skipping over existing job: ',farmfile])
        end
      else
        eval(cmd);
      end
      
    end % any(dofarm)
  end % loop over jacktypes
end % loop over multiple coaddopts / realizations

% the missing jacks have been concatenated for multiple realizations:
missingjacks=unique(missingjacks);
missingdp=unique(missingdp,'rows');

if OnlyMissing
  disp(sprintf('%d maps are missing',numel(missingmaps)));
end
  


return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ph=split_taglist(tags,nsplits)
%  function days=split_taglist(tags,nsplits)

% that is the current maximum deprojection timescale
[ph,idx]=taglist_to_phaselist(tags,'inclusive','raster');

if nsplits>length(ph.tags) | nsplits==1
  nsplits=length(ph.tags);
end
stepw= round(length(ph.tags)/nsplits);

if stepw>1
  count=1;
  start=1;
  stop=start+stepw-1;
  while stop<=length(ph.tags);
    ph2.name{count}=ph.name{start};
    ph2.tags{count}=tags(idx>=start & idx<=stop);
    start=stop+1;
    stop=start+stepw-1;
    count=count+1;
  end
  if stop~=length(ph.tags) & start<=length(ph.tags)
    ph2.name{count}=ph.name{start};
    ph2.tags{count}=tags(idx>=start & idx<=length(ph.tags));
  end
%    sumtags = 0;
%    for ii=1:length(ph2.tags)
%      sumtags=sumtags+length(ph2.tags{ii});
%    end
%    display(['Sum over tags in the split: ',num2str(sumtags)]);
  ph=ph2;
end

return


