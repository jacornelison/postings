function set_symlinks(experiment,dir_ext,localdirs,simflag)
% set_symlinks(experiment,dir_ext,localdirs,simflag)
% set up symlinks for running pipeline on B2 or Keck data
%
% the following links are created in the analysis directory:
% arc
% aux_data
% data
% log
% scratch
% input_maps
% pairmaps
% maps
% pntscrmask
% reducplots
% aps
% noispec
%
% experiment is either 'bicep2' or 'keck'
%
% simflag = {0 (default), 1}
%           if set, "data" symlink points to temporary file system on odyssey for speed
%           during simulations
%
% dir_ext: 
%   if specified symlink [linkname] will have as target ./[linkname]+dir_ext  
%   i.e. default behavior creates link pairmaps->/data/[experiment]/pipeline/pairmaps/,
%   setting dir_ext='_sandbox' creates link pairmaps->./pairmaps_sandbox
%
% localdirs:  
%   cell string array {[linkname1],[linkname2],...} specifying which symlinks have as
%   target ./[linkname]+dir_ext
%
%   default: {'scratch','pairmaps','maps','reducplots','aps'};
%
%   i.e setting local={'data','aps'} would make all links as normal except
%       aps->./aps_sandbox/
%       data->./data_sandbox/
%   which would allow testing of reduc_initial and testing of reduc_makeaps using official
%   maps 
%
% simflag = {0 (default), 1}
%           if set, "data" symlink points to temporary file system on odyssey for speed
%           during simulations

%
% ex.
% set_symlinks('keck')
% set_symlinks('bicep2')
% set_symlinks('bicep2','_sandbox')
% set_symlinks('bicep2','_local',{'aps'})
% set_symlinks('bicep2',[],[],1)

if(~exist('simflag','var'))
  simflag=[];
end
if(isempty(simflag))
  simflag=0;
end

if(~exist('dir_ext','var'))
  dir_ext=[];
end

if(~exist('localdirs','var'))
  localdirs=[];
end
if(isstr(localdirs))
  localdirs={localdirs};
end

% only run this function from bicep2_analysis or keck_analysis directories
p=pwd;
k=strfind(pwd,'/');
currentfolder=p(k(end)+1:end);
if(~(strcmp(currentfolder,'bicep2_analysis') | strcmp(currentfolder,'keck_analysis')));
  error('Run only from bicep2_analysis or keck_analysis');
end

% what machine are we running on?
host=whichhost;

% these are the symlink names we will create
linknames{1}='arc';
linknames{2}='aux_data';
linknames{3}='data';
linknames{4}='log';
linknames{5}='scratch';
linknames{6}='input_maps';
linknames{7}='pairmaps';
linknames{8}='maps';
linknames{9}='pntsrcmask';
linknames{10}='reducplots';
linknames{11}='farmfiles';
linknames{12}='tuning';
linknames{13}='tuning_plots';
linknames{14}='reduc_plotcuts_online';
linknames{15}='starpoint_plots';
linknames{16}='matrixdata';

% get username
username=whoami();

switch host
 
 case 'spud'
 
  switch experiment
   case 'keck'
    dataloc='data';
   case 'bicep2'
    dataloc='data2';
  end
  % directories (no trailing slashes)
  datdir{1}=sprintf('/%s/%sdaq/arc',dataloc,experiment); % arc files
  datdir{2}=sprintf('/%s/%s/aux_data',dataloc,experiment); % aux_data
  datdir{3}=sprintf('/%s/%s/pipeline/data',dataloc,experiment); % reduced .mat data
  datdir{4}=sprintf('/%s/%sdaq/log',dataloc,experiment); % gcp log files
  datdir{5}=sprintf('/%s/%s/pipeline/data',dataloc,experiment); % scratch
  datdir{6}=sprintf('/%s/%s/input_maps',dataloc,experiment); % simulation maps
  datdir{7}=sprintf('/%s/%s/pipeline/pairmaps',dataloc,experiment); % pairmaps
  datdir{8}=sprintf('/%s/%s/pipeline/maps',dataloc,experiment); % maps
  datdir{9}=sprintf('/%s/%s/pipeline/pntsrcmask',dataloc,experiment); % point source masks
  datdir{10}=sprintf('/%s/%s/pipeline/reducplots',dataloc,experiment); % reducplots
  datdir{11}=sprintf('/data/%s/farmfiles',username); % farmfiles
  
 case 'odyssey';

  switch experiment
   case 'keck'
    dataloc='bicepfs2';
   case 'bicep2'
    dataloc='bicepfs1';
  end
  % directories (no trailing slashes)
  datdir{1}=sprintf('/n/%s/%s/%sdaq/arc',dataloc,experiment,experiment); 
  datdir{2}=sprintf('/n/%s/%s/%s_aux_data',dataloc,experiment,experiment);
  datdir{3}=sprintf('/n/%s/%s/pipeline/data',dataloc,experiment);
  datdir{4}=sprintf('/n/%s/%s/%sdaq/log',dataloc,experiment,experiment);
  datdir{5}=sprintf('/n/panlfs2/bicep/%s/pipeline/data',experiment);
  datdir{6}=sprintf('/n/%s/%s/pipeline/input_maps',dataloc,experiment);
  datdir{7}=sprintf('/n/%s/%s/pipeline/pairmaps',dataloc,experiment);
  datdir{8}=sprintf('/n/%s/%s/pipeline/maps',dataloc,experiment);
  datdir{9}=sprintf('/n/%s/%s/pipeline/pntsrcmask',dataloc,experiment);
  datdir{10}=sprintf('/n/%s/%s/pipeline/reducplots',dataloc,experiment); 
  datdir{11}=sprintf('/n/panlfs2/bicep/%s/farmfiles',username);
  datdir{12}=sprintf('/n/%s/%s/%sdaq/tuning',dataloc,experiment,experiment); 
  datdir{13}=sprintf('/n/%s/www/%s/tuning_plots',dataloc,experiment); 
  datdir{14}=sprintf('/n/%s/www/%s/reduc_plotcuts',dataloc,experiment); 
  datdir{15}=sprintf('/n/%s/www/%s/starpoint_plots',dataloc,experiment);
  datdir{16}=sprintf('/n/%s/%s/pipeline/matrixdata',dataloc,experiment); 
  
  if(simflag)
    % for sims on odyssey use real data stored on temporary file system for speed
    datdir{3}=sprintf('/n/panlfs2/bicep/%s/pipeline/data',experiment);
  end
  
end

% point to local directories if dir_ext is specified
if(isstr(dir_ext));
  for i=1:numel(localdirs)
    ind=find(strcmp(linknames,localdirs{i}));
    if(~isempty(ind))
      datdir{ind}=[linknames{ind},dir_ext];
    end
  end
end

% remove existing symlinks and create new ones
for i=1:numel(datdir)

  % make double sure linknames don't have trailing slashes
  if(strcmp(linknames{i}(end),'/'))
    linknames{i}=linknames{i}(1:end-1);
  end

  % only delete symlinks, we don't want to wipe out anyone's data
  [dum,x]=unix(['ls -ld ',linknames{i}]);
  if(~isempty(strfind(x,'->')))
    % delete
    cmd=['rm -f ',linknames{i}];
    unix(cmd);
  end

  % create symlinks if there doesn't already exist a directory or file with the same name
  if exist(fullfile(pwd,linknames{i}),'file')
    disp(sprintf('%s exists and is not a symlink; cannot create symlink',linknames{i}));
  else
    cmd=['ln -s ',datdir{i},' ',linknames{i}];
    unix(cmd);
  end

end

return


