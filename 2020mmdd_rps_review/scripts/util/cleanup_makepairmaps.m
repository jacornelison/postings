function missind=cleanup_makepairmaps(n,refarm,checkmat,perphase)
% missind=cleanup_makepairmaps(n,refarm,checkmat,perphase)
% Identifies missing acsav's for simset n.  Will resubmit
% if refarm==1.
% Input args:
% n = simset of type XXXYYYZ
% Option refarm: Set to 1 to resubmit jobs to job queue.  Default is 0.
% Option to check .mat file for corruption... takes extra time at the moment.
% output missind is a vector of phases missing pairmaps
if nargin<4
  perphase=0;
end
if nargin<3
  checkmat=0;
end
if nargin<2
  refarm=0;
end

load(['simrunfiles/' n]);

% Need mapopt defaults if unspecified to get filename
if(~isfield(mapopt,'weight'))
  mapopt.weight=3;
end
if(~isfield(mapopt,'filt'))
  mapopt.filt='p3';
end
if(~isfield(mapopt,'coaddtype'))
  mapopt.coaddtype=0; % freq split
end
if(~isfield(mapopt,'proj'))
  mapopt.proj='radec';
end

% Get the filename extension
[gs,coaddtype,proj,filt]=gen_fname_ext(mapopt);

if(perphase)
  phases=taglist_to_phaselist(tags);
  list=phases.name;
else
  list=tags;
end

cc=0;
missind=[];

for i=1:length(list)
  mapopt.exmaptag=[list{i} '_'];
  mapopt.exname='';
  bad=0;
  filename=sprintf('%sfilt%s_weight%1d%s%s',mapopt.exmaptag,filt,mapopt.weight,gs,proj);
  % Check to see if pairmap is there:
  if ~exist(['pairmaps/' n '/' filename '.mat'],'file'); 
    disp(['pairmap ' int2str(i) ', ' filename ' does not exist.']); 
    cc=cc+1;
    bad=1;
    missind=[missind i];
  else
    if checkmat
      % if file exists, check to see if it is corrupted
      filestatus=check_matfile(['pairmaps/' n '/' filename '.mat']);
      if filestatus==0 
        disp(['pairmap ' int2str(i) ', ' filename ' is corrupted.']);
        cc=cc+1;
      end
    end
  end
  % if refarm == 1, send tag to be reprocessed.
  if bad & refarm
    disp('Refarming: ')
    farm_makepairmaps(n,1,i);
  end
end
disp([num2str(cc) ' of ' num2str(length(list)) ' are missing or corrupted.'])

function [gs,coaddtype,proj,filt]=gen_fname_ext(mapopt)

% append _gs if ground sub is on
if(mapopt.gs>0)
  gs='_gs';
  if(mapopt.gs>1)
    gs=[gs,num2str(mapopt.gs)];
  end
else
  gs='';
end

% only append the second digit to jacktype if non default coaddtype
if(mapopt.coaddtype>0)
  coaddtype=sprintf('%1d',mapopt.coaddtype);
else
  coaddtype='';
end

% append _ortho if sin projection
if(~strcmp(mapopt.proj,'radec'))
  proj=['_',mapopt.proj];
else
  proj='';
end

% concat sum/diff filts for output filename
if(iscell(mapopt.filt))
  filt=[mapopt.filt{:}];
else
  filt=mapopt.filt;
end

return
