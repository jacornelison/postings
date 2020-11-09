function [tags,scandecs]=get_scandecs(indec,sdopt)
% [tags,scandecs]=get_scandecs(indec,sdopt)
%
% OUTPUT ARGUMENTS
% if only tags output is specified, outputs only tags
% corresponding to days where indec declination is covered
%
% if scandecs option is specified, outputs scandecs struct for 'all' days
% scandecs.tags = days 
%         .dec  = corresponding dec
%         .el   = corresponding el
%
% INPUT ARGUMENTS
% indec specifies the dec required, empty uses util/get_src for dec
%
% sdopt.generate = 0 - load last saved scandecs.mat
%                = 1 - generate a new scandecs for 'all' days
%
%      .save     = 0 - uses new scandecs but doesnt save
%                = 1 - saves and exits after generating new scandecs
%                  sdopt.type string appended to filename
%
%      .hitchan  = 'center' - days indec is covered by only center pixel
%                = 'any'    - days indec is covered by any pixel
%      
%      .intags   = opt - same as get tags (ex. 'mel67')
%                  specifies which subset of days to query
%
%      .src      = 0 - use indec
%                = integer > 0 - get source coords from util/get_src
%
%      .type     ='quad05' - load scandecs file for all of quad 05 days
%                ='quad06' - quad 06 days (as of 10/16 - ~exist)
%                  really only used to specify which data/scandecs*.m
%                  file to load.  not the same as get_map_defn input 'type'
%
% 06/10/05 - data/scandecs_quad05.mat generated, all 04-05 season days
%


if(~exist('indec','var'))
  indec=[];
end
if(~exist('sdopt','var'))
  sdopt=[];
end
% set default sdopt flags
sdopt=set_sdopt(sdopt,indec);

% generate new scandecs array, 
% basically for one time use only or if scandecs file is lost
if(sdopt.generate==1)
  
  tags=get_tags('all');

  for i=1:length(tags)

    eval(sprintf('load data/%s_tod.mat',tags{i}));
    mapind=make_mapind(d,fs);

    scandecs.tag{i}=tags{i};
    scandecs.dec(i,:)=[min(d.dec(mapind)),max(d.dec(mapind))];
    scandecs.el(i,:)=[min(d.el(mapind)),max(d.el(mapind))];

  end
  
  if(sdopt.save==1)
    eval(sprintf('save data/scandecs_%s.mat',sdopt.type));
    return
  end  

% load existing scandecs array
else

  eval(sprintf('load data/scandecs_%s',sdopt.type));
  
end

% if specified, aquire source coords from util/get_src
% otherwise use existing indec
if(sdopt.src>0)
  src=get_src;
  indec=src.dec(sdopt.src);
end

%locate those days that covered the quasar range of dec
switch sdopt.hitchan
  case 'any'
    hitind=find(((scandecs.dec(:,1)-0.77)<indec)&((scandecs.dec(:,2)+0.77)>indec));
  case 'center'
    hitind=find((scandecs.dec(:,1)<indec)&(scandecs.dec(:,2)>indec));
end

tags=get_tags('all');
intags=get_tags(sdopt.intags);

tags=tags(hitind);

tags=setdiff(tags,setdiff(tags,intags));

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function sdopt=set_sdopt(sdopt,indec)

%set sdopt to defaults
if(isempty(indec))
  sdopt.src=1;
end
if(~isfield(sdopt,'generate'))
  sdopt.generate=0;
end
if(~isfield(sdopt,'save'))
  sdopt.save=0;
end
if(~isfield(sdopt,'hitchan'))
  sdopt.hitchan='any';
end
if(~isfield(sdopt,'intags'))
  sdopt.intags='mel67';
end
if(~isfield(sdopt,'src'))
  sdopt.src=0;
end
if(~isfield(sdopt,'type'))
  sdopt.type='quad05';
end
