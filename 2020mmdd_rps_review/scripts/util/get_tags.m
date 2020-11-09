function tags=get_tags(tagsetname,flag)
% tags=get_tags(tagsetname,flag)
%
% Get a list of usable tags according to specified option
% 
% if flag='new' then only tags for which no _tod_uncal.mat file exists will be
% returned - used to update the reduction
% if flag='has_tod' only tags for which the _tod.mat file exists will
% be returned - used to avoid tags which crashed reduc_initial
% flag can also be a cell array of multiple requirements
%      
% e.g. tags=get_tags('cmb2010')
%      tags=get_tags('cmb2010','new')
%      tags=get_tags('cmb2010',{'has_tod','has_cuts'})
% Valid data selections:
%      'allcmb', 'allgal', 'all'
%      ['cmb' YYYY], ['cmb' YYYY 'full'], ['gal' YYYY], ['all' YYYY]
%      subset2010b2, galsubset2010b2


if(~exist('tagsetname','var'))
  tagsetname='cmb2010';
end

if(~exist('flag','var'))
  flag='all';
end

if ischar(flag)
  flag={flag};
end

% get experiment name
expt=get_experiment_name;

auxpath='aux_data';
taglist='tag_list.csv';
[pt,kt]=ParameterRead(fullfile(auxpath,taglist));

% find desired year (assumed to be the only 4-digit sequence)
targetyearstr = regexp(tagsetname, '(\d{4})', 'match', 'once');
if ~isempty(targetyearstr)
  targetyear = str2double(targetyearstr);
else
  targetyearstr = '+Inf';
  targetyear = Inf;
end

% exclude bad tags as identified by the 'good' field in the tag list:
if isfield(pt,'good')
  % get values - set empty values to 1
  pt.good=strtrim(pt.good);
  good=strcmp(pt.good,'') | strcmp(pt.good,'1');
  % downselect to exclude tags for which the 'good' field is zero
  pt = structcut(pt, logical(good));
  clear good
end

% construct convenient phase/year/date/set var
[date,phase] = parse_tag(pt.tag);
date = cvec(date);
phase = cvec(phase);
year = fix(date / 1e4);
pt.set=strtrim(pt.set);
setn=zeros(size(pt.set));
setn(~strcmp(pt.set,''))=str2num(char(pt.set(~strcmp(pt.set,''))));
type=pt.type;

% Define season start dates - format as YYYYMMDD integers
%   These may need to depend on experiment (B2/Keck)
startdate = zeros(size(date));
startdate(year==2010) = 20100215;
startdate(year==2011) = 20110301;
startdate(year==2012) = 20120101;
startdate(year==2013) = 20130101;
startdate(year==2014) = 20140301;
switch expt
 case('keck')
  startdate(year==2015) = 20150301;
  startdate(year==2016) = 20160321;
  startdate(year==2017) = 20170315;
  startdate(year==2018) = 20180307;
 case('bicep3')
  startdate(year==2015) = 20150419;
  startdate(year==2016) = 20160201;
  startdate(year==2017) = 20170317;
  startdate(year==2018) = 20180306; %temp
  startdate(year==2019) = 20180301; %temp
end

% likewise, also define a season end date.
enddate = Inf(size(date));
enddate(year==2010) = 20101231;
enddate(year==2011) = 20111231;
enddate(year==2012) = 20121231;
enddate(year==2013) = 20131231;
enddate(year==2014) = 20141103;
switch expt
 case('keck')
  enddate(year==2015) = 20151115;
  enddate(year==2016) = 20161130;
  enddate(year==2017) = 20171109;
  enddate(year==2018) = 20181108;
 case('bicep3')
  enddate(year==2015) = 20151231;
  enddate(year==2016) = 20161231;
  enddate(year==2017) = 20171106;
  enddate(year==2018) = 20181231; %temp
  enddate(year==2019) = 20191231; %temp
end

% find the set of matching tags
switch tagsetname

  % Be sure to exclude A/X and 00 tags from all sets

  case ['cmb' targetyearstr 'full']
    sel=year==targetyear & date>=startdate ...
      & ~ismember(phase,{'A','X'}) & setn>0 & strncmp(type,'cmb',3);
    tag=pt.tag(sel)';

  case ['cmb' targetyearstr]
    sel=year==targetyear & date>=startdate & date<=enddate ...
      & ~ismember(phase,{'A','X'}) & setn>0 & strncmp(type,'cmb',3);
    tag=pt.tag(sel)';

  case {'allcmb'}
    sel=date>=startdate & ~ismember(phase,{'A','X'}) ...
      & setn>0 & strncmp(type,'cmb',3);
    tag=pt.tag(sel)';
   
  case {'allgal'}
    sel=date>=startdate & ~ismember(phase,{'A','X'}) ...
      & setn>0 & strncmp(type,'gal',3);
    tag=pt.tag(sel)';

  case ['gal' targetyearstr]
    sel=year==targetyear & date>=startdate & date<=enddate ...
      & strncmp(type,'gal',3) & ~ismember(phase,{'A','X'}) ...
      & setn>0;
    tag=pt.tag(sel)';

  case ['all' targetyearstr]
    sel=year==targetyear & date>=startdate ...
      & (strncmp(type,'gal',3) | strncmp(type,'cmb',3)) ...
      & ~ismember(phase,{'A','X'}) & setn>0;
    tag=pt.tag(sel)';

  case {'all'}
    sel=date>=startdate ...
      & (strncmp(type,'gal',3) | strncmp(type,'cmb',3)) ...
      & ~ismember(phase,{'A','X'}) & setn>0;
    tag=pt.tag(sel)';

  % JAB 2011-03-15
  % 8 B2 standard phases (4 days) from 2010, after TES/SQ bias changes 
  case 'subset2010b2'
    % this list drove my f'ing crazy - it is not a complete
    % coverage pattern!
    %taglist={'20100918I_dk113','20100918F_dk113','20100915C_dk068','20100915E_dk068','20100914H_dk293','20100914I_dk293','20100921B_dk248','20100921C_dk248'};
    % this one is...
    taglist={'20100914H_dk293','20100914I_dk293','20100915E_dk068','20100915F_dk068','20100918H_dk113','20100918I_dk113','20100921B_dk248','20100921C_dk248'};
    
    tl={};
    for jj=1:length(taglist);
      tagtemp=taglist{jj};
      indx=strfind(tagtemp,'_');
      for ii=1:10;
	tl{ii+10*(jj-1)}=[tagtemp(1:indx-1) sprintf('%.2d',ii) tagtemp(indx:end)];
      end
    end
    tag=tl;

  case 'galsubset2010b2'
    taglist={'20100912D_dk293','20100918D_dk113','20100915D_dk068','20100921D_dk248'};
    tl={};
    for jj=1:length(taglist);
      tagtemp=taglist{jj};
      indx=strfind(tagtemp,'_');
      for ii=1:7;
	tl{ii+10*(jj-1)}=[tagtemp(1:indx-1) sprintf('%.2d',ii) tagtemp(indx:end)];
      end
    end
    tag=tl;

  case 'subsetkeck2012'
    startdate = 20120606;
    enddate = 20120613;
    sel= date>=startdate & date<=enddate ....
      & setn>0 & strncmp(type,'cmb',3);
    tag=pt.tag(sel);

  case 'subsetkeck2013'
    startdate = 20130724;
    enddate = 20130807;
    sel= date>=startdate & date<=enddate ....
      & setn>0 & strncmp(type,'cmb',3);
    tag=pt.tag(sel);

  case 'cmbmirror2015'
    startdate = 20150401;
    enddate = 20150407;
    sel= date>=startdate & date<=enddate & setn>0 & strncmp(type,'cmb',3);
    tag=pt.tag(sel);

  case 'galmirror2015'
    startdate = 20150401;
    enddate = 20150407;
    sel= date>=startdate & date<=enddate & setn>0 & strncmp(type,'gal',3);
    tag=pt.tag(sel);

  case 'coldspotmirror2019'
    startdate = 20190101;
    enddate = 20190301;
    sel= date>=startdate & date<=enddate & setn>0 & strncmp(type,'cmb',3);
    tag=pt.tag(sel);

  otherwise
    tag={};
end

for i=1:length(flag)
switch flag{i}
  % If requested strip out all tags for which tod files already exist
  case 'new'
    tags=tag(~has_tod(tag) | ~has_tod(tag,[],'_cutparams.mat'));
    
  % If requested leave out all tags that don't have tod
  case 'has_tod'
    tags=tag(has_tod(tag));
    
  case 'no_cuts'
    tags=tag(~has_tod(tag,[],'_cutparams.mat'));
    
   case 'has_cuts'
    tags=tag(has_tod(tag,[],'_cutparams.mat'));
    
  otherwise
    tags=tag; 
end
tag=tags;
end
disp(sprintf('get_tags: you have chosen to analyze %d tags',length(tags)));

return
