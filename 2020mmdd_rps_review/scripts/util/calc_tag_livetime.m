% lt=calc_tag_livetime(tags,logdir)
%
% Parse the log files to calculate the time
% spent in various activities.  The tag(s)
% may be a single scanset or a cell array of
% tag names, and may include entire phases
% or schedules: for example,
%   calc_tag_livetime('20100605D')
% will sum up the times from all scansets
% within this D phase.
%
% logdir = directory with GCP log files
% (defaults to 'log')
%
% output:
%   lt(:,1) = total time
%   lt(:,2) = time on source
%   lt(:,3) = in el nods
%   lt(:,4) = in partial load curves
%   lt(:,5) = in fridge cycles
%   lt(:,6) = sky dip
%   lt(:,7) = full load curve

% RWO 20100609

function P=calc_tag_livetime(tags,logdir)

if ischar(tags)
  tags={tags};
end
if nargin<2 || isempty(logdir)
  logdir='log';
end

% Identify date ranges of all log files
if length(tags)>1
  disp(['Scanning log files in ' logdir]);
end
d=dir(fullfile(logdir,'20*_*.log'));
logfiles=[];
idx=0;
for j=1:length(d)
  dn=datenum(d(j).name(1:15),'yyyymmdd_HHMMSS');
  idx=idx+1;
  logfiles(idx).name=d(j).name;
  logfiles(idx).tstart=dn;
end

% Loop over tags
for j=1:length(tags)

  if length(tags)>1
    disp(['  ' tags{j}]);
  end

  % Expand any partial tag names
  % and sum up live times if needed
  tinfo=get_run_info(tags{j});
  subtags=tinfo.tag;
  tmp_P=[];
  for ontag=1:length(subtags)

    % Get tag start and end
    tinfo=get_run_info(subtags{ontag});
    t1=datenum(tinfo.tstart,'dd-mmm-yyyy:HH:MM:SS');
    t2=datenum(tinfo.tend,'dd-mmm-yyyy:HH:MM:SS');

    % Identify log files of interest
    k1=find([logfiles(:).tstart]<t1);
    k2=find([logfiles(:).tstart]<t2);
    k1=k1(end);
    k2=k2(end);

    % Loop over log files
    for k=1:30
      fbits(k).val=0;
      fbits(k).ton=[];
      fbits(k).toff=[];
    end
  
    tt=[];
    for k=k1:k2
      f=fopen(fullfile(logdir,logfiles(k).name),'rt');
      while ~feof(f)
        ll = fgetl (f);
        if ~ischar(ll) || isempty(ll)
          continue;
        end;
        llcell = textscan (ll, '%s %s %[^\n]');
        if (length(llcell) < 3) || isempty(llcell{3})
          continue;
        end;
        tmp_date = llcell{1}{1};
        tmp_time = llcell{2}{1};
        ll = strtrim(llcell{3}{1});
        tt=[];
  
        [fplus,fminus]=parse_mark_cmd(ll);
        if ~isempty(fplus) || ~isempty(fminus)
          tt=datenum([tmp_date ' ' tmp_time],'yymmdd HH:MM:SS');
          if tt<t1
            continue
          end
          if tt>t2
            break
          end
        end

        if ~isempty(fplus)
          for b=fplus
            if (b+1)>length(fbits)
              continue
            end
            if fbits(b+1).val==0
              fbits(b+1).val=1;
              fbits(b+1).ton=[fbits(b+1).ton tt];
            end
          end
        end
        if ~isempty(fminus)
          for b=fminus
            if (b+1)>length(fbits)
              continue
            end
            if fbits(b+1).val==1
              fbits(b+1).val=0;
              fbits(b+1).toff=[fbits(b+1).toff tt];
            end
          end
        end
  
      end
      fclose(f);
    end
  
    for k=1:length(fbits)
      if length(fbits(k).toff)<length(fbits(k).ton)
        fbits(k).toff(end+1:length(fbits(k).ton))=t2;
      end
    end

    tmp_P(ontag,1)=t2-t1;                                 % Total time
    tmp_P(ontag,2)=sum(fbits(1+1).toff-fbits(1+1).ton);   % On source
    tmp_P(ontag,3)=sum(fbits(3+1).toff-fbits(3+1).ton);   % El nod
    tmp_P(ontag,4)=sum(fbits(14+1).toff-fbits(14+1).ton); % Partial load curve
    tmp_P(ontag,5)=sum(fbits(2+1).toff-fbits(2+1).ton);   % Fridge cycle
    tmp_P(ontag,6)=sum(fbits(7+1).toff-fbits(7+1).ton);   % Sky dip
    tmp_P(ontag,7)=sum(fbits(6+1).toff-fbits(6+1).ton);   % Full load curve

if any(tmp_P(:)<0)
keyboard
end

  end

  P(j,:)=nansum(tmp_P,1);
end

% Note: for fridge cycle (feature bit 2), ignore
%       actual mark commands.  These get screwed
%       up by star pointing.  Instead, use sched
%       boundaries, and fake the f2 feature bits.
function [fplus,fminus]=parse_mark_cmd(ll)
  fplus=[];
  fminus=[];
  time=NaN;
  sss = 'mark add';
  if strncmp (ll, sss, length(sss))
    [tmp,tpl]=strtok(ll,',');
    tpl=strtrim(tpl);
    if tpl(1)==','
      tpl(1)='';
    end
    tpl=strtok(tpl,',');
    tpl(tpl=='f' | tpl=='+') = ' ';
    fplus=str2num(tpl);
    fplus(fplus==2) = [];
    return
  end
  sss = 'mark remove';
  if strncmp (ll, sss, length(sss))
    [tmp,tpl]=strtok(ll,',');
    tpl=strtrim(tpl);
    if tpl(1)==','
      tpl(1)='';
    end
    tpl=strtrim(tpl);
    if strcmp(tpl,'all')
      fminus=[0:30];
      fminus(fminus==2) = [];
      return
    end
    tpl=strtok(tpl,',');
    tpl(tpl=='f' | tpl=='+') = ' ';
    fminus=str2num(tpl);
    fminus(fminus==2) = [];
    return
  end
  s=check_for_cycle_boundary(ll);
  if s==1
    fplus=2;
    return
  elseif s==-1
    fminus=2;
    return
  end
  return

function s = check_for_cycle_boundary (ll)
  s = 0;
  sss = 'cycleFridge';
  if strncmp (ll, sss, length(sss))
    tmp = ll((1+length(sss)):end);
    tmp = strtrim (tmp);
    if isempty(tmp)
      s = 1;
    elseif strcmp (tmp, 'state=0')
      s = -1;
    end;
  end;

