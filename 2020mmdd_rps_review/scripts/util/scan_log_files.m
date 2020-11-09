function get_run_info(d)
% get_run_info()
%
% scan the log files looking for BICEP CMB schedules
%
% Schedules which do not start with phase A will have incorrect phase
% letter and hence won't match with Cynthia/Denis list coming from
% get_tags
% Very hard to solve this problem as log files does not contain clear
% messages as to which phase is which
%
% e.g.
% get_run_info('/data/bicepdaq/log/200*')

% get list of log files
l=dir(d);

% open output file
of=fopen('get_run_info.txt','w');

% for each log file

i=1; fp=0;
while(1)

  % fetch a line
  [fp,s,i,ex]=getline(fp,l,i);

  % if file list finished
  if(s==0)
    break;
  end
  
  % if not a schedule start line skip to next line
  if(isempty(strfind(s,'Starting schedule')))
    continue;
  end
  
  % if not a 48hr schedule schedule skip to next line
  k=strfind(s,'48hr');
  if(isempty(k))
    continue;
  end
    
  % only consider schedules with numeric version number
  if(isempty(sscanf(s(k+5:end),'%d')))
    continue;
  end
  
  % record schedule start line
  schsline=s;
  
  disp(schsline);
  
  % look for skydip to indicate start of phase B
  while(1)
    [fp,s,i,ex]=getline(fp,l,i);
    if(s==0) break; end;
    if(ex|~isempty(strfind(s,'Starting command do_skydip')))
      [sd,st]=strread(s,'%s %s',1);
      break;
    end
  end
  if(ex|s==0) continue; end;
    
  % look for skydip to indicate start of phase C
  while(1)
    [fp,s,i,ex]=getline(fp,l,i);
    if(s==0) break; end;
    if(ex|~isempty(strfind(s,'Starting command do_skydip')))
      [ed,et]=strread(s,'%s %s',1);
      writeblock(of,sd,st,ed,et,schsline,'B');
      [sd,st]=strread(s,'%s %s',1);
      break;
    end
  end
  if(ex|s==0) continue; end;
  
  % infuriatingly cmb12 schedules use do_CMBfield_set instead of
  % do_transgalacticfield_set to run phase D - this means we
  % must skip forward to the phase C do_CMBfield_set
  while(1)
    [fp,s,i,ex]=getline(fp,l,i);
    if(s==0) break; end;
    if(ex|~isempty(strfind(s,'do_CMBfield_set')))
      break;
    end
  end
  if(ex|s==0) continue; end;
  
  % look for start of phase D
  while(1)
    [fp,s,i,ex]=getline(fp,l,i);
    if(s==0) break; end;
    if(ex|~isempty(strfind(s,'do_CMBfield_set'))|...
	  ~isempty(strfind(s,'do_transgalacticfield_set')))
      [ed,et]=strread(s,'%s %s',1);
      writeblock(of,sd,st,ed,et,schsline,'C');
      [sd,st]=strread(s,'%s %s',1);
      break;
    end
  end
  if(ex|s==0) continue; end;
  
  % look for skydip to indicate start of phase E
  while(1)
    [fp,s,i,ex]=getline(fp,l,i);
    if(s==0) break; end;
    if(ex|~isempty(strfind(s,'Starting command do_skydip')))
      [ed,et]=strread(s,'%s %s',1);
      writeblock(of,sd,st,ed,et,schsline,'D');
      [sd,st]=strread(s,'%s %s',1);
      break;
    end
  end
  if(ex|s==0) continue; end;
  
  % look for skydip to indicate start of phase F
  while(1)
    [fp,s,i,ex]=getline(fp,l,i);
    if(s==0) break; end;
    if(ex|~isempty(strfind(s,'Starting command do_skydip')))
      [ed,et]=strread(s,'%s %s',1);
      writeblock(of,sd,st,ed,et,schsline,'E');
      [sd,st]=strread(s,'%s %s',1);
      break;
    end
  end
  if(ex|s==0) continue; end;
  
  % look for schedule exit to indicate end of phase F
  while(1)
    [fp,s,i,ex]=getline(fp,l,i);
    if(s==0) break; end;
    if(ex)
      [ed,et]=strread(s,'%s %s',1);
      writeblock(of,sd,st,ed,et,schsline,'F');
      break;
    end
  end
  
end

fclose(of);

return

%%%%%%%%%%%%%%%%%%%%%
function writeblock(of,sd,st,ed,et,schsline,phalet)

k=strfind(schsline,'48hr');

% get schedule start date
ssd=strread(schsline,'%s',1);

% decode schedule name
[dummy,v,galdk]=strread(schsline(k:end),'%s_%d_%s_',1,'whitespace','_');
switch galdk{1}(end)
  case 'A'
    dk=315;
  case 'B'
    dk=0;
  case 'C'
    dk=135;
  case 'D'
    dk=180;
  otherwise
    dk=999;
end

spn=datenum([sd{1},st{1}],'yymmddHH:MM:SS');
sp=datestr(spn,'dd-mmm-yyyy:HH:MM:SS');
epn=datenum([ed{1},et{1}],'yymmddHH:MM:SS');
ep=datestr(epn,'dd-mmm-yyyy:HH:MM:SS');

fprintf(of,'%s %s B20%s_cmb%02d_%s_%c_%03d %s %5.2f\n',sp,ep,ssd{1},v,galdk{1},phalet,dk,schsline(k-2:end),(epn-spn)*24);

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [fp,s,i,ex]=getline(fp,l,i)
  
% if end of file close old file
if(feof(fp))
  fclose(fp);
  fp=0;
end

% if no file open start new one
if(fp==0)
  i=i+1;
  
  % if done all files return with empty string
  if(i>length(l))
    s=0; ex=false;
    return
  end
  
  fp=fopen(['/data/bicepdaq/log/',l(i).name]);
end

% get string
s=fgetl(fp);

% flag if schedule exit detected
ex=~isempty(strfind(s,'Exiting schedule'));

return
