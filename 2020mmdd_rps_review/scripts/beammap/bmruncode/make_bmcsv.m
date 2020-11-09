function make_bmcsv(year, exp)
% Script to generate a nice archival .csv for the Keck or BICEP3 beam maps.
% year should be 2014 or 2015
% exp should be 'keck' or 'bicep3'
% The script is general enough to work for Keck/B3 experiments. 
% You must run it from a directory that has the
% correct symlinks (ie beammaps, log, arc, aux_data, plots)
% Don't overwrite a different experiment!
  
if ~strcmp(get_experiment_name(),exp);
  error('Get into the right folder for %s!',exp)
end

savefile = sprintf('beammaps/bmrunlist_%s_tmp.csv',year);

bm = bm_info(year);

for ii = 1:length(bm)
  
  p.number(ii) = bm(ii).number;
  p.t1{ii} = bm(ii).t1;
  p.t2{ii} = bm(ii).t2;
  p.dk(ii) = bm(ii).dk;
  p.sched{ii} = bm(ii).sched;
  p.filename{ii} = bm(ii).filename;
  p.notes{ii} = bm(ii).notes;
  
  % For Keck, update this when we switch the mirror position
  switch year
    case 2014
      mirrornum = 25;
    case 2015
      mirrornum = 29;
  end
  
  if ii < mirrornum
    p.mirror{ii} = 'back';
  else
    p.mirror{ii} = 'forward';
  end
  
end

k.comments{1} = [sprintf('# BEAM MAPPING RUNS: %s  %s',upper(exp),year)];
k.fields = {'number','t1','t2','dk','mirror','sched','filename','notes'};
k.formats = ...
    {'integer','string','string','double','string','string','string','string'};
k.units = ...
    {'unitless','time','time','deg','unitless','unitless','time','unitless'};

ParameterWrite(savefile,p,k);

return

function bm = bm_info(year)
  
  % Get relevant log info
  disp('### Parsing log files for fff beam map data .... ')
  % Change this for BICEP3 (mapomast?)
  [a,b] = system(sprintf('grep 7_ffflat_dslmast log/%s0[1-4]*.log',year));
  tok1 = strfind(b,'Starting');
  tok2 = strfind(b,'Exiting');
    
  % check that there is a matching Exiting for each starting.
  for i = 1: min(numel(tok1), numel(tok2))
    delta = tok2(i)-tok1(i);
    if delta > 150
      tok1(i)=[];
    end  
  end
  
  k = 0
  % Loop over found schedules
  for i = 1:length(tok1)
        
    % find the t1
    p1 = tok1(i);
    p = p1;
    hr = b(p-9:p-8);
    mn = b(p-6:p-5);
    se = b(p-3:p-2);
    yr = b(p-36:p-33);
    mo =  b(p-32:p-31);
    da =  b(p-30:p-29);
    t1str =  datestr([str2num(yr),str2num(mo),str2num(da),str2num(hr),str2num(mn),str2num(se)],...
                        'yyyy-mmm-dd:HH:MM:SS');
    
    % to list the arc files right after t1
    t1_tmp =  datestr([str2num(yr),str2num(mo),str2num(da),str2num(hr)+1,str2num(mn),str2num(se)],...
                        'yyyy-mmm-dd:HH:MM:SS');
    % find the t2
    p2 = tok2(i);
    p=p2;
    hr = b(p-9:p-8);
    mn = b(p-6:p-5);
    se = b(p-3:p-2);
    yr = b(p-36:p-33);
    mo =  b(p-32:p-31);
    da =  b(p-30:p-29);
    t2str =  datestr([str2num(yr),str2num(mo),str2num(da),str2num(hr),str2num(mn),str2num(se)],...
                        'yyyy-mmm-dd:HH:MM:SS');
    
    % check that the schedule lasts for at least 1 hour
    t2 = datenum(t2str,'yyyy-mmm-dd:HH:MM:SS');
    t1 = datenum(t1str,'yyyy-mmm-dd:HH:MM:SS');
    if (t2-t1)*24*3600 < 3600
      continue
    end
   
    k = k+1;
    bm(k).number=k ;
    bm(k).t1 = t1str;
    bm(k).t2 = t2str;
                  
    % find the sch name and dk angle
    sch = b(p1+51:p1+78);
    sched_type = sch(20);
    bm(k).dk = schedule_name_2_dk_angle(sched_type);
    bm(k).sched = sch;
    bm(k).notes = '';

    bm(k)
    
    % find the filename
    disp('### Parsing arc files to find arc filename .... ')
    files = list_arc_files('arc/',t1str,t1_tmp,'.dat.gz');
    
    % Loop over arc files:
    for fl = 1:length(files)
      file = files{fl}
      load_file = char(strcat('arc/',file));
      % Load the data
      d = load_arc(load_file);
      disp('File loaded');
      
      %Snippet of code from lines 117 to l158 of bm_demod.m to check if arc file contains valid data
      % Data valid flag:
      mark = d.array.frame.features;
      
      % Check the data mode
      dm = d.mce0.rc1.data_mode(logical(mark));  
      % Check to see that there is valid data:
      if numel(dm) == 0
        disp(['File ' load_file ' is empty'])
        continue 
      end

      % Generate field scans
      % Fast/slow ratio
      sampratio = ...
          length(d.antenna0.pmac.fast_az_pos)/length(d.antenna0.pmac.az_pos); 
      try
        % When is feature bit 1 on?
        fs = find_scans(d.antenna0.tracker.scan_off(:,1),...
                        bitand(mark,2^1),sampratio,'keepend');
      catch
        warning('No field scans found in this file.  Moving to next file.')
        continue
      end
      bm(k).filename = file(1:end-7);
      disp('Found scans');
      break
      
    end
    
  end
  
return

function dk = schedule_name_2_dk_angle(name)
  if ischar(name)
    switch name
     case 'a'
      dk = -122;
     case 'b'
      dk = -50;
     case 'c'
      dk = 22;
     case 'd'
      dk = 94;
     case 'e'
      dk = 166;
     case 'f'
      dk = -86;
     case 'g'
      dk = -14;
     case 'h'
      dk = 58;
     case 'i'
      dk = 130;
     case 'j'
      dk = -158;
    end
  end
    
return

