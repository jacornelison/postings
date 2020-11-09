function [p,k]=ParameterRead(filename)
% [p,k]=ParameterRead(filename)
%
% Matlab version of Ganga's IDL ParameterRead.pro
%                   Kiwon's IDL aux_data_read.pro
%
% e.g.: p=ParameterRead('fp_pars_20050108_with_notes.csv')
%
%  p is parameters, which are stored in the tabular body of the file
%  k is key, containing info from the header and table column headings
%  Modification:
%  May 6 2007: complete rewrite by dB to improve speed ( get rid of
%  double read) and to avoid problems when reading in  files with
%  multiple ,,, in a row.
% Apr 28 2008, fixed the  error when multiple occurance try ot open same
% file (dB)
  
fid = fopen(filename);
if(fid==-1)
  % if file exists, then keep trying
  if exist(filename,'file') == 2 
    while fid == -1
      [fid,errmsg] = fopen(filename);
      disp(['Error message "' errmsg '" opening file ' filename '.  Will try again.']);
      pause(0.5);
    end
  else
     error(sprintf('Cannot open file because file %s does not exist',filename))
  end
end
  

k=[];
k.comments = {};

while(1)
  line=fgetl(fid);
  if(strncmp(line,'BODY',4))
    break;
  end;
  if any(line(1:2)=='#')   % record comment lines that begin with '#' or ' #'
    k.comments = {k.comments{:} line};
  else                     % record header lines that contain a keyword and a value
    [key,val]=strtok(line,','); val=strtrim(val(2:end));
    v=str2double(val);
    if(~isnan(v))
      val=v;
    end
    k=setfield(k,lower(key),val);
  end
end

fields  =splitstr(lower(fgetl(fid)));   % force field names to lowercase
units  =splitstr(fgetl(fid));
formats=splitstr(fgetl(fid));

% JMK change, 10/06
 % if these lines start with extra cells "Units" and "Format" (an historical problem for QUAD fp_pars files), remove them
 if length(units)==(length(fields)+1)
   units = units(2:end);
 end
 if length(formats)==(length(fields)+1)
   formats = formats(2:end);
 end

 k.fields  = fields;
 k.units   = units;
 k.formats = formats;
 
%;create formatstring to read in the BODY of the data
 formatstring='';
 for i= 1:length(formats)
   if  strcmp(char(formats{i}),'string')
     formatstring=strcat(formatstring,' %s');
   else
     formatstring=strcat(formatstring,' %n');
   end
 end
 
data=textscan(fid, formatstring ,'delimiter', ',', 'commentStyle','#');
fclose(fid);
 
p=cell2struct(data, fields, 2);

% Put data in structure
% for j=1:length(data)
%   eval(['p.' fields{j} ' = data{j};']);
% end
 

 return

 %%%%%%%%%%%%%%%%%%%%%
 function spl=splitstr(s)

 spl={};
 i=1;
 while(length(s)>1)
   [sp,s]=strtok(s,',');
   spl{i}=strtrim(sp);
    i=i+1;
end

return

