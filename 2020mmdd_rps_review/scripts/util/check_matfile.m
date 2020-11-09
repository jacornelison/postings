function status=check_matfile(filename)
% status=check_matfile(filename)
%
% status:  1=ok, 0=corrupt
%
% This function will check to see if a mat-file
% is corrupted due to incomplete saving.  On occasion
% farmed jobs that create mat-files end during the
% saving operation and the files are corrupted.  To
% avoid failures using reduc_coaddpairmaps.m check to
% make sure the files are complete.  
%
% At the moment this function just checks to see if the
% saved variables are complete.  It also assumes v5 mat-files.
%
% e.g. check_matfile('pairmaps/03030045/20100317C_dk310_filtp3_weight3_gs.mat')

% open the file
f=fopen(filename,'rb');

% skip the header
fseek(f,128,'bof');

% if ftell(f)==0, then file is empty
if ftell(f)==0
  status=0;
  return
end

varinfo=1;
while ~isempty(varinfo) & varinfo(1)~=0 
  
  % get the variable info 
  % varinfo(1)=data type, varinfo(2)=# bytes
  varinfo=fread(f,2,'int32');
  
  % skip through the variable data
  if ~isempty(varinfo)
    dat=fread(f,varinfo(2),'uint8');
  end

end

% if it's empty, you've gone through
% all the variables and they're intact
if isempty(varinfo)
  status=1;
else
  % if the variable size is 0, then
  % the file is probably corrupted
  if varinfo(1)==0
    status=0;
  end
end

fclose(f);

return
