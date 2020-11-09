function cleanupfarmfiles(nbase_in,redo)
% cleanupfarmfiles(nbase)
% 
% nbase is first four digits of the sim serial number
%
% redo: 'wait10' do it every 10 minutes
%
% removes the .err and .out files in the directory farmfiles/nbase that lack a
% corresponding .mat file, i.e. have successfully completed
%
% ex. cleanupfarmfiles(1205)
%     cleanupfarmfiles(1205:1210)

if ~exist('redo','var')
  redo='';
end

while 1
  for k=1:numel(nbase_in)
    
    nbase=nbase_in(k);
    
    farmfiledir=sprintf('farmfiles/%04d/',nbase);
    
    % if the directory doesn't exist, exit
    if ~exist(farmfiledir,'dir')
      return
    end
    
    % These .mat files exist and thus are jobs that have not completed successfully
    % the .mat files need to be check last, since there might be more jobs
    % created during runtime.
    derr=dir([farmfiledir '/*.err']);
    dout=dir([farmfiledir '/*.out']);
    dmat=dir([farmfiledir '/*.mat']);
      
    % Turn into a cell array and shave off the filename extension
    fmat=dir2cellstr(dmat);
    ferr=dir2cellstr(derr);
    fout=dir2cellstr(dout);
    
    % These .err, .out and .mat that lack a partner
    completed_err=setxor(ferr,fmat);
    completed_out=setxor(fout,fmat);
    
    % These .err and .out files that lack a .mat file
    completed_err=intersect(ferr,completed_err);
    completed_out=intersect(fout,completed_out);
    
    % Remove files
    remove_files(farmfiledir,completed_err,'.err');
    remove_files(farmfiledir,completed_out,'.out');
    
  end
  
  if strfind(redo,'wait')
    disp(['cleaning ',num2str(nbase_in),' every ',redo(5:end),' minutes'])
    pause(60*str2num(redo(5:end)));
  else
    break
  end
  
end
  
  
return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function f=dir2cellstr(d)

for i=1:numel(d)
  f{i}=d(i).name(1:end-4);
end

if ~exist('f')
  f=[];
end

return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function remove_files(farmfiledir,files,fext)

% There are often too many files to rm them all in one command, so make a file with all
% their filenames
fname=sprintf('%s/filestoremove',farmfiledir);
fid=fopen(fname,'w');
for i=1:numel(files)
  fprintf(fid,'%s/%s%s\n',farmfiledir,files{i},fext);
end
fclose(fid);
cmd=sprintf('xargs rm < %s; rm -f %s',fname,fname);
system(cmd);

return
