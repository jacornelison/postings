function rm_corrupt_matfiles(p,chunk,nchunks)
% rm_corrupt_matfiles(dir)
%  or
% rm_corrupt_matfiles(dir,chunk,nchunks)
%
% i.e. 
% for i=1:100
%  farmit('farmfiles/',sprintf('rm_corrupt_matfiles(''./'',%d,100)',i));
% end


% defaults
if ~exist('chunk','var')
  chunk=1;
  nchunks=1;
end
if ~exist('p','var')
  p='./';
end

% get the mat files
d=dir([p,'/*.mat']);

% break up, useful for farming
nfiles=numel(d);
chunksize=ceil(nfiles/nchunks);

% check these files
ind=(1:chunksize) + (chunk-1)*chunksize;

for i=ind
  
  if i>nfiles
    break
  end
  
  if mod(i,100)==1
    disp(sprintf('loading file %d of %d',i,nfiles));
  end
  
  fn=fullfile(p,d(i).name);

  try
    load(fn);
  catch
    cmd=sprintf('rm -f %s',fn);
    disp(cmd);
    system(cmd);
  end

end

return

