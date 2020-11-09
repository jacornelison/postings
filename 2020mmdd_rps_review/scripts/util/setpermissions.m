function setpermissions(fname,perm)
% functionsetpermissions(filename,perm)
%
% Set the permissions of a file.  The filename may contain
% wildcards.  By default the permissions are 664, that is:
% readable by everyone, and writeable by the owner and members
% of the group.

if nargin<2 || isempty(perm)
  perm='664';
end

if exist(fname,'file')
  [stat,res]=system(['chmod ' num2str(perm) ' ' fname]);
else
  [fdir fbase fext]=fileparts(fname);
  if isempty(fdir)
    fdir='./';
  end
  d=dir(fname);
  if isempty(d)
    d=dir([fname '.mat']);
  end
  for i=1:length(d)
    [stat,res]=system(['chmod ' num2str(perm) ' ' fullfile(fdir,d(i).name)]);
  end
end

