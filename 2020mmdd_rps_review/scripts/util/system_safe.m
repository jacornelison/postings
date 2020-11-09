% system_safe(...) is the same as system() or unix(),
% but it makes sure that user input doesn't get mixed in with
% the output from the shell command.  For a discussion of the
% problem, see:
% http://www.mathworks.com/matlabcentral/newsreader/view_thread/257052
function [s,b]=system_safe(cmd,varargin)

% arbitrary magic string lets us know when echo is definitely off
% and the real output will begin
STARTSTR='xyz1xyz2xyz3';
% turn echo off, then print magic string, then run real command
[s,btmp]=system(['stty -echo ; /bin/echo "' STARTSTR '" ; ' cmd],varargin{:});
if s
  nind = strfind(btmp,sprintf('\n'));
  nind = nind(1)+1;
  error(btmp(nind:end))
end
% for output, grab everything after the magic string
while ~strncmp(btmp,STARTSTR,length(STARTSTR))
  btmp(1)='';
end
btmp(1:length(STARTSTR))='';
btmp(1)='';  % for the newline
b=btmp;

return

