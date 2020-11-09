function pid=myprocesses(str)
% pid=myprocesses(string)
%
% Search your processes for string.  Search is case sensitive. String can be a regular
% expression. 
% Returns PIDs of processes
% 
% ex. pid=myprocesses('matlab');
% ex. pid=myprocesses(''); % lists all processes owned by you

username=whoami();

% list process ID because when calling from matlab, the pgrep process itself is listed
% and we need to remove it
searchstr=['pgrep -f -U ' username ' "' str '"'];
[dum,p]=system(searchstr);
pid=str2num(p);
pid=pid';

% if process does not exist, it was probably the pgrep process itself.  don't return it
keepind=false(size(pid));
for i=1:numel(pid)
  [dum,result]=system(sprintf('ps %d',pid(i)));
  keepind(i)=~isempty(strfind(result,num2str(pid(i))));
end

pid=pid(keepind);

return
