function saveandtest(varargin)
% saveandtest(varargin)
%
% Call just like matlab save function. 
%
% Save temporary file, then attempt to load the smallest variable. If this fails, attempt to
% resave up to 2 times. If the save does not work, remove any (corrupt) file that might
% have been generated and throw an error. If the save does work, move temporary file to
% the real file name.
%
% i.e. saveandtest('filename.mat','ac','m','mapopt','-v7');

% Construct the file name
fname=varargin{1};
[a,b,c]=fileparts(fname);
if isempty(c) & isempty(intersect(varargin,'-ASCII'))
  c='.mat';
end
fname=fullfile(a,[b c]);

% Temporary file name to use before intactness has been verified
tmpfname=fullfile(a,[b '.' gen_stamp() 'tmp' c]);

% Construct the save command
cmd=['save(''',tmpfname,''','];
for i=2:numel(varargin)
  cmd=[cmd,'''',varargin{i},''','];
end
cmd=[cmd(1:end-1),');'];

% Clear the last warning so we can unambiguously identify a warning issued
% by this invocation to save.
lastwarn('')
% Evaluate the save command in the workspace above
evalin('caller',cmd);
[warnmsg,warnid] = lastwarn();
% If this is a warning about upgrading to v7.3 format, upgrade to an error
% since we definitely want to keep the data.
if strcmp(warnid,'MATLAB:save:sizeTooBigForMATFile')
  error(warnid, warnmsg)
end

% Now load the resulting file and resave if it breaks. Attempt to load the smallest
% variable to check goodness. If it loads then the file ought to be fine. This is
% both for speed and to avoid loading two copies of very large variables.

if isempty(intersect(lower(varargin),'-struct'))
  % For a save which did not have the -struct parameter included, read the
  % sizes from memory to minimize the amount of filesystem I/O.
  vars=sprintf('''%s'',',varargin{2:end});
  vars=vars(1:end-1); % remove unwanted trailing comma after last argument
  S=evalin('caller',['whos(',vars,');']);
  b=horzcat(S(:).bytes);
else
  % If the contents of a struct where saved, then read the sizes of the
  % constituent parts directly from the file since whos can't introspect
  % the sizes of struct fields.
  S=whos('-file',tmpfname);
  b=horzcat(S(:).bytes);
end

% Smallest variable
[dum,ind]=nanmin(b);

% Now load the variable. If it doesn't work, resave.
saveok=0;
Nsaves=0;
maxNsaves=2;
while ~saveok & Nsaves<maxNsaves
  try
    lastwarn('')
    load(tmpfname,S(ind).name);
    % Trying to load a variable which doesn't exist actually only emits a
    % warning, so make this an error instead.
    [warnmsg,warnid] = lastwarn();
    if strcmp(warnid, 'MATLAB:load:variableNotFound')
      error(warnid, warnmsg)
    end
    saveok=1;
  catch
    disp(sprintf('File %s did not save properly: resaving.',fname));
    evalin('caller',cmd);
    Nsaves=Nsaves+1;
  end
end

% If the save didn't work, remove the file and throw an error.
if Nsaves==maxNsaves
  system_safe(sprintf('rm -f %s',tmpfname));
  error(sprintf('File %s did not save',fname));
end

% File is good, so move to real file name
system_safe(['mv ',tmpfname,' ',fname]);

return

