function compile_code(codename)
%%%
% compile_code(codename) uses the MATLAB Compiler Runtime to compile codename.
% Currently, codename can be reduc_makesim(.m) or reduc_coaddpairmaps(.m)
% when compiling code, any addpath commands in your startup.m file should be inside an
% if ~isdeployed block.  If this is not done, the compiled code should still execute,
% but there will be extra errors/warning messages at the start of execution.

warnmsg=['Please make sure the addpath commands in your startup.m file are inside an if ' ...
         '~isdeployed statement.'];
% display warning message for user to check their startup.m file for addpath commands
disp(warnmsg)

% check that codename is a string
if ~ischar(codename)
  error('codename must be a string')
end

% change codename to end with .m so that there are less options to check
codename=strread(codename,'%s');
codename=codename{1};
if ~strcmp(codename(end-1:end),'.m')
  codename=[codename '.m'];
end

% which directory is pipeline code directory?
pipedir=fileparts(which('reduc_coaddpairmaps.m'));

% change to the pipeline directory if not already there.  Save current directory to get
% back
oldpwd=pwd;

if ~isequal(pipedir,oldpwd)
  cd(pipedir)
  % set switcheddir to true if directory was changed
  switcheddir=true;
else
  % if directory was not changed, switcheddir=false;
  switcheddir=false;
end

% compile code if it is one of the acceptable pipeline codes
switch(codename)
 case 'reduc_coaddpairmaps.m'
  compile(codename)
 case 'reduc_makesim.m'
  compile(codename)
 otherwise,
  error(['"' codename '"' ' is not a currently accepted function to compile.'])
end

setpermissions(strrep(codename,'.m',''),'770')

% change back to old pwd if directory was changed
if switcheddir
  cd(oldpwd)
end

return

function compile(codename)
% subfunction to perfom the matlab compiling function
disp(['Compiling ' codename])
eval(['mcc -R -nojvm -R -nodesktop -R -nosplash -R -singleCompThread -m ' codename])

return

