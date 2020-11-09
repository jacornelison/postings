function libsphere()
% libsphere()
%
% Loads the C-library 'libsphere' into Matlab.
%

  libname = 'libsphere';

  if ~libisloaded(libname)
    here = fileparts(mfilename('fullpath'));
    lib = fullfile(here, libname);
    hfile = fullfile(here, libname, [libname '.h']);
    % Generate the <libname>_proto.m file side-by-side with the library
    % by cd-ing to the source directory
    oldcwd = cd(here);
    if ~exist([libname '_proto.m'], 'file')
      % Even though Matlab finds the library here, it ends up throwing errors
      % about not finding system libraries. If we continue using absolute
      % paths, though, no problems. Don't know why, but this is still easy
      % enough to do.
      loadlibrary(lib, hfile, 'mfilename',[libname '_proto']);
    else
      loadlibrary(lib, str2func([libname '_proto']));
    end
    cd(oldcwd)
  end
end
