function files = listfiles(ls_string,sec_pause,n_try)
  if ~exist('sec_pause','var')
    sec_pause=5;
  end
  if ~exist('n_try','var')
    n_try=0;
  end

  files = {};

  % Previously ls was used to do this, but that results in the shell doing any
  % glob expansions, putting the entire file listing into the executed command,
  % and then ls redoing stat() calls on all the files. Using find avoids that
  % by being able to parse the glob expansions on its own.

  % Find needs a path and file spec in different locations, so we have to do a
  % bit of work to separate the two.
  if ~isempty(strfind(ls_string, '/'))
    [basedir,pattern,ext] = fileparts(ls_string);

  % If no other path components exist, just the current working directory as
  % the base and then pattern match the whole passed in string.
  else
    basedir = pwd();
    pattern = ls_string;
    ext = '';
  end
  
  % Fall back to using ls if the basedir contains any glob characters. This
  % should or could be smarter, but until the need arises, use this simple
  % solution. Also note that this is still a kludge; technically we should
  % then verify that the * and ? aren't escaped, too.
  %
  % Also fallback if using bash features like {a,b} and [0-9] in the pattern.
  % They should still work in the base since bash will expand them as multiple
  % search paths for find.
  if ~isempty(strfind(basedir, '*')) || ~isempty(strfind(basedir, '?')) ...
  || (~isempty(strfind(pattern,'[')) && ~isempty(strfind(pattern,']'))) ...
  || (~isempty(strfind(pattern,'{')) && ~isempty(strfind(pattern,'}')))
    files = lsfallback(ls_string, sec_pause, n_try);
    return
  end

  cmd = sprintf('/usr/bin/find -H %s -maxdepth 1 -iname ''%s'' -print', ...
    basedir, [pattern ext]);
  try
    [ret,msg] = system_safe(cmd);
  catch
    % Also fallback if find fails, which may happen if command line arguments
    % for ls are present in ls_string.
    warning('Call to find failed. Falling back to using ls...')
    files = lsfallback(ls_string, sec_pause, n_try);
    return
  end

  files = strread(msg, '%s', 'delimiter',sprintf('\n'));
end

function files=lsfallback(ls_string,sec_pause,n_try)
  files={};
  n=0;
  while 1
    % system_safe now throws an exception if s~=0, but it still uses the
    % command output as error message.
    try
      [s,d]=system_safe(['/bin/ls -1 ',ls_string]);
    catch ex
      s = 1;
      d = ex.message;
    end

    if ~s
      d=strip_nonsense(d);
      files=strread(d,'%s', 'delimiter', sprintf('\n'));
      break
    else
      if ~isempty(strfind(d,'No match'))
        break
      elseif ~isempty(strfind(d,'No such file or directory'))
        allFiles=strread(d,'%s', 'delimiter', sprintf('\n'));
        emptyCells = strfind(allFiles,'No such file or directory');
        emptyCells = cellfun(@isempty,emptyCells);
        files = allFiles(emptyCells);
        break
      else
        n=n+1;
        if (n_try>0 & n>n_try)
          error(['directory listing ',ls_string,' failed ',num2str(n_try),' times'])
        end
        display(['directory listing ',ls_string,' failed, try again'])
        pause(sec_pause)
      end
    end
  end
end
