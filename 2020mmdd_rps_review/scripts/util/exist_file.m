% exist_file is the same as exist(...,'file')
% except that it doesn't look on the Matlab
% search path.  Only full paths and paths
% relative to the current working directory
% will be found.  This is faster than using
% exist() and avoids rehashing the Matlab
% search path.
function s=exist_file(fname)

s=(1==length(dir(fname)));

return

