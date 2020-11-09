function d=strip_nonsense(d)
%   d=strip_nonsense(d)
%
% remove the non printing nonsense from a string
% specifically used when reading strings using
% system('ls -1 something_??.mat');
%
% the non-printing is there due to color ls which may be aliased as
% default option to ls. Can therefore get rid of it by simply doing
% system('/bin/ls -1 something_??.mat');

  d=strrep(d,'[27]','');
  d=strrep(d,'[0m','');
  d=strrep(d,'[00;','');    
  d=strrep(d,'35m','');
  d=strrep(d,'[00m','');
  d=strrep(d,'[m','');
  d=strrep(d,'^[0C','');
  d=strrep(d,'[01;34m','');
  d=strrep(d,'32m','');
  d=strrep(d,'36m','');
  d=strrep(d,'34m','');
  d=strrep(d,'@','');
  d=strrep(d,'\033[','');

return