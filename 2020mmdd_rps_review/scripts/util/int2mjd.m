%+
% NAME:
%    int2mjd
%
% PURPOSE:
%   Converts read_interval-style 64-bit date/time integer to MJD.
%
% CALLING SEQUENCE:
%    result = int2mjd(in)
%
% INPUTS:
%     int - 64-bit date/time integer, e.g. 20060401000000
%
%
% OUTPUTS:
%   MJD corresponding to the input date/time integer
%
% MODIFICATION HISTORY
%    2007-04-27 Written. (dB)

function mjd = int2mjd(in)
  in=double(in);
  yr = fix(mod(in / 10000000000,10000));
  mo = fix(mod(in / 100000000, 100));
  dy = fix(mod(in / 1000000, 100));
  hh = fix(mod(in / 10000, 100));
  mm = fix(mod(in / 100, 100));
  ss = fix(mod(in,100));

  mjd = date2mjd(yr,mo,dy, hh, mm, ss);
  
  
  
