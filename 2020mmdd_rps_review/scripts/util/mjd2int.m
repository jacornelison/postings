%+
% NAME:
%    mjd2int
%
% PURPOSE:
%    Converts MJD to read_interval-style 64-bit date/time integer.
%
% CALLING SEQUENCE:
%    result = mjd2int(mjd)
%
% INPUTS:
%    mjd
%
% KEYWORD PARAMETERS:
%    (none)
%
% OUTPUTS:
%    MJD
%
% MODIFICATION HISTORY
%    2007-04-30 Written. (dB)
%-


function int = mjd2int(mjd)
  
  [yr,mo,dy,hh,mm,ss]=mjd2date(mjd);
  
  ss = floor(ss + 0.5);
  mm = mm * 100;
  hh = hh * 10000;
  dy = dy * 1000000;
  mo = mo * 100000000;
  yr = yr * 10000000000;

int = yr + mo + dy + hh + mm + ss;
