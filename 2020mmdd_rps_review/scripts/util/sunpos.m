function [ra,dec,sunelon]=sunpos(mjd)
% [ra,dec,sunelon]=sunpos(mjd)
%
% Calculate the ra and dec of the sun for a given date.
% Also calculate the solar elongation.
% 

% convert from mjd to jd
jd=mjd+2400000.5;

% call the mex function
[ra,dec,sunelon]=sunposc(jd);

return
