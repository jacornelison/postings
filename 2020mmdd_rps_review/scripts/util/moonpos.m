function [ra,dec]=moonpos(mjd)
% [ra,dec]=moonpos(mjd)
%
% Calculate the ra and dec of the moon for a given date
% 

% convert from mjd to jd
jd=mjd+2400000.5;

% call the mex function
[ra,dec]=moonposc(jd);

return
