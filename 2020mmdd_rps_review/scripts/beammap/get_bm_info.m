function p = get_bm_info(year,filename)
% p = get_bm_info(year)
%
% Gets information on Pole far field beam maps 
%
% INPUTS
%       year:     looks for a .csv file based on the year
%                 Currently: beammaps/bmrunlist_201?.csv
%       filename: Optional - filename if don't want default
%
% OUTPUTS
%       p:     struct containing all the fields in the .csv
%

if exist('filename','var')
  beamfile = filename;
else
  beamfile = ['beammaps/bmrunlist_' num2str(year) '.csv'];
end

[p,k] = ParameterRead(beamfile);

return