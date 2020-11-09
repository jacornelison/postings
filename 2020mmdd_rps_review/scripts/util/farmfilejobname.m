function [farmfile,jobname]=farmfilejobname(sernum,varargin)
% [jobfile,jobname]=farmfilejobname(sernum,varargin)
%
% Create babysitjobs-compatible farm files and job names using the typical
% combination of serial number and an identifier string.
%
% INPUTS
%   sernum      Serial number of the form 'NNNNXXXY' or 'NNNNreal' or the
%               generic descriptor '*BICEP', where '*' can be any arbitrary
%               string which is prepended to the farmfile path (useful for
%               putting farmfiles in a subdirectory).
%
%   varargin    A sequence of strings which will be appended to the farm file
%               and job name, separated by underscores.
%
% OUTPUTS
%   farmfile    Path name to a farm file formatted according to
%
%                   farmfiles/NNNN/XXXY_ARG1_ARG2_...
%                   farmfiles/ARG1_ARG2_...
%
%               for sernum a serial number or 'BICEP', respectively.
%
%   jobname     Corresponding job name which is compatible with the farmfile
%               for use with babysitjobs. Job names will be, respectively,
%
%                   NNNN_XXXY_ARG1_ARG2_...
%                   BICEP_ARG1_ARG2_..
%
% SEE ALSO
%   farmfile2jobname
%   jobname2farmfile

  % Simply concatenate all arguments together as strings, separated by _
  fullid = sprintf('_%s', varargin{:});

  if length(sernum)>=5 && strcmp(upper(sernum(end-4:end)),'BICEP')
    if length(sernum) > 5
      prefix = sernum(1:(length(sernum)-5));
      if prefix(end) ~= filesep()
        prefix = [prefix filesep()];
      end
      pathid = fullid(2:end);
    else
      prefix = '';
      pathid = fullid(2:end);
    end
    farmfile = sprintf('farmfiles/%s%s.mat', prefix, pathid);
    jobname = sprintf('%sBICEP%s', prefix, fullid);

  else
    % Split the serial number into the user portion and the sim portion
    simnum = sernum(5:8);
    sernum = sernum(1:4);

    is4num = @(s) strcmp(s, sprintf('%04d', str2num(s)));
    if ~is4num(sernum) || ( ~strcmp(simnum, 'real') && ~is4num(simnum) ) 
      error('farmfilejobname: Invalid sernum.')
    end

    % Then build the outputs
    farmfile = sprintf('farmfiles/%s/%s%s.mat', sernum, simnum, fullid);
    jobname  = sprintf('%s_%s%s', sernum, simnum, fullid);
  end
end
