function [apsfiles,inmaps]=get_apsfilename(apsopt,mapname,wildcard)
% [apsfiles,inmaps]=get_apsfilename(apsopt,mapname,wildcard)
%
% Generates the list of all input maps and resultant APS files for the given
% APS options and map name (wildcards permitted) + replacement patterns.
%
% INPUTS
%   apsopt    An apsopt structure, as passed to reduc_makeaps().
%
%   mapname   Cell array of initial map name and replacments; see
%             reduc_makeaps() for full description.
%
%   wildcard  Defaults to true. If false, no wildcard expansion is performed
%             (via listing files on disk), in which case the maps will not
%             have to necessarily exist.
%
% OUTPUTS
%   apsfiles  Cell array of all output APS file names (or a non-cell character
%             array if only a single output is generated).
%
%   inmaps    2D cell array of input map filenames after wildcard expansion
%             and pattern substitution. Dimension 1 runs through independent
%             APS inputs, and dimension 2 is each explicit map to be read in
%             for a single APS.
%
% EXAMPLE
%   apsopt = get_default_apsopt();
%   mapname = {'1459/???5_aabde_filtp3_weight3_gs_dp1100_jack0.mat', ...
%              {'5_','6_'}, {'1459','1615', '5_', '6_'}};
%   [apsfiles,inmaps] = get_apsfilename(apsopt, mapname);
%

if ~exist('wildcard', 'var')
  wildcard = true;
end

% guarantee all relevant fields are defined
apsopt = get_default_apsopt(apsopt);

if ~iscell(mapname)
  mapname = {mapname};
end

% this makes it backward compatible. If second part in mapname is not
% a cell, make it a cell of total replacements; this indicates that just
% two (or more) maps are being read.
if length(mapname)>1 && ~iscell(mapname{2})
  mapname = {mapname{1}, {mapname{2:end}}};
end

if ~strncmp(mapname{1}, 'maps/', 5)
  mapname{1} = ['maps/' mapname{1}];
end

% this is the version which concatenates maps along the first dimension
% and uses the string replacement pattern
% it also works if just one map is handed in
% first fetch the maps according to the wildcard:
if wildcard
  maps1 = sort(listfiles(mapname{1}));
else
  maps1 = mapname(1);
end
nmaps = length(maps1);

apsfiles = cell(nmaps, 1);
inmaps = cell(nmaps, length(mapname));
inmaps(:,1) = maps1;
for ii=1:nmaps
  % do string replacement to create additional filenames:
  % each cell array in mapname is one more filename.
  % in the cell array has the string replacement pattern:
  % the first string is replaced with the second string,
  % the third with the forth and so on.
  % {'B2map*',{'B2','Keck','jack0','jack01'},{'B2','B1'},...}.
  % loop over the cells in mapname
  for jj=2:length(mapname)
    % the initial map name is going to be tweaked to produce the other map names:
    map_next = inmaps{ii,1};
    % loop through the replacement in key/values manner
    for kk = 1:2:length(mapname{jj})-1
      map_next = regexprep(map_next, mapname{jj}{kk}, mapname{jj}{kk+1}, 'once');
    end
    inmaps{ii,jj} = map_next;
  end

  apsfiles{ii} = combine_mapnames(apsopt, inmaps(ii,:));
end

if length(apsfiles) == 1
  apsfiles = apsfiles{1};
end

end

function apsfile=combine_mapnames(apsopt,mapname)

%assemble the filename by summing the maps that go into it:
apsfile = '';

if iscell(mapname)

  % We have two maps to do either a cross spectrum or the experiment jack

  for ii=1:length(mapname)
    cmapname=mapname{ii};
    [p,f]=fileparts(mapname{ii});
    mapdirs{ii} = strrep(p,'maps/','');
    apsfile = [apsfile, '_', f];
  end

  %remove leading '_'
  apsfile = apsfile(2:end);

  % append .mat
  apsfile = [apsfile, '.mat'];

  % use unique basenumber to create aps subdir, with x inbetween
  % the sorting will be in the order of appearance:
  [dum,im] = unique(mapdirs);
  mapdirs=mapdirs(sort(im));

  apsdir = mapdirs{1};
  for ii=2:length(mapdirs)
    apsdir = [apsdir 'x' mapdirs{ii}];
  end
  apsfile = ['aps/' apsdir '/' apsfile];

else
  
  % Just a single map
  apsfile=strrep(mapname,'maps/','aps/');
  
end

% if freq jack we mod output filename
jackn=strfind(apsfile,'jack')+4;
switch(apsopt.howtojack)
  case 'dim1'
   if(apsfile(jackn)~='0')
     error('fed me a non jack0 map to freq diff');
   end
   apsfile(jackn)='f';
end

% append "pureB" if using pure-B estimator
if(strcmp(apsopt.pure_b,'kendrick'))
  apsfile=strrep(apsfile,'.mat','_pureB.mat');
end

% append "matrix" if using matrix proj
if isfield(apsopt,'purifmat')
  apsfile=strrep(apsfile,'.mat','_matrix.mat');
else
  if isfield(apsopt,'purifmatname') && ~isempty(apsopt.purifmatname)
    apsfile=strrep(apsfile,'.mat','_matrix.mat');
  end
end

% common mask filename extension
if apsopt.commonmask
  apsfile=strrep(apsfile,'.mat','_cm.mat');
end

% append overrx if coadding over rx
if apsopt.coaddrx==1
  apsfile=strrep(apsfile,'.mat','_overrx.mat');
end

% append overfreq if coadding over freq
if apsopt.coaddrx==2
  apsfile=strrep(apsfile,'.mat','_overfreq.mat');
end

% append overall if coadding overall
if (iscell(apsopt.overall) || apsopt.overall) && ~isempty(apsopt.overall)
  apsfile=strrep(apsfile,'.mat','_overall.mat');
end

% append a daughter aps set if requested
if(isfield(apsopt,'daughter'))
  apsfile=strrep(apsfile,'.mat',['_' apsopt.daughter '.mat']);
end

%append bintype if it is not bicep_norm
if(~strcmp(apsopt.bintype,'bicep_norm'))
  apsfile=strrep(apsfile,'.mat',['_' apsopt.bintype '.mat']);
end

% append noisemask to filename
if(isfield(apsopt,'noisemask'))
  apsfile=strrep(apsfile,'.mat','_t.mat');
end

end
