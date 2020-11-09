function parent=explode_setfield(parent,name,varargin)
% function parent=explode_setfield(parent,name,varargin)
%
% If value is a non-scalar input, then each element in value is added to
% parent as a unique field like [name '_0001', name '_0002', ...].
%
% The motivating use case is for very large data structures which are saved
% to disk. A subset of elements in an array cannot be loaded in using load()
% since it can only read an entire top-level variable.
%
%
% INPUTS
%   parent    The parent structure that exploded elements will be appended
%             to.
%
%   name      Field which should be set. The syntax can be as complicated as
%
%                 field(?).sub.member
%
%             The '(?)' allows for an array of inputs to be exploded, and can
%             occur at the end of any valid structure member name. This is
%             useful when a structure array member needs to be exploded. In
%             this syntax, 'member' is the large value to be exploded to the
%             top level, but when the structure is imploded, the 'field' array
%             will be reconstructed as well. Only one such array qualifier is
%             permitted.
%
%             The '.' allows values to be nested into structs upon implosion.
%             Only Matlab limits to nesting level apply.
%
%   varargin  Value(s) to be inserted into parent.
%
% OUTPUTS
%   parent    Modified to include newly exploded fields as well as the
%             'explodeinfo' field to record the original structure.
%
%
% EXAMPLES
%
%   % Collect variables to save into struct
%   sav.metadata = metadata;
%   sav.auxinfo = auxinfo;
%   % Explode big data field
%   sav = explode_setfield(sav, 'parts', parts);
%
%   save(filename, '-v7.3', '-struct', 'sav')
%
%
% SEE ALSO
%   implode_struct
%

  nvarargin = length(varargin);

  % explode_struct can only operate on a scalar struct which contains array
  % structs/cells. Make sure this is the case.
  if ~isstruct(parent)
    if isempty(parent)
      parent = struct();
    else
      error('explode_setfield: parent must be a structure');
    end
  end
  if ~ischar(name)
    error('explode_setfield: name must be a string');
  end

  if (numel(parent) ~= 1)
    error('explode_setfield: parent must be a scalar struct');
  end

  if ~isfield(parent,'explodeinfo')
    parent.explodeinfo.name = {};
    parent.explodeinfo.size = {};
  end

  % Multiple values assumes that one of the members itself is an array which
  % should be filled (probably a struct array with a large member variable that
  % needs to be exploded).
  arrayloc = strfind(name, '(?)');
  if nvarargin > 1 || ~isempty(arrayloc)
    % Make sure we have the info we need to explode a member array.
    if isempty(arrayloc)
      error('explode_setfield: setting struct arrays requires a ''(?)'' name qualifier.');
    elseif length(arrayloc) ~= 1
      error('explode_setfield: only one struct array nesting is supported');
    end

    % For each array member, just recursively explode with an appropriately
    % constructed name.
    for ii=1:nvarargin
      modname = strrep(name, '(?)', sprintf('_%05d', ii));
      parent = explode_setfield(parent,modname,varargin{ii});
    end

    % Then register that this field has been "exploded" as well, and therefore
    % implode_struct will do the right thing automatically (since it processes
    % the expandinfo list in order).
    parent.explodeinfo.name{end+1} = name(1:arrayloc(1)-1);
    parent.explodeinfo.size{end+1} = [1 nvarargin];

  % For the simple (or base) case, explode the input into parent structure
  % members.
  else
    % Generate a unique name for each new element being exploded
    fieldname = strrep(name,'.','_');

    newnames = arrayfun(@(i) sprintf('%s_%05d', fieldname, i), ...
        1:numel(varargin{1}), 'UniformOutput', false);

    % Save this info into the explodinfo structure
    parent.explodeinfo.name{end+1} = name;
    parent.explodeinfo.size{end+1} = size(varargin{1});

    % Start appending children
    for ii=numel(newnames):-1:1
      parent.(newnames{ii}) = varargin{1}(ii);
    end
  end

end
