function parent=implode_struct(parent)
% function parent=implode_struct(parent)
%
% If parent has fields which were created by explode_setfield, the elements
% will be imploded back.
%
% The motivating use case is for very large data structures which are saved
% to disk. A subset of elements in an array cannot be loaded in using load()
% since it can only read an entire top-level variable. Using setfield_expand,
% though, a large array can be saved as many variables (using the '-struct'
% option in save()), the desired entries be explicitly loaded (much more
% quickly), and then finally put back into a canonical form with
% implode_struct().
%
%
% NOTES
%   Any missing expected entries are silently skipped, and the field is
%   automatically expanded to the ideal size. For example, if there should an
%   a field expanded to 100 fields (parent.field is [1x100]) but only
%   parent.field_00001 and parent.field_00010 are present, those to entries
%   will be placed into parent.field(1) and parent.field(10), and then
%   parent.field will be extended to have size 1 x 100; the unset entries will
%   all use Matlab's default fill values.
%
% INPUTS
%   parent    The parent structure that needs to have fields imploded. This
%             *must* contain the an 'explodeinfo' field, which will usually
%             need to be specified in the corresponding load() call.
%
%
% OUTPUTS
%   parent    Modified structure with no exploded fields.
%
%
% EXAMPLES
%
%   % Enumerate parts(22:32) names
%   parts = cellfun(@(n) sprintf('part_%05d',n), 22:32, 'UniformOutput',false);
%   % Load relevant data
%   data = load(filename, 'metadata', 'auxinfo', 'explodeinfo', parts{:});
%   % Implode structure
%   data = implode_struct(data);
%
%
% SEE ALSO
%   explode_setfield
%

  if (~isstruct(parent))
      error('implode_struct: parent must be a structure');
  end
  
  if (numel(parent) ~= 1)
    error('implode_struct: parent must be a scalar struct');
  end

  % If there is no explodeinfo member field, just pass through with no
  % modifications
  if (~isfield(parent,'explodeinfo'))
    return
  end

  % Go through each entry in the explodeinfo registry
  for ii=1:length(parent.explodeinfo.name)
    fullname = parent.explodeinfo.name{ii};
    membname = strrep(fullname, '.', '_');
    numelems = prod(parent.explodeinfo.size{ii});

    % Figure out the deepest-nested name of the field we need to fill
    fieldparts = regexp(fullname, '\.', 'split');

    % Enumerate the mangled names to be removed from parent and put into the
    % new temporary structure.
    manglenames = arrayfun(@(i) sprintf('%s_%05d', membname, i), ...
        1:numelems, 'UniformOutput', false);
    remfields = {};

    % Now copy in the relevant data
    working = struct();
    fname = fieldparts{end};
    for jj=1:numelems
      if isfield(parent, manglenames{jj})
        working.(fname)(jj) = parent.(manglenames{jj});
        remfields{end+1} = manglenames{jj};
      end
    end
    % Free up a bunch of memory while also cleaning up parent
    parent = rmfield(parent, remfields);

    % If the field does not exist in working, then this iteration became a
    % noop and we don't have to continue.
    if ~isfield(working, fname)
      continue
    end

    % Fix-up the size of the field if only a subset actually existed
    if length(working.(fname)) ~= numelems
      % Copy the first to one past the end to have Matlab fill all intermediate
      % entries without us having to worry about the actual form of each
      % element.
      working.(fname)(numelems+1) = working.(fname)(1);
      % Now truncate to the proper length
      working.(fname) = working.(fname)(1:end-1);
    end

    % Make sure we maintain the proper shape
    working.(fname) = reshape(working.(fname), parent.explodeinfo.size{ii});

    % Now nest the working structure as necessary. We've already consumed the
    % last one.
    for kk=(numel(fieldparts)-1):-1:1
      newworking.(fieldparts{kk}) = working;
      working = newworking;
    end

    % If the destination field already exists in the parent, then we can't
    % clobber it by just setting the value. Instead, merge the struct by setting
    % new fields as necessary.
    if isfield(parent, fieldparts{1})
      subfields = fieldnames(working.(fieldparts{1}));
      for jj=1:length(subfields)
        parent.(fieldparts{1}).(subfields{jj}) = working.(fieldparts{1}).(subfields{jj});
      end

    % Otherwise if it doesn't exist, no merging of fields is necessary
    else
      % Since we've consumed the top level name in both the nested and unnested
      % cases, simply move the data directly into parent.
      parent.(fieldparts{1}) = working.(fieldparts{1});
    end
  end

  % Once everything has been imploded, we can safely remove the explodeinfo
  % registry.
  parent = rmfield(parent,'explodeinfo');
end

