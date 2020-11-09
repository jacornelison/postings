function B=zerodup(A)
% function B=zerodup(A)
%
% Create a duplicate of A which has zeroed all elements.
%
% INPUTS
%   A    A structure, cell array, numeric array, string, or any combination
%        thereof which is to be the template for a zeroed-out copy to be
%        returned.
%
% OUTPUT
%   B    Data which has the same format as A, but all contents are zeroed
%        (or for strings, empty). The dimensions are maintained, and sparse
%        inputs remain sparse. Unrecognized data types emit a warning and
%        are turned into the empty vector.
%
% EXAMPLE
%
%   acz = zerodup(ac(1));
%   acp = repmat(acz, [1 length(ind.a)]);
%

  B = [];

  % Descend into structs and work on the fields
  if isstruct(A)
    fnames = fieldnames(A(1));
    for ii=1:length(fnames)
      % Do all the work of handling fields in a struct by just recursively
      % calling ourself.
      B.(fnames{ii}) = zerodup(A(1).(fnames{ii}));
    end

  % Cell arrays need to be handled element-by-element since types internally
  % can differ.
  elseif iscell(A)
    for ii=1:length(A)
      % Again, just call ourself recursively and let the top-level cases
      % handle the contents
      B{ii} = zerodup(A{ii});
    end

  % Sparse arrays are more special than a dense array, so capture that case
  % first.
  elseif issparse(A)
    ss = size(A);
    B = sparse(ss(1), ss(2));

  % Zeros don't make sense for character strings, so handle it specially
  elseif ischar(A)
    B = '';

  % Final handled case is just a regular, numeric array.
  elseif isnumeric(A)
    B = zeros(size(A), class(A));

  % Other things should emit a warning about unrecognized data types. Default
  % to just giving back an empty array in this case.
  else
    B = [];
    warning('zerodup:DataTypeError', sprintf(...
      'Encountered unhandled data type `%s`; setting to empty list.', ...
      class(A)))
  end

end
