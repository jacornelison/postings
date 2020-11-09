function st=desparsify(st)
% st=desparsify(st)
%
% Converts everything within the given structure to dense arrays, including
% recursivley descending into child structures.
%
% INPUTS
%   st    A generic structure which can contain arrays, cell arrays, or
%         further nested structures.
%
% OUTPUTS
%   st    Equivalent to input structure except elements have been turned into
%         dense arrays, including elements within nested structures.
%
% EXAMPLE
%
%   ac = desparsify(ac);
%

  % Recursively call this function for all elements of a structure array.
  if numel(st) > 1
    for ii=1:numel(st)
      st(ii) = desparsify(st(ii));
    end
    return
  end

  fnames = fieldnames(st);

  % We'll only get here for a "scalar" structure
  for f=1:length(fnames)
    ff=fnames{f};

    if isstruct(st.(ff))
      st.(ff) = desparsify(st.(ff));

    elseif iscell(st.(ff))
      for cc=1:length(st.(ff))
        if issparse(st.(ff){cc})
          st.(ff){cc} = full(st.(ff){cc});
        end
      end

    elseif issparse(st.(ff))
      st.(ff) = full(st.(ff));
    end
  end

end
