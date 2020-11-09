function st=sparsify(st)
% st=sparsify(st)
%
% Converts all matrices within the given structure to sparse matrices,
% including recursivley descending into child structures/cells.
%
% INPUTS
%   st    A generic structure which can contain arrays, cell arrays, or
%         further nested structures.
%
% OUTPUTS
%   st    Equivalent to input structure except elements have been turned into
%         sparse arrays, including elements within nested structures.
%
% EXAMPLE
%
%   ac = sparsify(ac);
%

  % Recursively call this function for all elements of a structure array.
  if numel(st) > 1
    for ii=1:numel(st)
      st(ii) = sparsify(st(ii));
    end
    return
  end

  fnames=fieldnames(st);

  % We'll only get here for a "scalar" structure
  for f=1:length(fnames)
    ff=fnames{f};

    if isstruct(st.(ff))
      st.(ff)=sparsify(st.(ff));

    elseif iscell(st.(ff))
      for cc=1:length(st.(ff))
        if isfloat(st.(ff){cc}) && ndims(st.(ff){cc})==2
          st.(ff){cc}=sparse(st.(ff){cc});
        end
      end

    elseif isfloat(st.(ff)) && ndims(st.(ff))==2
      st.(ff)=sparse(st.(ff));
    end
  end

end
