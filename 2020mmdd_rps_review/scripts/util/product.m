function combos=product(varargin)
% combos=product(varargin)
%
% Returns a cell-of-cells vector where each element cell array is a unique
% combination of one element from each of the input arguments.
%
% This is similar to the Cartesian product of iterators in languages which
% have iterators.
%
% EXAMPLE
%
%   colors = {'black','red','green','blue'};
%   tf = [true, false];
%   optcombos = combinations(colors, tf, tf, tf);
%   for ii=1:length(optcombos)
%     [color,bold,underline,italic] = deal(optcombos{ii}{:});
%     fancyhtml(color, bold, underline, italic);
%   end
%

  % The number of input arguments determines the number of dimensions in the
  % output, and the size of each dimension corresponds to the number of
  % values given in each argument.
  ndims = numel(varargin);

  % Encapsulate each element in an array in it's own cell. This is required
  % for the loop later, and it puts cell arrays and numeric arrays on the
  % same footing so that extra branches can be avoided later.
  for ii=1:ndims
    if iscell(varargin{ii})
      varargin{ii} = cellfunc(@(v) {v}, varargin{ii});
    else
      varargin{ii} = arrayfunc(@(v) {v}, varargin{ii});
    end
    % Force everything to be a row vector for consistency.
    varargin{ii} = reshape(varargin{ii}, 1, []);
  end

  % Initialize the combinations output with just the entries of the first
  % argument.
  combos = varargin{1}';

  % Then for all remaining dimensions, the basic scheme is:
  %
  %   1. Replicate the existing combinations from a column vector to a
  %      matrix, where the width matches the number of choices in the
  %      next option.
  %
  %   2. Similarly, make the next option row-vector into a matrix, matching
  %      the number of rows in the existing combinations matrix.
  %
  %   3. Then do a cell-by-cell concatenation of the two cell-matrices to
  %      extend the combinations to include the next option.
  %
  %   4. Reshape the combination matrix back into a vector so that we can
  %      return to step 1 of the loop (or return this to the user).
  for ii=2:ndims
    combos = repmat(combos, 1, numel(varargin{ii}));
    nextc  = repmat(varargin{ii}, size(combos,1), 1);
    combos = cellfunc(@horzcat, combos, nextc);
    combos = reshape(combos, [], 1);
  end
end

