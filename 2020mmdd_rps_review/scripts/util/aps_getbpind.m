function bpind = aps_getbpind(field1, field2, iscross)
% bpind = aps_getbpind(field1, field2, iscross)
%
% Get the index of a particular spectrum in pipeline bandpower
% arrays. For a single experiment (like B2xB2), there are six
% unique spectra, ordered TT, TE, EE, BB, TB, EB. For the cross
% between two different experiments, there are nine unique spectra
% (alt-crosses).
%
% [Input]
%   field1   String (or cell array) specifying one or more CMB
%            fields: 'T', 'E', or 'B' 
%   field2   Same as field1 argument, except specifying the second
%            field in the bandpower.
%   iscross  0 (false) for single experiments, 1 (true) for
%            experiment cross-spectra
%
% [Output]
%   bpind  Index of the desired spectrum, in range [1:6] for
%          experiment auto or range [1:9] for experiment cross.
%          Unknown spectra are assigned bpind of 0.

% If field1, field2 inputs are not cell arrays, convert them.
if ~iscell(field1)
  field1 = {field1};
end
if ~iscell(field2)
  field2 = {field2};
end

% If field1, field2 inputs have different lengths, expand the
% shorter one by repeating its last value.
nbp = max(numel(field1), numel(field2));
if numel(field1) < nbp
  field1(end+1:nbp) = field1(end);
end
if numel(field2) < nbp
  field2(end+1:nbp) = field2(end);
end

% Pipeline ordering of spectra for experiment auto and cross.
spectra = {'TT', 'TE', 'EE', 'BB', 'TB', 'EB', 'ET', 'BT', 'BE'};
if iscross
  order = [1:9];
else
  order = [1:6, 2, 5, 6];
end

% Determine index for each (field1, field2) pair.
bpind = zeros([1, nbp]);
for ii=1:nbp
  ind = find(strcmp(spectra, [field1{ii}, field2{ii}]));
  if ~isempty(ind)
    bpind(ii) = order(ind);
  end
end
