function out = reformMulti(in, n_samples)

% NAME:
%    reformMult
%
% PURPOSE:
%  function to reform a strcuture of M repeated arrays. (ie
% d.antenna0.bolo.mag0,1,2,3...each double[NX x 1]) into a single array that
% has size double[Nx x M]
%
% CALLING SEQUENCE:
%    out =  reformMulti(in, n_samples)
%
% INPUTS:
%      - in: a structure made of M repeated arrays.
%      - n_samples: number of arrays to make the final structure in ie ( M=n_samples)
%
% DEPENDENCIES:
%    - assumes the M repeated input arrays start with a name with first
%    index 0. (ie: mag0 mag1 mag2)
%
% MODIFICATION HISTORY
%    2007-05-10 Written. (dB)
%    2008-01-25 Bug fix (added possibility of non equal n_rows and made
%    the final array  n_samples x M insead. n_samples is n_elements of utcfast



names= fieldnames(in);
n_elements=size(names, 1);
n_rows=zeros(n_elements,1);
for i=1:n_elements
n_rows(i,1)  = eval(['size(in.' names{i} ',1); ']);
end

out =  zeros(n_samples, n_elements);

for r = 1:n_elements
     eval(['out(1: n_rows(r,1) ,r) = in.' char(names{r}) ';']);
end



