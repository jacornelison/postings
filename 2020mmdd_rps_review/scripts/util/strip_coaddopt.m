function coaddopt=strip_coaddopt(coaddopt)
% coaddopt=strip_coaddopt(coaddopt)
%
% Removes diagnostic fields from coaddopt and coaddopt.mapopt as a
% space-saving technique.
%

  % Deprojection:        b, bi, bw
  % Cut statistics:      c
  % Halfscan statistics: devhist, hsmax, whist
  % Poly filter stats:   filtc
  % Scan motion info:    traj
  % Jackknife mask:      jackmask
  rmflds = intersect(fieldnames(coaddopt), ...
    {'b','bi','bw','c','devhist','filtc','hsmax','jackmask','traj','whist'});
  if isempty(rmflds)
    return
  end

  coaddopt = rmfield(coaddopt, rmflds);

  % Also strip out any diagnostic information from mapopt.
  if isfield(coaddopt, 'mapopt')
    if ~iscell(coaddopt.mapopt)
      coaddopt.mapopt{1} = coaddopt.mapopt;
    end
    coaddopt.mapopt = cellfunc(@strip_mapopt, coaddopt.mapopt);
  end

end

