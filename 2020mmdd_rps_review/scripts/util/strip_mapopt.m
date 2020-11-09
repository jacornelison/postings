function mapopt=strip_mapopt(mapopt)
% mapopt=strip_mapopt(mapopt)
%
% Removes diagnostic fields from mapopt as a space-saving technique.
%

  % Cut statistics:        c
  % Halfscan statistics:   devhist, hs
  % Poly filter stats:     filtc
  % Telescope description: p
  % Scan motion info:      traj
  rmflds = intersect(fieldnames(mapopt), ...
    {'c','devhist','filtc','hs','p','traj'});
  if isempty(rmflds)
    return
  end
  mapopt = rmfield(mapopt, rmflds);
end

