function d=strip_tod(d)
% d=strip_tod(d)
% remove unnecessary fields from the TOD data structure

for ff={'array','scancheck','rawpoint',...
         'relgains','relgains_err','scan','tracker'}
    if(isfield(d,ff))
      d=rmfield(d,ff);
    end
end

if(isfield(d.pointing,'cel'))
  d.pointing=rmfield(d.pointing,'cel');
end

%for ff={'az','el'}
for ff={'el'}
  if(isfield(d.pointing.hor,ff))
    d.pointing.hor=rmfield(d.pointing.hor,ff);
  end
end

for ff={'time','tracker'}
  if(isfield(d.antenna0,ff))
    d.antenna0=rmfield(d.antenna0,ff);
  end
end

return