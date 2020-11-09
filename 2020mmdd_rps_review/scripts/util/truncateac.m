function ac=truncateac(m1,m2,ac)
% ac=truncateac(m1,m2,ac)
%
% Undo the action of expandac().
%
% Take ac with the expanded size of defn m2 and truncate back down to the map
% size of defn m1.

if(m1.pixsize~=m2.pixsize)
  error('pixsize must be equal');
end

if isfield(ac, 'wsum')
  sz = size(ac(1).wsum);
else
  sz = size(ac(1).T); % map instead
end
if any(sz ~= [m2.ny,m2.nx])
  error('ac does not match m1');
end

% check m1 lies wholly inside m2
%if(length(intersect(m1.x_tic,m2.x_tic))~=m1.nx|length(intersect(m1.y_tic,m2.y_tic))~=m1.ny)
%  error('m1 not wholly inside m2');
%end

% allow for machine precision differences in the m1.x_tic and m2.x_tic
xpre = find(abs(m1.x_tic(1)-m2.x_tic)<1e-10);
xend = find(abs(m1.x_tic(end)-m2.x_tic)<1e-10);

% allow for machine precision difference in the m1.y_tic and m2.y_tic
ypre = find(abs(m1.y_tic(1)-m2.y_tic)<1e-10);
yend = find(abs(m1.y_tic(end)-m2.y_tic)<1e-10);

% truncate
for i=1:numel(ac)
  for fi=fieldnames(ac)'
    % special treatmeant for deprojection templates
    if iscell(ac(i).(fi{1}))
      ac(i).(fi{1}) = cellfunc(@(acc) acc(ypre:yend,xpre:xend),...
                               ac(i).(fi{1}));
    else
      ac(i).(fi{1}) = ac(i).(fi{1})(ypre:yend,xpre:xend);
    end
  end
end

return
