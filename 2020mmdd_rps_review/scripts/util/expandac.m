function ac=expandac(m1,m2,ac,padval)
% ac=expandac(m1,m2,ac,padval)
%
% Take ac with map defn m1 and expand to ac with map defn m2 by
% padding around the edge. padval defaults to 0.
%
% See also truncateac().

if ~exist('padval','var') || isempty(padval)
  padval = 0;
end

if(m1.pixsize~=m2.pixsize)
  error('pixsize must be equal');
end

if isfield(ac, 'wsum')
  sz = size(ac(1).wsum);
else
  sz = size(ac(1).T); % map instead
end
if any(sz ~= [m1.ny,m1.nx])
  error('ac does not match m1');
end

% check m1 lies wholly inside m2
%if(length(intersect(m1.x_tic,m2.x_tic))~=m1.nx|length(intersect(m1.y_tic,m2.y_tic))~=m1.ny)
%  error('m1 not wholly inside m2');
%end

% allow for machine precision differences in the m1.x_tic and m2.x_tic
xpre = find(abs(m1.x_tic(1)-m2.x_tic)<1e-10)-1;
xpost= m2.nx-find(abs(m1.x_tic(end)-m2.x_tic)<1e-10);

% allow for machine precision difference in the m1.y_tic and m2.y_tic
ypre=find(abs(m1.y_tic(1)-m2.y_tic)<1e-10)-1;
ypost=m2.ny-find(abs(m1.y_tic(end)-m2.y_tic)<1e-10);

% expand
for i=1:numel(ac)
  for fi=fieldnames(ac)'
    % special treatmeant for deprojection templates
    if iscell(ac(i).(fi{1}))
      ac(i).(fi{1}) = cellfunc(@(acc) ...
                               padarray(acc,[ypre,xpre],padval,'pre'),...
                               ac(i).(fi{1}));
      ac(i).(fi{1}) = cellfunc(@(acc) ...
                               padarray(acc,[ypost,xpost],padval,'post'),...
                               ac(i).(fi{1}));
    else
      ac(i).(fi{1})=padarray(ac(i).(fi{1}),[ypre,xpre],padval,'pre');
      ac(i).(fi{1})=padarray(ac(i).(fi{1}),[ypost,xpost],padval,'post');
    end
  end
end

return
