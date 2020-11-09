function mapname_new = expandmap(mapname,m_new_name,dry_run)
% mapname_new=expandmap(mapname,m_new,dry_run)
%
% Pad a map to a new map with mapadefn m_new_name,
% and save it to a new map file.
% 
% This is supposed to be used together with reduc_plotcomap_pager
% e.g. 
% mapname_new = expandmap('maps/1351/real_f_filtp3_weight3_gs_dp1100_jack01.mat','bicepext',0);
% reduc_plotcomap_pager(mapname_new,...);

% New map name
x = load(mapname);
x.coaddopt.daughter = [x.coaddopt.daughter '_' m_new_name];
mapname_new = get_map_filename(x.coaddopt);

if dry_run
  return
end

m_new = get_map_defn(m_new_name);

% Padding process
if iscell(x.ac)
  for ii=1:numel(x.ac)
    x.ac{ii} = expandac(x.m,m_new,x.ac{ii});
  end
else
  x.ac = expandac(x.m,m_new,x.ac);
end

x.m = m_new;

% Save
saveandtest(mapname_new,'-struct','x','-v7.3');

end
