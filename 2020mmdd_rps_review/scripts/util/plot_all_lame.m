function plot_all_lame(mapfile, outdir, rx, clim)
  if ~exist('clim', 'var')
    clim = [-1 4];
  end

  mkdir(outdir);
  mkdir([outdir '/full']);
  [p ind]=get_array_info('20131231');
  load(mapfile);

  holdmap.data = permute(map, [2 1 3]);
  holdmap.theta = y_bin; %theta values of map bins
  holdmap.phi = x_bin;

  figure('Visible','off')
  close all
  colormap('JET')
  for i=1:528
    ip = i + 528 .* rx;
    clf
    lame_viewer2(holdmap, i, 'Log', 'Graticule', [10, 60], 'NPix', 101, 'CLim', clim)
    ro=int2str(p.det_row(ip));
    co=int2str(p.det_col(ip));
    ti=int2str(p.tile(ip));
    img_str=strcat(outdir,'/full/map_row_',ro,'_col_',co,'_tile_',ti,'_',p.pol(ip),'.png');
    title_str=strcat('Row ',ro,' Col ',co,' Tile ',ti,' Pol ',p.pol(ip),' Log Scale');
    title(title_str{1},'fontsize',18)
    print('-dpng', '-r100', img_str{:})
  end
end
