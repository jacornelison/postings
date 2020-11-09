function plot_all_small_lame(mapfile, outdir, rx, clim)
  if ~exist('clim', 'var')
    clim = [-1 4];
  end

  mkdir(outdir);
  mkdir([outdir '/small']);
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
    lame_viewer2(holdmap, i, 'Log', 'NPix', 101, 'CLim', clim, 'ThetaMax', 30)
    colorbar('off');
    ro=int2str(p.det_row(ip));
    co=int2str(p.det_col(ip));
    ti=int2str(p.tile(ip));
    img_str=strcat(outdir,'/small/map_row_',ro,'_col_',co,'_tile_',ti,'_',p.pol(ip),'.png');
    print('-dpng',img_str{:})
  end
  system(['mogrify -resize ''20%x20%'' ' outdir '/small/*'])

end
