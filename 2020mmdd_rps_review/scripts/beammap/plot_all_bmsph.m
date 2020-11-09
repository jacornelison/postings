function plot_all_bmsph(mapfile, outdir, rx, clim)
%Plot beam maps in Lambert equal area projection
%
%mapfile - .mat file created by reduc_*_bmsph
%outdir - directory to place map images
%rx - receiver # (0--4)
%clim - optional, color axis limits for plot
  if ~exist('clim', 'var')
    clim = [-1 4];
  end

  mkdir(outdir);
  mkdir([outdir '/full']);
%  [p ind]=get_array_info('20131231'); %Get from file
  m = load(mapfile);
  p = m.p;

  figure('Visible','off')
  close all
  colormap('JET')
  for i=1:528
    ip = i + 528 .* rx;
    clf
    viewer_bmsph(m, i, 'Log', 'Graticule', [10, 90], 'NPix', 101, 'CLim', clim)
    ro=int2str(p.det_row(ip));
    co=int2str(p.det_col(ip));
    ti=int2str(p.tile(ip));
    img_str=strcat(outdir,'/full/map_row_',ro,'_col_',co,'_tile_',ti,'_',p.pol(ip),'.png');
    title_str=strcat('Row ',ro,' Col ',co,' Tile ',ti,' Pol ',p.pol(ip),' Log Scale');
    title(title_str{1},'fontsize',18)
    print('-dpng', '-r100', img_str{:})
  end
end
