% PLOT_TILES  Plot detector parameters in tile view, as laid
% out on the focal plane, using fp_data file for tile layout information.
%
% PLOT_TILES(DATA,P) plots the values in the vector DATA, using fp_data
% information in P.  P can be a loaded P structure from get_array_info;
% a tag or date to pass to get_array_info; or the name of a CSV file
% to be loaded.  If P is not specified, get_array_info will be called
% with no argument.
%
% Additional arguments are optional and can be specified as parameter-
% value pairs:
%
%   PLOT_TILES(...,'fig', gcf()) -- reuse an existing figure
%   PLOT_TILES(...,'clim',[-1 1]) -- set color scale
%   PLOT_TILES(...,'title','Tile view') -- set figure title
%   PLOT_TILES(...,'clab','uK_{CMB}') -- set colorbar title
%   PLOT_TILES(...,'tilenames',{'Tile 1','Tile 2','Tile 3','Tile 4'})
%   PLOT_TILES(...,'rx',0) -- select receiver(s).  One figure per rx.
%   PLOT_TILES(...,'rot',90) -- overall rotation of focal plane view
%   PLOT_TILES(...,'scale',1) -- tweak size of indiviual tile axes
%   PLOT_TILES(...,'pair','diff') -- take pair differences
%       Can be 'ab' (default), 'sum', 'diff', or 'mean'.
%   PLOT_TILES(...,'tile',1) -- Isolate a single tile. 
%
% Note: sometimes Matlab has trouble rendering colors with this
% function.  Try either a later version of Matlab, or another renderer:
% >> set(0,'DefaultFigureRenderer','zbuffer') % or painters, opengl
function h = plot_tiles (dat, p, varargin)

% Load P structure if needed
if ~exist('p','var') || isempty(p)
  p=get_array_info();
elseif ischar(p)
  if exist(p,'file')
    csvfname=p;
    try
      p=ParameterRead(csvfname);
    catch
      error(['Failed to load file ' csvfname ' as valid fp_data CSV file.']);
    end
  else
    tagdate=p;
    p=get_array_info(tagdate);
  end
elseif ~isstruct(p)
  error(['Bad data type for fp_data input P.']);
end

% Default options
S.clim=[];
S.title='';
S.clab='';
S.tilenames={};
S.rx=[];
S.pair='ab';
S.rot = 90;
S.fig = [];
S.tile = [];

% Default axis size tweaks
switch(p.expt)
  case 'bicep3', S.scale = 1.3;
  otherwise, S.scale = 1.1;
end

% Handle any options specified as parameter-value pairs
if length(varargin)==1 && isstruct(varargin{1})
  S2=varargin{1};
  fnames=fieldnames(S2);
  for i=1:length(fnames)
    if isfield(S,fnames{i})
      S.(fnames{i})=S2.(fnames{i});
    else
      error(['Unrecognized option ' fnames{i} '.']);
    end
  end
else
  for i=1:2:length(varargin)
    pname=lower(varargin{i});
    if i+1>length(varargin)
      error(['No value specified for parameter ' pname '.']);
    end
    pval=varargin{i+1};
    if ~isfield(S,pname)
      error(['Unrecognized option ' pname '.']);
    end
    S.(pname)=pval;
  end
end

% Try just cutting all channels that aren't the tile we want.
if ~isempty(S.tile)
    if length(unique(p.tile)) > 1
        if length(dat) == length(p.gcp)
            dat = dat(p.tile==S.tile);
        end
        
        if isfield(p,'expt')
            p = rmfield(p,'expt');
        end
        p = structcut(p,p.tile==S.tile);
    end
end

% Multi-rx logic (loop if more than one rx requested)
if isempty(S.rx)
  S.rx=unique(p.rx(isfinite(p.rx)));
  if length(S.rx)>1
    S2=S;
    h=[];
    for irx=1:length(S.rx)
      ttl = ['rx' num2str(S.rx(irx))];
      if ~isempty(S.title) 
        ttl = [S.title ' ' ttl];
      end
      S2.title=ttl;
      S2.rx=S.rx(irx);
      S2.fig=[];
      tmph= plot_tiles(dat,p,S2);
      h=[h,tmph];
    end
    return
  end
end
% Need to pad data vector?
if length(dat)<length(p.rx)
  if length(dat)==sum(p.rx==S.rx)
    tmpdat=zeros(size(p.rx));
    tmpdat(p.rx==S.rx)=dat;
    dat=tmpdat;
  else
    error(['Data structure has wrong length.']);
  end
end

% Color stretch if not given
if isempty(S.clim)
  S.clim=[min(dat(p.rx==S.rx)),max(dat(p.rx==S.rx))];
  if diff(S.clim)==0
    S.clim=S.clim(1)+[-1 1];
  end
  if any(~isfinite(S.clim))
    S.clim=[0 1];
  end
end

% Font size heuristics
FS1 = get (0, 'DefaultAxesFontSize');
FS2 = floor (FS1 / 2);
if (FS2 < 10)
  if (FS1 >= 10)
    FS2 = 10;
  else
    FS2 = FS1;
  end;
end;
FS3 = floor (FS1 * 2/3);
if (FS3 < 12)
  if (FS1 >= 12)
    FS3 = 12;
  else
    FS3 = FS1;
  end;
end

% Draw title and colorbar, and calculate
% how much space we have left to work with.
if isempty(S.fig)
  S.fig = figure();
  % Make reasonably big by default
  setwinsize(gcf(), 900, 800);
end
h = S.fig;
figure(h);
clf();
% Work around a bug in *at least* 2009a where after plotting 4 tiles, the
% figure suddenly closes, reopens, redraws empty tiles, and then fails to
% correctly draw all later tiles. The default painter is 'opengl', and
% apparently switching to 'painters' fixes the issue.
if ~isempty(strfind(version(), 'R2009a'))
  set(h, 'Renderer', 'painters')
end

if ~isempty(S.title)
  title(S.title);
end
set (gca, 'XTick', []);
set (gca, 'YTick', []);
hcol = colorbar;
if ~isempty(S.clim)
  caxis(S.clim);
end
if ~isempty(S.clab)
  % set (hcol, 'XAxisLocation', 'top');
  title (hcol, S.clab);
  % set (get (hcol, 'XLabel'), 'String', S.clab);
end;
pos = get (gca, 'Position');
max_h = 1 - (1 - (pos(2) + pos(4))) * 0.75;
% Use OuterPosition for < R2014b. Newer only has Position property.
if isprop(hcol, 'OuterPosition')
  pos = get(hcol, 'OuterPosition');
else
  pos = get(hcol, 'Position');
end
% max_w = 1 - pos(3);
max_w = pos(1)*0.95;
axis off;
% tmp_px = p.pix_phys_x;
% tmp_py = p.pix_phys_y;
tmp_px = p.r .* cosd(p.theta + p.drumangle);
tmp_py = p.r .* sind(p.theta + p.drumangle);
px = tmp_px * cosd(S.rot) - tmp_py * sind(S.rot);
py = tmp_py * cosd(S.rot) + tmp_px * sind(S.rot);
pix_phys_box = [...
  min(px(p.rx==S.rx)), ...
  max(px(p.rx==S.rx)), ...
  min(py(p.rx==S.rx)), ...
  max(py(p.rx==S.rx)), ...
];
px = (px - pix_phys_box(1))/diff(pix_phys_box(1:2));
py = (py - pix_phys_box(3))/diff(pix_phys_box(3:4));

% Loop over tiles

tiles=unique(p.tile(p.rx==S.rx));
tiles=tiles(isfinite(tiles));
for i=1:length(tiles)
  ct = (p.rx==S.rx) & (p.tile==tiles(i));
  nrow = max(p.det_row(ct))-min(p.det_row(ct))+1;
  ncol = max(p.det_col(ct))-min(p.det_col(ct))+1;
  cdet = isfinite(px) & isfinite(py) & isfinite(p.det_row) & isfinite(p.det_col);

  % Determine tile orientation
  A = [p.det_col(ct&cdet),p.det_row(ct&cdet)];
  A(:,1)=A(:,1)-nanmean(A(:,1));
  A(:,2)=A(:,2)-nanmean(A(:,2));
  B = [px(ct&cdet),py(ct&cdet)];
  B(:,1)=B(:,1)-nanmean(B(:,1));
  B(:,2)=B(:,2)-nanmean(B(:,2));
  k = A\B;
  k = sign(round(k/max(abs(k(:)))));
  if all(k==[1,0;0,-1])
    rot = 0;
  elseif all(k==[-1,0;0,1])
    rot = 180;
  elseif all(k==[0,-1;-1,0])
    rot = +90;
  elseif all(k==[0,1;1,0])
    rot = -90;
  else
    error(['Unable to determine orientation of tile ' num2str(tiles(i)) '.']);
  end
  
  % Axis outer position for this tile
  lower_margin=max_h/20;
  left_margin=max_w/10;
  max_use_h = max_h - lower_margin;
  max_use_w = max_w - left_margin;
  pos(1)=max_use_w*min(px(ct));
  pos(3)=max_use_w*max(px(ct)) - pos(1);
  pos(2)=max_use_h*min(py(ct));
  pos(4)=max_use_h*max(py(ct)) - pos(2);
  pos = pos + [left_margin, lower_margin, 0, 0];
  
  % disp(['Tile ' num2str(tiles(i)) ' : ax = ' num2str(pos) ', rot = ' num2str(rot)]);
  istop = abs((pos(2)+pos(4))-max_h)<0.1;
  isbot = abs(pos(2)-lower_margin)<0.1;
  isright =  abs((pos(1)+pos(3))-max_w)<0.1;
  isleft = abs(pos(1)-left_margin)<0.1;
  
    
  % Apply scaling tweak
  uu = [max_use_w, max_use_h, max_use_w, max_use_h];
  pos = (pos - [left_margin, lower_margin, 0, 0]) ./ uu;
  pos([1,3])=pos([1,3]) + pos(3)*[-0.5,1] * (S.scale-1);
  pos([2,4])=pos([2,4]) + pos(4)*[-0.5,1] * (S.scale-1);
  uu = 1 - (1 - uu) * S.scale;
  pos = pos .* uu + [left_margin,lower_margin,0,0];

  ax = axes('OuterPosition', pos, 'FontSize', FS2);

  
  % Set up this axis
  switch(rot)
    case   0, nx=ncol; xx=1:ncol;    ny=nrow; yy=nrow:-1:1;
    case  90, nx=nrow; xx=nrow:-1:1; ny=ncol; yy=ncol:-1:1;
    case -90, nx=nrow; xx=1:nrow;    ny=ncol; yy=1:ncol;
    case 180, nx=ncol; xx=ncol:-1:1; ny=nrow; yy=1:nrow;
  end

  % Set up this axis
  axis([0 nx 0 ny]);
  axis square;
  set(ax, 'XTick', 0.5:1:(nx-0.5));
  set(ax, 'YTick', 0.5:1:(ny-0.5));
  if isbot
    set(ax,'XAxisLocation','top');
  end
  if isleft
    set(ax,'YAxisLocation','right');
  end
  set(ax,'XTickLabel', num2str(xx'));
  set(ax,'YTickLabel', num2str(yy'));

  % Draw the title
  use_title = '';
  if length(S.tilenames)>length(tiles)
    use_title = S.tilenames(tiles(i));
  elseif ~isempty(S.tilenames)
    use_title = S.tilenames(i);
  end
  if isempty(use_title)
    use_title = ['Tile ' num2str(tiles(i))];
  end
  hut = text (nx/2, ny+0.5 - (ny+1.5)*(isbot), use_title, ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'FontSize', FS3);
  
  % Define patch polygons
  % LL + UR triangles
  if strcmpi(S.pair,'ab')
    xxll = [[1; 0; 0], [1; 1; 0]];
    yyll = [[0; 0; 1], [0; 1; 1]];
  else
    xxll = [[1; 0; 0; 1]];
    yyll = [[0; 0; 1; 1]];
  end
  % first column
  xx1 = [];
  yy1 = [];
  for (ir = (1:ny)-1)
    xx1 = [xx1, xxll];
    yy1 = [yy1, yyll+(ny-1-ir)];
  end;
  % all columns
  xx = [];
  yy = [];
  for (ic = (1:nx)-1)
    xx = [xx, xx1+ic];
    yy = [yy, yy1];
  end;
  dr = repmat (1:ny, 1, nx);
  dc = repmat (1:nx, ny, 1);
  switch(rot)
    case 90, [dr,dc]=deal(nrow-dc+1,dr);
    case -90, [dr,dc]=deal(dc,ncol-dr+1);
    case 180, [dr,dc]=deal(nrow-dr+1,ncol-dc+1);
  end
  dr = dr(:);
  dc = dc(:);
  
  
  % Indexing
  adat = NaN * zeros(size(dr));
  bdat = NaN * zeros(size(dr));
  for j=1:length(dr)
    ka = find(ct & cdet & p.det_row==dr(j) & p.det_col==dc(j) & strcmp(p.pol,'A'));
    if ~isempty(ka)
      adat(j) = dat(ka(1));
    end
    kb = find(ct & cdet & p.det_row==dr(j) & p.det_col==dc(j) & strcmp(p.pol,'B'));
    if ~isempty(kb)
      bdat(j) = dat(kb(1));
    end
  end
  switch(lower(S.pair))
    case 'diff', abdat=adat-bdat;
    case 'sum', abdat=adat+bdat;
    case 'mean', abdat=(adat+bdat)/2;
    case 'ab', abdat=[adat,bdat]';
  end
  
  patch(xx, yy, abdat(:)', 'EdgeColor', 'black');

  % Color stretch
  if ~isempty(S.clim)
    caxis(S.clim);
  end
end

return
