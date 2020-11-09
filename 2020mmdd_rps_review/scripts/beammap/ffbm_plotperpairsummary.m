function ffbm_plotperpairsummary(compopt)
% ffbm_plotperpairsummary(compopt)
%
% Function to plot components and fits all in one
% Loads in component maps and fits
% Must use cuts, only uses common AB component maps
% Can also plot a thumbnail with a differential parameter
% 
% INPUTS (should all be sent in with compopt)
%
%   expt:              'keck','b3'
%   year:              keck: 2014, 2015
%
% OPTIONAL INPUTS
%
%   plotdir:         directory in which to save plots
%   componentdir:    directory from which to load component struct
%                    (default beammaps/maps_component)
%   componentfile:   name of file in componentdir to load
%                    (default ffbm_year_allcomp)
%   cutdir:          directory from which to load the cut structure
%                    (default beammaps/cuts)
%   cutfile:         name of cut file to load
%                    (default ffbm_cuts_year)
%   cutlist:         criteria on which to cut 
%                    (default all in cutfile)
%   beamparamdir:    directory from which to load beam params (for thumbnail)
%                    (default beammaps/beamparams)
%   beamparamfile:   name of file in beamparamdir to load
%                    (default beamparams_year)
%   thumbnail:       'dsigma','dellip','dpointing'

% Parse compopt
expt = compopt.expt;
year = compopt.year;
if ~isfield(compopt,'componentdir')
  compopt.componentdir = 'beammaps/maps_component';
  componentdir = compopt.componentdir;
else
  componentdir = compopt.componentdir;
end
if ~isfield(compopt,'componentfile')
  compopt.componentfile = ['ffbm_' num2str(year) '_allcomp'];
  componentfile = compopt.componentfile;
else
  componentfile = compopt.componentfile;
end
if ~isfield(compopt,'cutdir')
  compopt.cutdir = 'beammaps/cuts';
  cutdir = compopt.cutdir;
else
  cutdir = compopt.cutdir;
end
if ~isfield(compopt,'cutfile')
  compopt.cutfile = ['ffbm_cuts_' num2str(year)];
  cutfile = compopt.cutfile;
else
  cutfile = compopt.cutfile;
end
if ~isfield(compopt,'cutlist')
  compopt.cutlist = {'mirror','notlight','nanfit','peak',...
      'sigma','ellip','hand'};
  cutlist = compopt.cutlist;
else
  cutlist = compopt.cutlist;
end
if ~isfield(compopt,'beamparamdir')
  compopt.beamparamdir = 'beammaps/beamparams';
  beamparamdir = compopt.beamparamdir;
else
  beamparamdir = compopt.beamparamdir;
end
if ~isfield(compopt,'beamparamfile')
  compopt.beamparamfile = ['ffbm_' num2str(year)];
  beamparamfile = compopt.beamparamfile;
else
  beamparamfile = compopt.beamparamfile;
end
if ~isfield(compopt,'plotdir')
  compopt.plotdir = 'plots/';
  plotdir = compopt.plotdir;
else
  plotdir = compopt.plotdir;
end
if ~isfield(compopt,'thumbnail')
  compopt.thumbnail = [];
  thumbnail = compopt.thumbnail;
else
  thumbnail = compopt.thumbnail;
end

switch year
  case 2014
    [p ind] = get_array_info(20140201);
  case 2015
    [p ind] = get_array_info(20150201);
end

% Load up cuts
load([cutdir '/' cutfile]);
cut = cuts.cuts;

% Load up component map file
load([componentdir '/' componentfile]);

% Load up beam param file
load([beamparamdir '/' beamparamfile]);
params = ffbm_evalbeamparam(compopt);

for ii = 1:length(ind.la)
  
  if mod(ii,10) == 0
    disp(['Plotting pair ' int2str(ii) ' / ' int2str(length(ind.la))]);
  end

  % Make a total cut mask
  totmask_a = ones(length(comp.bm.number),1);
  totmask_b = ones(length(comp.bm.number),1);
  
  for jj = 1:length(comp.bm.number);
    a = structcut(cut(jj),ind.la(ii)); % all cuts for this detector
    b = structcut(cut(jj),ind.lb(ii)); 
    totcut_a = 0;
    totcut_b = 0;
    for kk = 1:length(cutlist)
      totcut_a = totcut_a + getfield(a,cutlist{kk});
      totcut_b = totcut_b + getfield(b,cutlist{kk});
    end
    if totcut_a == 0
      totmask_a(jj) = 0;
    end
    if totcut_b == 0
      totmask_b(jj) = 0;
    end
  end
    
  % Indices of good maps
  plotinds_a = find(totmask_a == 0);
  plotinds_b = find(totmask_b == 0);
  
  % Must use common AB
  plotinds_a = intersect(plotinds_a,plotinds_b);
  plotinds_b = plotinds_a;

  ro = p.det_row(ind.la(ii));
  co = p.det_col(ind.la(ii));
  ti = p.tile(ind.la(ii));
  filebase = [plotdir '/full/det_row_' num2str(ro) '_col_' num2str(co) ...
	'_tile_' num2str(ti) '_dif'];
  
  % Plot 6-panel for components
  for jj = 1:length(plotinds_a)
    p0 = structcut(p,[ind.la(ii) ind.lb(ii)]);
    plot_compsub(comp.map{ind.la(ii)}.component{plotinds_a(jj)},...
	comp.map{ind.lb(ii)}.component{plotinds_b(jj)},...
	comp.map{ind.la(ii)}.fit{plotinds_a(jj)},...
	comp.map{ind.lb(ii)}.fit{plotinds_b(jj)},...
	p0,-comp.bm.dk(plotinds_a(jj))-5.25,...
	comp.x_bin{plotinds_a(jj)},comp.y_bin{plotinds_a(jj)});
    %disp([int2str(ii) ' ' int2str(jj)])

    savename = [filebase '_bm' num2str(plotinds_a(jj)) ];
    print('-dpng',savename);
  end

  filestosmoosh =  [filebase '_bm*'];
  % Montage all panels together
  system(['montage ' filestosmoosh ' -geometry 600x500+10+10 ' filebase '.png']);

  if ~isempty(thumbnail)
    savename = [plotdir '/small/det_row_' num2str(ro) '_col_' num2str(co) ...
	  '_tile_' num2str(ti) '_dif_sm.png'];
    param_onepair = structcut(params,ind.la(ii));
    plot_thumb(thumbnail,param_onepair);
    print('-dpng',savename);
    system(['mogrify -crop ''500x500+95+70'' ' savename]);
    system(['mogrify -resize ''10%x10%'' ' savename]);
  end
  

end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plot_compsub(mapA,mapB,fitA,fitB,p,dk,x_bin,y_bin)

% Position vector is [left,bottom,width,height]
% in units of the respective dimension
% Following is appropriate for a 2rowX3col panel of square
width = 0.333;
height = 0.5;

% Beam params
% Fits are evaluated for dk = 0
% Rotate to x'/y' (rotangle = p.drumangle + dk) 
[A.sig A.e_p A.e_c] = ...
    ffbm_findingparam(fitA,p,dk,1,'b3');
[B.sig B.e_p B.e_c] = ...
    ffbm_findingparam(fitB,p,dk,1,'b3');

% Prep for fits
ZA = gauss2d(fitA',x_bin,y_bin,'std');
ZB = gauss2d(fitB',x_bin,y_bin,'std');
% Map amplitudes (a la bm_plotter);
mA = fitA(1);
mB = fitB(1);

% Rotate maps to x'/y' for plotting
% Raw maps are at dk=0 so add back in the drumangle
%mapA = imrotate(mapA,p.drumangle(1)+180,'bilinear','crop');
%mapB = imrotate(mapB,p.drumangle(1)+180,'bilinear','crop');
%ZA = imrotate(ZA,p.drumangle(1)+180,'bilinear','crop');
%ZB = imrotate(ZB,p.drumangle(1)+180,'bilinear','crop');

fh = figure(1); clf;
setwinsize(fh,500,333)
set(gcf,'Visible','off')

% Top left: A map
axes('Position',[0,height,width,height])
imagesc(x_bin,y_bin,log10(abs(mapA)));
set(gca,'YDir','normal')
caxis([log10(mA)-3 log10(mA)])
axis off

% Bottom left: B map
axes('Position',[0,0,width,height])
imagesc(x_bin,y_bin,log10(abs(mapB)));
set(gca,'YDir','normal')
caxis([log10(mB)-3 log10(mB)])
axis off

% Top middle: A residual
axes('Position',[width,height,width,height])
imagesc(x_bin,y_bin,(mapA - ZA)./mA);
set(gca,'YDir','normal')
caxis([-0.1 0.1])
axis off

% Bottom middle: B residual
axes('Position',[width,0,width,height])
imagesc(x_bin,y_bin,(mapB - ZB)./mB);
set(gca,'YDir','normal')
caxis([-0.1 0.1])
axis off

% Top right: Fit contours and values
axes('Position',[width*2,height,width,height])
[XX YY] = meshgrid(x_bin,y_bin);
contour(XX,YY,ZA./mA,[0.68 0.34 0.05],'LineColor','r'); hold on;
contour(XX,YY,ZB./mB,[0.68 0.34 0.05],'LineColor','b');
plot(fitA(2),fitA(3),'r+');
plot(fitB(2),fitB(3),'b+');
text(-1.9+mean(x_bin),1.7+mean(y_bin),sprintf('A Sig: %.3f',A.sig));
text(-1.9+mean(x_bin),1.4+mean(y_bin),sprintf('B Sig: %.3f',B.sig));
text(-1.9+mean(x_bin),1.1+mean(y_bin),sprintf('d Sig: %.3f',A.sig - B.sig));
%text(0.3,1.7,sprintf('dX: %.3f',A.x - B.x));
%text(0.3,1.4,sprintf('dY: %.3f',A.y - B.y));
text(-1.9+mean(x_bin),-1.1+mean(y_bin),sprintf('A E_p: %.3f',A.e_p));
text(-1.9+mean(x_bin),-1.4+mean(y_bin),sprintf('B E_p: %.3f',B.e_p));
text(-1.9+mean(x_bin),-1.7+mean(y_bin),sprintf('d E_p: %.3f',A.e_p - B.e_p));
text(0.3+mean(x_bin),-1.1+mean(y_bin),sprintf('A E_c: %.3f',A.e_c));
text(0.3+mean(x_bin),-1.4+mean(y_bin),sprintf('B E_c: %.3f',B.e_c));
text(0.3+mean(x_bin),-1.7+mean(y_bin),sprintf('d E_c: %.3f',A.e_c - B.e_c));
set(gca,'XTick',[])
set(gca,'YTick',[])

% Bottom right: Difference residual
axes('Position',[width*2,0,width,height])
imagesc(x_bin,y_bin,mapA./mA - mapB./mB);
set(gca,'YDir','normal')
caxis([-0.1 0.1])
colorbar('Position',[width*2 0.05 0.02 0.4]);
text(0.3+mean(x_bin),-1.7+mean(y_bin),['dk: ' num2str(dk)]);
axis off

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plot_thumb(thumbnail,param_onepair)

fh = figure(2); clf;
setwinsize(fh,400,400)
set(fh,'Visible','off')

switch thumbnail
  case 'dsigma'
    scale = 15000;
    if isnan(param_onepair.dsig_med)
      param_onepair.dsig_med = 1e-9;
    end
    q = plot(0,0,'o','MarkerSize',abs(param_onepair.dsig_med)*scale);
    if param_onepair.dsig_med > 0
      set(q,'MarkerFaceColor','r')
    else
      set(q,'MarkerFaceColor','b')
    end
    xlim([-1 1])
    ylim([-1 1])
    axis equal
  case 'dpointing'
  case 'dellip'
  case 'dsigma_dpointing'
    scale = 15000;
    if isnan(param_onepair.dsig_med)
      param_onepair.dsig_med = 1e-9;
    end
    q = plot(0,0,'o','MarkerSize',abs(param_onepair.dsig_med)*scale);
    hold on;
    if param_onepair.dsig_med > 0
      set(q,'MarkerFaceColor','r')
    else
      set(q,'MarkerFaceColor','b')
    end
    scale = 60;
    q = quiver(0,0,param_onepair.aboffset_x*scale,...
	param_onepair.aboffset_y*scale,0,'filled','AutoScale','off');
    set(q,'MarkerSize',4);
    set(q,'LineWidth',6);
    set(q,'ShowArrowHead','on');
    set(q,'MaxHeadSize',4);
    set(q,'Color','black');
    set(gca,'YDir','normal');
    xlim([-1 1])
    ylim([-1 1])
    axis equal
    ylim([-1 1])
end

return