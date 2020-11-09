function ffbm_makehandcuts(compopt)
% function ffbm_makehandcuts(compopt)
%
% Function to plot up windowed maps det-by-det and visually inspect the
% beam, its fit, and residuals.  Can then choose to keep or discard.
% Requires a cut struct to have been made, and modifies the 'hand' field 
%
% As always, relies on bmrunlist_201?.csv
%
% INPUTS (should all be sent in with compopt)
%        Works best if you use the compopt saved with the cut struct
%
%   expt:            'b2','keck','b3'
%   year:            b2: 2012 only
%                    keck: 2014, 2015
%                    b3:   2015
% OPTIONAL INPUTS
%   
%   mapcoord:        'azel_ap' (default), 'raw' (used to get map filename)
%   demodtype:       'square' (default), 'choplet'
%   size:            map size in deg (2 default), used for map filename
%   suffix:          optional string to add to map name, '' default
%   bmdir:           directory from which to load beammaps
%                    (default beammaps/maps_windowed)
%   mapstoload:      vector of beam map run numbers to look at, i.e. 1:10
%                    (default is everything in bmrunlist.csv)
%   cutdir:          directory from which to load the cut structure 
%                    (default beammaps/cuts)
%   cutfile:         name of cut file to load
%                    (default ffbm_cuts_year)
%   cutlist:         criteria on which to cut
%
% This will load up the cut structure and evaluate the cuts before showing
% you per-detector maps which pass all cuts.  You can then assign each map a
% number corresponding to what you think it is:
%   0 = good beam and good fit
%   1 = good beam and bad fit - note that this means we might still use
%       it in a composite if the fit just chose a weird value - so be
%       conservative here!  If in doubt give it a 2
%   2 = bad beam, who cares about the fit (this is the case if there is a
%       missing stripe of data)
%   3 = no beam
% 
% Then the hand cut field is saved to the file you loaded.  If the field
% already exists for a particular map, it is overwritten.  Note that each
% map will have you looking at 2-3 rxs worth (for Keck) of detectors, so you
% may not be able to do all maps in one sitting!  Control this with the
% mapstoload input.
% 
% This should be fun...like some sort of super nerdy video game
% To make hand cuts for beam map 1 of 2014:
% >> load('beammaps/cuts/ffbm_cuts_2014.mat');
% >> compopt = cuts.compopt;
% >> compopt.mapstoload = 1;
% >> ffbm_makehandcuts(compopt);

% Parse compopt
expt = compopt.expt;
if strcmp('expt','b2') 
  year = 2012;
else
  year = compopt.year;
end
if ~isfield(compopt,'mapcoord');
  mapcoord = 'azel_ap';
else
  mapcoord = compopt.mapcoord;
end
if ~isfield(compopt,'demodtype');
  demodtype = 'square';
else
  demodtype = compopt.demodtype;
end
if ~isfield(compopt,'size');
  sizestr = '2deg';
else
  sizestr = [strrep(num2str(compopt.size),'.','p') 'deg'];
end
if ~isfield(compopt,'suffix')
  suffix = '';
else
  suffix = compopt.suffix;
end
if ~isfield(compopt,'bmdir')
  bmdir = 'beammaps/maps_windowed/';
else
  bmdir = compopt.bmdir;
end
% Which maps do we want?
% Get beam map run info
bm = get_bm_info(year);
if ~isfield(compopt,'mapstoload')
  compopt.mapstoload = bm.number;
else
  bm = structcut(bm,compopt.mapstoload);
end
if ~isfield(compopt,'cutdir')
  compopt.cutdir = 'beammaps/cuts/';
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

[p ind] = get_array_info([num2str(year) '0201']);

% Load up cut struct
cutname = [cutdir '/' cutfile];
load(cutname);
cut = cuts.cuts;

% No need for the cut struct and mapstoload to have the same length

% Load up each beam map requested
for ii = 1:length(bm.number)

  % Annoying indexing since cuts struct may not have all maps
  mapnum = bm.number(ii);
  cutind = find(cuts.bm.number == mapnum);
  % Need to use cuts.bm here, since bm has already been structcut
  filename = cuts.bm.filename{cutind};
  disp(['Loading map ' num2str(ii) ': ' filename]);
  wmapname = [bmdir,'/mapwin_',filename,'_',demodtype,...
	'_',mapcoord,'_',sizestr,suffix];
  w = load(wmapname);
  w = w.bm;

  % Add up all the cuts
  % If we think a lot of (probably good) beams are unnecessarily 
  % penalized by bad fits, in compopt probably want to send in 
  % something like {'mirror','notlight','nanfit'} instead of everything
  totcut = zeros(length(p.gcp),1);
  for jj = 1:length(cutlist)
    % If hand cuts already exist, ignore them
    if ~strcmp(cutlist{jj},'hand')
      totcut = totcut + getfield(cut(bm.number(ii)),cutlist{jj});
    end
  end
  goodind = find(totcut == 0);
  
  % Accept only pairs in which both A and B pass
  good_a = find(ismember(ind.la,goodind));
  good_b = find(ismember(ind.lb,goodind));
  good_pair = intersect(good_a,good_b);
  
  % Plot all pairs passing these cuts...mostly copied from bm_plotter
  F = w.A;
  x_win = w.x_bin - mean(w.x_bin);
  y_win = w.y_bin - mean(w.y_bin);
  
  fh = figure(1); clf;
  set(fh,'Position',[10 10 1400 750]);
  
  for jj = 1:length(good_pair)

    det_a = ind.la(good_pair(jj));
    det_b = ind.lb(good_pair(jj));
    
    disp(['Beam map ' int2str(bm.number(ii)) ': Pair ' int2str(jj) ...
	  ' / ' int2str(length(good_pair))]);
    ro = int2str(p.det_row(det_a));
    co = int2str(p.det_col(det_a));
    ti = int2str(p.tile(det_a));
    rx = int2str(p.rx(det_a));
    gcp_a = int2str(p.gcp(det_a));
    gcp_b = int2str(p.gcp(det_b));
    
    Z_A = gauss2d(F(:,det_a),w.x_bin,w.y_bin,'std');
    Z_B = gauss2d(F(:,det_b),w.x_bin,w.y_bin,'std');
    W_A = w.map(:,:,det_a);
    W_B = w.map(:,:,det_b);
    m_a = w.A(1,det_a);
    m_b = w.A(1,det_b);
    sig_a = sqrt(((w.A(4,det_a))^2 + (w.A(5,det_a))^2)/2);
    sig_b = sqrt(((w.A(4,det_b))^2 + (w.A(5,det_b))^2)/2);
    ell_a = ((w.A(4,det_a))^2 - (w.A(5,det_a))^2) / ...
            ((w.A(4,det_a))^2 + (w.A(5,det_a))^2);
    ell_b = ((w.A(4,det_b))^2 - (w.A(5,det_b))^2) / ...
            ((w.A(4,det_b))^2 + (w.A(5,det_b))^2);
    
    clim_sum  = get_clim_sum(expt,p,det_a,jj,filename);
    
    % A detector

    % Map, dB scale
    subplot(3,5,1)
    imagescnan(x_win,y_win,10*log10(abs(W_A)/abs(m_a)));
    title(['Rx ' rx ' Ti ' ti ' Ro ' ro ' Co ' co ' GCP ' gcp_a ...
	  ' dB'])
    prettify_large([-30 0])
    if m_a < 0
      text(-2,0,'NEGATIVE AMPLITUDE')
    end
    
    % Map, linear scale
    subplot(3,5,2)
    imagescnan(x_win,y_win,W_A);
    title('A Linear, fit peak')
    prettify_large([0 abs(m_a)]);
    
    % Beam fit
    subplot(3,5,3)
    imagescnan(x_win,y_win,Z_A)
    prettify_large([0 abs(m_a)])
    title('A Fit, fit peak')
    
    % Residual
    subplot(3,5,4)
    imagescnan(x_win,y_win,(W_A-Z_A)./m_a)
    prettify_large([-0.05 0.05])
    title('A Residual')
    
    % Map, linear scale
    subplot(3,5,5)
    imagescnan(x_win,y_win,W_A);
    title('A Linear, abs scale')
    prettify_large(clim_sum)

    % B detector
    
    % Map, log scale
    subplot(3,5,6)
    imagescnan(x_win,y_win,10*log10(abs(W_B)/abs(m_b)));
    title(['Rx ' rx ' Ti ' ti ' Ro ' ro ' Co ' co ' GCP ' gcp_b ...
	  ' dB'])
    prettify_large([-30 0])
    if m_b < 0
      text(-2,0,'NEGATIVE AMPLITUDE')
    end
    
    % Map, linear scale
    subplot(3,5,7)
    imagescnan(x_win,y_win,W_B);
    title('B Linear, fit peak')
    prettify_large([0 abs(m_b)])
    
    % Beam fit
    subplot(3,5,8)
    imagescnan(x_win,y_win,Z_B)
    prettify_large([0 abs(m_b)])
    title('B Fit, fit peak')
    
    % Residual
    subplot(3,5,9)
    imagescnan(x_win,y_win,(W_B-Z_B)./m_b)
    prettify_large([-0.05 0.05])
    title('B Residual')
    
    % Map, linear scale
    subplot(3,5,10)
    imagescnan(x_win,y_win,W_B);
    title('B Linear, abs scale')
    prettify_large(clim_sum)
    
    % Diff
    
    subplot(3,5,12)
    imagescnan(x_win,y_win,W_A./m_a - W_B./m_b)
    prettify_large([-0.02 0.02])
    title('Real Difference')
    
    subplot(3,5,13)
    imagescnan(x_win,y_win,Z_A./m_a - Z_B./m_b)
    prettify_large([-0.02 0.02])
    title('Fit Difference')
    
    % Statistics
    
    subplot(3,5,11)
    axis([-2 2 -2 2])
    text(-1.8,1.3,['A amp: ' sprintf('%3.0f',m_a)]);
    text(0.2,1.3,['B amp: ' sprintf('%3.0f',m_b)]);
    text(-1.8,0,['A sig: ' sprintf('%0.3f',sig_a)]);
    text(0.2,0,['B sig: ' sprintf('%0.3f',sig_b)]);
    text(-1.8,-1.3,['A ell: ' sprintf('%0.2f',ell_a)]);
    text(0.2,-1.3,['B ell: ' sprintf('%0.2f',ell_b)]);
    axis off
    axis tight

    subplot(3,5,14)
    hist(log10(abs(W_A(:))),20);
    title('log10(abs(A)) histogram')
    axis square
    axis tight
    
    subplot(3,5,15)
    hist(log10(abs(W_B(:))),20);
    title('log10(abs(B)) histogram')
    axis square
    axis tight
    
    x = input(['0 good beam/fit, '...
               '1 good beam/bad fit (can use in composite), '...
               '2 bad beam, 3 no beam: ']);
    if isempty(x)
      x = 0;
    end
    cut(cutind).hand(det_a) = x;
    cut(cutind).hand(det_b) = x;
    
    clf;
  end
  
end

% Prep output struct
cuts.cuts = cut;

save(cutname,'cuts')

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function clim_sum = get_clim_sum(expt,p,det,ii,filename)
% Color axes are unfortunately experiment and frequency-specific, so choose
% them here

switch expt
  case 'b2'
    clim_sum = [0 300];
  case 'keck'
    switch filename(1:4)
      case '2012'
        clim_sum = [0 150];
      case '2013'
        clim_sum = [0 300];
      case {'2014','2015'}
        switch p.band(det)
          case {100,220}
            clim_sum = [0 200];
          case 150
            clim_sum = [0 300];
        end
      case '2016'
        switch p.band(det)
          case {100,220}
            clim_sum = [0 700];
          case 150
            clim_sum = [0 900];
        end
    end
  case 'b3'
    switch filename(1:4)
      case '2015'
        % Default color limit based on chopper type - before 2015-03-17 we
        % had the 24" and after the 18"
        bmtime = date2mjd(str2num(filename(1:4)),str2num(filename(5:6)),...
                          str2num(filename(7:8)),str2num(filename(10:11)),...
                          str2num(filename(12:13)),str2num(filename(14:15)));
        if bmtime < date2mjd(2015,03,17,00,00,00)
            clim_sum = [0 200];
        else
            clim_sum = [0 100];
        end
      case '2016'
        clim_sum = [0 150];
    end
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function prettify_large(cax)

set(gca,'YDir','normal')
colorbar('location','eastoutside')
caxis(cax)
axis equal
axis tight

return