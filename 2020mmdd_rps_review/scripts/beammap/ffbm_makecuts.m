function ffbm_makecuts(compopt)
% function ffbm_makecuts(compopt)
%
% Function to generate a binary cut structure for a particular year of
% beam map runs.  Relies on bmrunlist.csv for info.  Based on
% ffbm_cutting_pixels.  This loads up windowed and fitted maps one by one
% and cuts based on the fits. 
%
% Current possible cuts:
%
%   mirror:          cuts if main beam falls of mirror 
%                    (Keck only, does nothing for BICEPs)
%   notlight:        cuts if not in ind.l
%   nanfit:          cuts if the fitter failed
%   peak:            cuts fits with extreme peak amplitudes
%   sigma:           cuts fits with extreme beamwidths
%   ellip:           cuts fits with extreme ellipticities
%   abspoint:        cuts fits with beam center too far from origin
%   hand:            cuts based on visual inspection.  This function does
%                    not generate these, but the field is created anyway 
%                    (use ffbm_makehandcuts after running this)
% 
% 20180717 TSG: added abs pointing cut, changed mirror cut to be
% calculated based on whether main beam is cut by mirror timestream
% mask.  We want to move away from hand cuts as well, by making
% auto cuts more robust, but I'll leave the functionality in for now.
%
% INPUTS (should all be passed in through compopt)
%
%   expt:            'bicep2','keck','bicep3'
%   year:            bicep2: 2012 only
%                    keck: 2012,...,2018
%                    bicep3: 2016,...,2018
%
% OPTIONAL INPUTS
%
%   mapcoord:        coordinate of windowed maps to load
%                    'xpyp' (default), 'azel_ap'
%   demodtype:       'square' (default), 'choplet'
%   size:            windowed map size in deg (8 default)
%   suffix:          optional string to add to map name, '' default
%   bmdir:           directory from which to load beammaps
%                    (default beammaps/maps_windowed)
%   mapstoload:      vector of beam map run numbers to look at, i.e. 1:10
%                    (default is everything in bmrunlist.csv)
%   cutdir:          directory in which to save the cut structure 
%                    (default beammaps/cuts)
%   cutfile:         name of cut file to save to
%                    (default ffbm_cuts_year)
%   cutlist:         list of things to cut on, i.e.
%                    {'mirror','nanfit','peak'} (default all)
%   pklim:           limits for peak fit cuts 
%                    (default is experiment/frequency dependent) 
%   siglim:          limits for sigma cuts 
%                    (default is experiment/frequency dependent) 
%   elliplim:        limits for ellipticity
%                    (can be experiment/frequency dependent)
%   pointlim:        upper limit for abs pointing in xp/yp coords (deg)
%                    (default is 0.05, not freq dependent)
%
%   For multifrequency limits, send in a cell with one entry per frequency
%   in ascending order, i.e. for 100, 150, 220 GHz (as in p.band):
%     pklim{1} = [100 200]; pklim{2} = [100 300]; pklim{3} = [100 200]; 
%
% OUTPUTS
%
%   Saves a "cut" struct to cutdir of the same length as the number of beam
%   map runs processed, with corresponding fields for each type of cut
%   To generate a "standard" cut struct, run
%   >> compopt.expt = 'keck'; compopt.year = 2014;
%   >> ffbm_makecuts(compopt)
%   Note that as of now, if you run this and already have a structure with
%   hand cuts, they will be overwritten!!!

% Parse compopt
expt = compopt.expt;
if strcmp('expt','bicep2') 
  year = 2012;
else
  year = compopt.year;
end
if ~isfield(compopt,'mapcoord');
  mapcoord = 'xpyp';
else
  mapcoord = compopt.mapcoord;
end
if ~isfield(compopt,'demodtype');
  demodtype = 'square';
else
  demodtype = compopt.demodtype;
end
if ~isfield(compopt,'size');
  sizestr = '8deg';
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
      'sigma','ellip','abspoint','hand'};
  cutlist = compopt.cutlist;
else
  cutlist = compopt.cutlist;
end

% Set default fit limits here
[p ind] = get_array_info([num2str(year) '0201']);

for ii = 1:length(cutlist)
  switch cutlist{ii}
    case 'peak' % Eventually want this to be chopperdiam-based
      if ~isfield(compopt,'pklim')
	switch expt
          case 'bicep2'
            compopt.pklim{1} = [50 700];
	  case 'keck'
	    bands = nonzeros(unique(p.band));
	    for ii = 1:length(bands)
	      switch bands(ii)
		case 100
		  compopt.pklim{ii} = [10 500];
		case 150
                  switch year
                    case {2014,2015}
                      compopt.pklim{ii} = [25 500];
                    case 2016
                      compopt.pklim{ii} = [25 1100];
                  end
                case 210
                  compopt.pklim{ii} = [10 800];
		case 220
                  switch year
                    case 2015
                      compopt.pklim{ii} = [10 500];
                    case {2016,2017}
                      compopt.pklim{ii} = [10 800];
                  end
                case 270
                  compopt.pklim{ii} = [30 800];
	      end
	    end
	  case 'bicep3'
            switch year
              case {2016,2017}
                compopt.pklim{1} = [30 250];
            end
	end
	pklim = compopt.pklim;
      else
	pklim = compopt.pklim;
      end
    case 'sigma'
      if ~isfield(compopt,'siglim')
	switch expt
          case 'bicep2'
            compopt.siglim{1} = [0.15 0.3];
	  case 'keck'
	    bands = nonzeros(unique(p.band));
	    for ii = 1:length(bands)
	      switch bands(ii)
		case 100
		  compopt.siglim{ii} = [0.22 0.37];
		case 150
		  compopt.siglim{ii} = [0.15 0.3];
		case {210,220}
		  compopt.siglim{ii} = [0.1 0.28];
                case 270
                  compopt.siglim{ii} = [0.08 0.25];
	      end
	    end
	  case 'bicep3'
	    compopt.siglim{1} = [0.13 0.21]; 
	end
	siglim = compopt.siglim;
      else
	siglim = compopt.siglim;
      end
    case 'ellip'
      % Not frequency dependent now, but keep functionality
      if ~isfield(compopt,'elliplim')
	switch expt
          case 'bicep2'
            compopt.elliplim{1} = [0 0.2];
	  case 'keck'
	    bands = nonzeros(unique(p.band));
	    for ii = 1:length(bands)
	      switch bands(ii)
		case 100
		  compopt.elliplim{ii} = [0 0.2];
		case 150
		  compopt.elliplim{ii} = [0 0.2];
		case {210,220,270}
		  compopt.elliplim{ii} = [0 0.2];
	      end
	    end
	  case 'bicep3'
	    compopt.elliplim{1} = [0 0.2]; 
	end
	elliplim = compopt.elliplim;
      else
	elliplim = compopt.elliplim;
      end
    case 'abspoint'
      if ~isfield(compopt,'pointlim')
        compopt.pointlim = 0.05;
        pointlim = compopt.pointlim;
      else
        pointlim = compopt.pointlim;
      end
  end
end

% If file already exists, check that you really want to overwrite it
savename = [cutdir '/' cutfile '.mat'];
if exist(savename,'file')
  s = input([savename ' exists.  Overwrite?'],'s');
  if strcmp(s,'n') || strcmp(s,'no')
    return
  end
end

% Load up each beam map.  Only use windowed map!
for ii = 1:length(bm.number)
  
  filename = bm.filename{ii};
  disp(['Loading map ' num2str(ii) ': ' filename]);
  wmapname = [bmdir,'/mapwin_',filename,'_',demodtype,...
	'_',mapcoord,'_',sizestr,suffix];
  w = load(wmapname);
  w = w.bm;
  
  % Calculate beam parameters 
  [sig e_p e_c] = ffbm_findingparam(w.A,p,w.dk,1,expt);
  
  % Loop over cuts
  % General strategy: set most things to 1 (cut) and only allow beams that
  % pass each cut to go to 0 
  for jj = 1:length(cutlist)
    switch cutlist{jj}
      case 'mirror'
	switch expt
	  case 'keck'
	    %cut(ii).mirror = mirrorcut(p,bm.mirror{ii},bm.dk(ii));
            cut(ii).mirror = ~w.onmirror;
	  case {'bicep2','bicep3'}
	    cut(ii).mirror = zeros(size(p.gcp)); % Everything good here!
	end
      case 'notlight'
	tmp = ones(size(p.gcp));
	tmp(ind.l) = 0;
	cut(ii).notlight = tmp;
      case 'nanfit'
	tmp = ones(size(p.gcp));
	tmp(~isnan(w.A(1,:))) = 0;
	cut(ii).nanfit = tmp;
      case 'peak'
	%cut(ii).peak = peakcut(p,w.A(1,:),pklim);
	cut(ii).peak = cut_freq_dependent(p,w.A(1,:),pklim);
      case 'sigma'
	cut(ii).sigma = cut_freq_dependent(p,sig,siglim);
      case 'ellip'
	ellip = sqrt(e_p.^2 + e_c.^2);
	cut(ii).ellip = cut_freq_dependent(p,ellip,elliplim);
      case 'abspoint'
        tmp = ones(size(p.gcp));
        passcond = (abs(w.A(2,:)) < pointlim) & (abs(w.A(3,:)) < pointlim);
        tmp(passcond) = 0;
        cut(ii).abspoint = tmp;
      case 'hand'
	% Take a different tack with hand cuts here: assume none, then use
	% ffbm_makehandcuts  
	cut(ii).hand = zeros(size(p.gcp));
    end
  end

end

% Prep output struct
cuts.cuts = cut;
cuts.compopt = compopt;
cuts.bm = bm;

% Save
save(savename,'cuts');
disp(['Saved ' savename])

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cutarr = mirrorcut(p,mirrorpos,dk)

cutarr = ones(length(p.gcp),1);

% Get the general mask for each of the 10 potential rx positions
[mask100 mask150] = ffbm_maskMaker(mirrorpos);

% We don't care about beams that exist but which are distorted - only keep
% good ones!
mask100(mask100 == 2) = 1;
mask150(mask150 == 2) = 1;

% Multifrequency...great
for ii = 1:length(p.gcp)
  rxpos = dk_rx_pos_ffbm(p.rx(ii),dk);
  switch p.band(ii)
    case 100
      cutarr(ii) = mask100(rxpos,p.gcp(ii) + 1);
    case {150,210,220,270}
      cutarr(ii) = mask150(rxpos,p.gcp(ii) + 1);
  end
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cutarr = cut_freq_dependent(p,cutvals,cutlim)
% General cutting scheme when we have low/high limits which may depend on
% frequency
% cutvals should be a vector of the same length as p.gcp
% cutlim should have the same number of fields as unique bands

% Be ruthless!  1 means cut
cutarr = ones(size(p.gcp));

% Multifrequency...great
for ii = 1:length(p.gcp)

  bands = nonzeros(unique(p.band));
  % Choose appropriate limits for band
  limind = find(bands == p.band(ii));
  if ~isempty(limind) % There are some 0s in p.band...
    low = cutlim{limind}(1);
    high = cutlim{limind}(2);
    if cutvals(ii) > low & cutvals(ii) < high
      cutarr(ii) = 0;
    end
  end

end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cutarr = peakcut(p,fitamps,pklim)

cutarr = ones(size(p.gcp));

% Multifrequency...great
for ii = 1:length(p.gcp)

  bands = nonzeros(unique(p.band));
  % Choose appropriate pklim for band
  pklimind = find(bands == p.band(ii));
  if ~isempty(pklimind)
    pklow = pklim{pklimind}(1);
    pkhigh = pklim{pklimind}(2);
    if fitamps(ii) > pklow & fitamps(ii) < pkhigh
      cutarr(ii) = 0;
    end
  end

end

return
