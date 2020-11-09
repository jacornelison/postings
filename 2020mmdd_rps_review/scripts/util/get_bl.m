function [l,Bl]=get_bl(band,lin)
% [l,bl]=get_bl(band,lin)
%
% Get 1d B(l), e.g. beam transfer function that multiplies Fourier plane.
% B^2(l) is the beam transfer function for the power spectrum.
%
% band = 'nXXXX' where X are digits; HEALPix pixel windows
%      = numeric; returns Gaussian B(l), specified in arcmin FWHM
%      = any band name recognized by get_bl_filename()
%      = a file path to a FITS file
%
% lin  = Optional. If provided, beam profile is interpolated to the given
%        multipoles in `lin`.

  % Gaussian B_l, input is FWHM in arcmin
  if ~isstr(band)
    sigma = (band/60)*pi/180 / 2.355;
    Bl = exp(-0.5*lin.*(lin+1).*sigma^2);
    l=lin;
    return
  end

  % HEALPix pixel window functions for nNNNN where N are numbers, 'n0512'
  % for instance.
  if strncmp(band,'n',1) && all(isstrprop(band(2:end),'digit'))
    % Should these also be regested in get_bl_filename() ??
    filename=sprintf('/n/home03/sfliescher/Programs/healpix/data/pixel_window_%s.fits',band);
    data=fitsread(filename,'BinTable');
    x(:,2)=data{1};
    x(:,1)=0:size(x,1)-1;

  else
    % Not a pixel window function

    % Handle previous case where get_bl originally recognized 'bbns' instead
    % of 'B2bbns'.
    if strcmp(band,'bbns')
      band = 'B2bbns';
    end

    if exist_file(band)
      blfile = band;
    else
      % Get registered files
      blfile = get_bl_filename(band);
    end

    % Then load the data. Case-by-case for different types of inputs.

    % WMAP7 and WMAP9 have different number of rows of headers, so need to
    % handle the .txt files separately.
    if strncmp(band, 'wmap7', length('wmap7'))
      if ~iscell(blfile)
        blfile = {blfile};
      end
      x = cellfunc(@(f) dlmread(f,'',8,0), blfile);
      % Files contain matrices, so concatenate along 3rd dimension so we can
      % take a mean along dimension 3 to average the differencing assemblies'
      % responses.
      x = mean(cat(3, x{:}), 3);
    end
    if strncmp(band, 'wmap9', length('wmap9'))
      if ~iscell(blfile)
        blfile = {blfile};
      end
      x = cellfunc(@(f) dlmread(f,'',13,0), blfile);
      x = mean(cat(3, x{:}), 3);
    end

    % FITS files can be loaded with common handling
    if ~iscell(blfile) && strcmp(blfile(end-4:end), '.fits')
      finfo = fitsinfo(blfile);
      if any(ismember('Binary Table', finfo.Contents))
        data = fitsread(blfile, 'BinTable');
      else
        data = fitsread(blfile, 'Table');
      end
      x(:,2) = data{1};
      x(:,1) = 0:size(x,1)-1;
    end

    % Sanity check. In case a new band is defined in get_bl_filename but uses
    % a file type not recognized yet here, throw an error.
    if ~exist('x','var') || isempty(x)
      error('B_l file %s was not loaded in get_bl()', band)
    end
  end

  l=x(:,1);
  Bl=x(:,2);

  if exist('lin','var') && ~isempty(lin)
    Bl=interp1(l,Bl,lin);
    Bl(lin>max(l))=0;
    l=lin;
  end

end

