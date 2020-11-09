function libsphere_perf()
  libsphere()
  compare_distance();
  compare_cosdistance();
  compare_bearing();
  compare_polbearings();
end


function compare_distance()
  % Assumes BK pipeline is available

  nside = 512;
  pix = cvec(0:(nside2npix(nside)-1));
  ipix = ang2pix(nside, 0, -57.5);

  % Pixel locations:
  %
  % Unit vectors, for libsphere's distance
  rvec = pix2vec(nside, pix);
  % Lat-Lon pairs, for Matlab's distance
  [lat,lon] = pix2ang(nside, pix);
  clear pix;

  if true
    % Convert to pointer buffers for speed when using libsphere
    rvecptr = libpointer('doublePtr', rvec);
    distptr = libpointer('doublePtr', zeros(size(rvec,1), 1));
  else
    % This way is allowed, but we sacrifice speed due to Matlab creating
    % copies that protect the inputs.
    rvecptr = rvec;
    distptr = zeros(size(rvec,1), 1);
  end

  function d=libsphere()
    [ret,d] = calllib('libsphere', 'distance', size(rvec,1), distptr, rvecptr, ipix);
    if isa(d, 'lib.pointer')
      d = d.Value;
    end
  end
  function d=mlnative()
    d = distance(lon, lat, lon(ipix+1), lat(ipix+1));
  end

  if verLessThan('matlab', '8.2') % R2013b
    libsphere_time = naivetimeit(@libsphere, 1);
    mlnative_time  = naivetimeit(@mlnative, 1);
  else
    libsphere_time = timeit(@libsphere, 1);
    mlnative_time  = timeit(@mlnative, 1);
  end

  ratio = mlnative_time / libsphere_time;
  fprintf(1, ['distance:\n' ...
              '    libsphere: %0.3g sec (%0.2f× speedup)\n' ...
              '    ML native: %0.3g sec\n' ...
              '\n'], ...
      libsphere_time, ratio, mlnative_time);

  if false
    d1 = libsphere();
    d2 = mlnative();
    clf(); plot(d1 - d2); title('delta between outputs');
    keyboard
  end
end

function compare_cosdistance()
  % Assumes BK pipeline is available

  nside = 512;
  pix = cvec(0:(nside2npix(nside)-1));
  ipix = ang2pix(nside, 0, -57.5);

  % Pixel locations:
  %
  % Unit vectors, for libsphere's distance
  rvec = pix2vec(nside, pix);
  % Lat-Lon pairs, for Matlab's distance
  [lat,lon] = pix2ang(nside, pix);
  clear pix;

  if true
    % Convert to pointer buffers for speed when using libsphere
    rvecptr = libpointer('doublePtr', rvec);
    distptr = libpointer('doublePtr', zeros(size(rvec,1), 1));
  else
    % This way is allowed, but we sacrifice speed due to Matlab creating
    % copies that protect the inputs.
    rvecptr = rvec;
    distptr = zeros(size(rvec,1), 1);
  end

  function d=libsphere()
    [ret,d] = calllib('libsphere', 'cosdistance', size(rvec,1), ...
        distptr, rvecptr, ipix);
    if isa(d, 'lib.pointer')
      d = d.Value;
    end
  end
  function d=mlnative()
    d = distance(lon, lat, lon(ipix+1), lat(ipix+1));
    d = cosd(d);
  end

  if verLessThan('matlab', '8.2') % R2013b
    libsphere_time = naivetimeit(@libsphere, 1);
    mlnative_time  = naivetimeit(@mlnative, 1);
  else
    libsphere_time = timeit(@libsphere, 1);
    mlnative_time  = timeit(@mlnative, 1);
  end

  ratio = mlnative_time / libsphere_time;
  fprintf(1, ['cosdistance:\n' ...
              '    libsphere: %0.3g sec (%0.2f× speedup)\n' ...
              '    ML native: %0.3g sec\n' ...
              '\n'], ...
      libsphere_time, ratio, mlnative_time);

  if false
    d1 = libsphere();
    d2 = mlnative();
    clf(); plot(d1 - d2); title('delta between outputs');
    keyboard
  end
end

function compare_bearing()
  % Assumes BK pipeline is available

  nside = 512;
  pix = cvec(0:(nside2npix(nside)-1));
  ipix = ang2pix(nside, 0, -57.5);

  % Pixel locations:
  %
  % Unit vectors, for libsphere's distance
  rvec = pix2vec(nside, pix);
  % Lat-Lon pairs, for Matlab's distance
  [lat,lon] = pix2ang(nside, pix);
  clear pix;

  if true
    % Convert to pointer buffers for speed when using libsphere
    rvecptr = libpointer('doublePtr', rvec);
    aptr = libpointer('doublePtr', zeros(size(rvec,1), 1));
  else
    % This way is allowed, but we sacrifice speed due to Matlab creating
    % copies that protect the inputs.
    rvecptr = rvec;
    aptr = zeros(size(rvec,1), 1);
  end

  function a=libsphere()
    [ret,a] = calllib('libsphere', 'bearing', size(rvec,1), ...
        aptr, rvecptr, ipix);
    if isa(a, 'lib.pointer')
      a = a.Value;
    end
  end
  function a=mlnative()
    a = azimuth(lon(ipix+1), lat(ipix+1), lon, lat);
    a = unwrap(deg2rad(a));
  end

  if verLessThan('matlab', '8.2') % R2013b
    libsphere_time = naivetimeit(@libsphere, 1);
    mlnative_time  = naivetimeit(@mlnative, 1);
  else
    libsphere_time = timeit(@libsphere, 1);
    mlnative_time  = timeit(@mlnative, 1);
  end

  ratio = mlnative_time / libsphere_time;
  fprintf(1, ['bearing:\n' ...
              '    libsphere: %0.3g sec (%0.2f× speedup)\n' ...
              '    ML native: %0.3g sec\n' ...
              '\n'], ...
      libsphere_time, ratio, mlnative_time);

  if false
    a1 = libsphere();
    a2 = mlnative();
    clf();
    plot(a1 - a2);
      ylim([-1,1] .* 1e-13); title('delta between angle outputs');
    keyboard
  end
end

function compare_polbearings()
  % Assumes BK pipeline is available

  nside = 512;
  pix = cvec(0:(nside2npix(nside)-1));
  ipix = ang2pix(nside, 0, -57.5);

  % Pixel locations:
  %
  % Unit vectors, for libsphere's distance
  rvec = pix2vec(nside, pix);
  % Lat-Lon pairs, for Matlab's distance
  [lat,lon] = pix2ang(nside, pix);
  clear pix;

  if true
    % Convert to pointer buffers for speed when using libsphere
    rvecptr = libpointer('doublePtr', rvec);
    cijptr = libpointer('doublePtr', zeros(size(rvec,1), 1));
    sijptr = libpointer('doublePtr', zeros(size(rvec,1), 1));
    cjiptr = libpointer('doublePtr', zeros(size(rvec,1), 1));
    sjiptr = libpointer('doublePtr', zeros(size(rvec,1), 1));
  else
    % This way is allowed, but we sacrifice speed due to Matlab creating
    % copies that protect the inputs.
    rvecptr = rvec;
    cijptr = zeros(size(rvec,1), 1);
    sijptr = zeros(size(rvec,1), 1);
    cjiptr = zeros(size(rvec,1), 1);
    sjiptr = zeros(size(rvec,1), 1);
  end

  function [cij,sij,cji,sji]=libsphere()
    [ret,cij,sij,cji,sji] = calllib('libsphere', 'polbearings', ...
        size(rvec,1), cijptr, sijptr, cjiptr, sjiptr, rvecptr, ipix);
    if isa(cij, 'lib.pointer')
      cij = cij.Value;
      sij = sij.Value;
      cji = cji.Value;
      sji = sji.Value;
    end
  end
  function [cij,sij,cji,sji]=mlnative()
    a = azimuth(lon(ipix+1), lat(ipix+1), lon, lat);
    cij = cosd(2*a);
    sij = sind(2*a);
    a = azimuth(lon, lat, lon(ipix+1), lat(ipix+1));
    cji = cosd(2*a);
    sji = sind(2*a);
  end

  if verLessThan('matlab', '8.2') % R2013b
    libsphere_time = naivetimeit(@libsphere, 4);
    mlnative_time  = naivetimeit(@mlnative, 4);
  else
    libsphere_time = timeit(@libsphere, 4);
    mlnative_time  = timeit(@mlnative, 4);
  end

  ratio = mlnative_time / libsphere_time;
  fprintf(1, ['polbearings:\n' ...
              '    libsphere: %0.3g sec (%0.2f× speedup)\n' ...
              '    ML native: %0.3g sec\n' ...
              '\n'], ...
      libsphere_time, ratio, mlnative_time);

  if false
    [cij1,sij1,cji1,sji1] = libsphere();
    [cij2,sij2,cji2,sji2] = mlnative();
    clf();
    subplot(4,1,1); plot(cij1 - cij2);
      ylim([-1,1] .* 1e-13); title('delta between cij outputs');
    subplot(4,1,2); plot(sij1 - sij2);
      ylim([-1,1] .* 1e-13); title('delta between sij outputs');
    subplot(4,1,3); plot(cji1 - cji2);
      ylim([-1,1] .* 1e-13); title('delta between cji outputs');
    subplot(4,1,4); plot(sji1 - sji2);
      ylim([-1,1] .* 1e-13); title('delta between sji outputs');
    keyboard
  end
end

function t=naivetimeit(func, nout)
  if ~exist('nout', 'var') || isempty(nout)
    nout = 1;
  end

  nrep = 10;
  t = tic();
  for ii=1:nrep
    [outargs{1:nout}] = func();
  end
  t = toc(t) / nrep;
end

