function varargout=caxis_map(map,varmap)
% clim=caxis_map(map,varmap)
%
% Chooses a color axis range based on the provided map. If no assignment is
% made (i.e. nargout==0), then caxis(clim) is called automatically.
%
% If varmap is provided, it is used to limit which region of the map that is
% examined --- only pixels that surpass 25% of the peak weight are
% considered.
%
% If varmap is not provided, a box 60-by-40 pixels in size at the center of
% the map is used. (This is compatible with how reduc_plotcomap_pager() has
% historically computed the free color scale.)
%

  trunc = NaN(size(map,1), size(map,2));

  if ~exist('varmap','var') || isempty(varmap)
    center = fix(size(map) / 2);
    lox = center(1) - 30;
    hix = center(1) + 30;
    loy = center(2) - 20;
    hiy = center(2) + 20;
    varmap = NaN(size(map));
    varmap(lox:hix,loy:hiy,:) = 1.0;
  end

  if size(varmap,3) ~= size(map,3)
    if size(varmap,3) == 1
      varmap = repmat(varmap, [1, 1, size(map,3)]);
    else
      error('map and varmap have incompatible shapes')
    end
  end

  for ii=1:size(map,3)
    % Peak normalize the apodization mask
    apmask = (1./varmap(:,:,ii)) .* min(cvec(varmap(:,:,ii)));
    % Based on the region which reaches at least 25% of the maximum depth
    trunc(apmask>=0.25) = 1.0;
    npixap = sum(isfinite(trunc(:)));

    map(:,:,ii) = squeeze(map(:,:,ii)) .* trunc;
    if false
      % Traditional way, but may allow outlier "hot" pixel to survive
      clim(ii) = max(abs(cvec(map(:,:,ii))));
    else
      % Use quantiles with high fraction to avoid hot pixels
      sigma3 = 1 - 2*(1 - normcdf(3, 0, 1));
      clim(ii) = quantile(abs(cvec(map(:,:,ii))), sigma3);
    end
    if ~isfinite(clim(ii))
      clim(ii) = NaN;
    end
  end
  clim = max(clim);
  if ~isfinite(clim)
    clim = 1;
  end

  clim = clim .* [-1,1];

  if nargout == 0
    caxis(clim);
  else
    varargout{1} = clim;
  end
end
