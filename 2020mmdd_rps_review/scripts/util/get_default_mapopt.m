function mapopt=get_default_mapopt(mapopt)
% mapopt=get_default_mapopt(mapopt)
%
% Fills in the mapopt structure with defaults for use with reduc_makepairmaps.
%
% If mapopt.sernum is not specified on input, mapopt.sernum and mapopt.realpairmapset
% are not defined on output.
% 
% ex. mapopt.filt='p0';
%     mapopt.sernum='12050011';
%     mapopt=get_default_mapopt(mapopt);
%               OR
%     mapopt=get_default_mapopt;

if(~exist('mapopt','var'))
  mapopt=[];
end
if(~isfield(mapopt,'type'))
    mapopt.type='bicep';
end
if(~isfield(mapopt,'proj'))
  mapopt.proj='radec';
end
if(~isfield(mapopt,'binning'))
  mapopt.binning='ra';
end
if(~isfield(mapopt,'filt'))
  mapopt.filt='p3';
end
if(~isfield(mapopt,'filttype'))
  mapopt.filttype='normal';
end
if(~isfield(mapopt,'gs'))
  mapopt.gs=1; % default is now gs on
end
if(~isfield(mapopt,'weight'))
  mapopt.weight=3; % reciprocal variance over scanset post filtering
end
if(~isfield(mapopt,'beamcen'))
  mapopt.beamcen='obs';
end
if(~isfield(mapopt,'chi'))
  mapopt.chi='obs';
end
if(~isfield(mapopt,'epsilon'))
  mapopt.epsilon='obs';
end
if(~isfield(mapopt,'usepsm'))
  mapopt.usepsm=0;
end
if(~isfield(mapopt,'resrelgain'))
  mapopt.resrelgain=0;
end
if(~isfield(mapopt,'sernum'))
  disp(['No serial number (mapopt.sernum) specified. '...
           'mapopt.sernum and mapopt.realpairmapset will not be defined.'])
elseif ~isfield(mapopt,'realpairmapset')
  % if not specified look for matching real set
  mapopt.realpairmapset=['pairmaps/' mapopt.sernum(1:4),'/real'];
  if ~exist(mapopt.realpairmapset,'dir') & isempty(findstr(mapopt.sernum,'real'))
    warning([mapopt.realpairmapset ' does not exist. If you ' ...
             'are using weights from the real data, select an existing pairmap ' ...
             'set using mapopt.realpairmapset']);
  end
end
if(~isfield(mapopt,'cut'))
  % if not specified get standard cuts
  mapopt.cut=get_default_round1_cuts;
end

% Sparse pairmaps by default. This is very desirable in the coadd stage.
if(~isfield(mapopt,'acsparse'))
  mapopt.acsparse=true;
end

% Do not pack pairmaps by default. Size decrease is negligible when saving pairmaps as
% -v7 but the slowdown in the coadd stage is significant over many tags.
if(~isfield(mapopt,'acpack'))
  mapopt.acpack=false;
end

if(~isfield(mapopt,'deproj'))
  mapopt.deproj=true;
end

if(~isfield(mapopt,'deproj_map'))
  %mapopt.deproj_map='input_maps/wmap_derivs/wmap_band_iqumap_r9_7yr_V_37_v4.fits';
  mapopt.deproj_map='input_maps/planck/planck_derivs_nopix/synfast_deproj_143_nominal_B2.fits';
end

if ~isfield(mapopt,'deproj_pix_offset')
  mapopt.deproj_pix_offset=0;
end

if ~isfield(mapopt,'mapdarksquids')
  mapopt.mapdarksquids=0;
end

if ~isfield(mapopt,'mapfpntds')
  mapopt.mapfpntds=0;
end

if ~isfield(mapopt,'update')
  mapopt.update=false;
end

if ~isfield(mapopt,'abmap')
  mapopt.abmap = false;
end

if ~isfield(mapopt,'resmult') || isempty(mapopt.resmult)
  mapopt.resmult = 1;
end

return
