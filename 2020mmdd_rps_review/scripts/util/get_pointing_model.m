function pmi=get_pointing_model(mjd,mirror,d,fs)
% pmi=interp_pointing_model(mjd,mirror,d,fs)
%
% interpolate the pointing model to mjd
%
% e.g. pm=get_pointing_model(mean(d.t(d.t~=0)))

% fetch the pointing model data
pm=ParameterRead('aux_data/pointingmodels/pointingmodel_complete_v1.csv');

if ~exist('mirror','var')
  mirror=0;
end

% get field names
fn=fieldnames(pm);

% if mirror is present, just get the dummy pointing model:
if mirror
  if ~exist('fs','var')
    pmi.az_zero=double(d.antenna0.tracker.encoder_off(1,1))/3.6e6;
    pmi.el_zero=double(d.antenna0.tracker.encoder_off(1,2))/3.6e6;
  else
    pmi.az_zero=double(d.antenna0.tracker.encoder_off(fs.s(1),1))/3.6e6;
    pmi.el_zero=double(d.antenna0.tracker.encoder_off(fs.s(1),2))/3.6e6;
  end
  pmi.el_tilt=0;
  pmi.az_tilt_lat=0;
  pmi.az_tilt_ha=0;
  pmi.bracket=0;
  return
end

% warn user if we don't have bracketing star-pointing schedules
bracket = min(pm.mjd) < mjd & mjd < max(pm.mjd);
if ~bracket
  disp(['warning: using a date not bracketed by a boresight pointing' ...
        ' model measurement. re-run tod when model is updated.'])
end

% for each numerical field perform interp to required mjd
pmi = struct();
if mjd < min(pm.mjd)
  fbidx = 1;
elseif max(pm.mjd) < mjd;
  fbidx = numel(pm.mjd);
end
for i=1:length(fn)
  if isreal(pm.(fn{i}))
    if ~bracket
      pmi.(fn{i}) = pm.(fn{i})(fbidx);
    else
      pmi.(fn{i}) = interp1(pm.mjd, pm.(fn{i}), mjd, 'linear');
    end
  end
end
pmi.bracket = bracket;

return

