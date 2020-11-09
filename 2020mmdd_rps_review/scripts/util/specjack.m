function r=specjack(r,a,b,c,simtype)
% r=specjack(r,a,b,c,simtype)
%  
% difference spectra and sims
%  
% simtype = 'sim' or 0 : use sim  field
%          = 'simr' or 1: use simr field if available
%          = 'simd'     : use simd field if available
%          = '...'      : string indicating whatever fields are available

if(~exist('simtype','var'))
  simtype='sim';
else
  if isnumeric(simtype) 
    switch simtype
      case 1
        simtype='simr';
      case 0
        simtype='sim';
    end
  end
end

% fall back
if ~isfield(r(a),simtype)
  display(['''',simtype,''' field not available for spectral jack, use ''sim'' instead'])
  simtype='sim';
end

% this selects the 6 spectra which are available both in auto and cross,
% removes the additional three that come in the cross (ET,BT,BE):
r(a).real=r(a).real(:,1:6); r(a).(simtype)=r(a).(simtype)(:,1:6,:); r(a).noi=r(a).noi(:,1:6,:);
r(b).real=r(b).real(:,1:6); r(b).(simtype)=r(b).(simtype)(:,1:6,:); r(b).noi=r(b).noi(:,1:6,:);
r(c).real=r(c).real(:,1:6); r(c).(simtype)=r(c).(simtype)(:,1:6,:); r(c).noi=r(c).noi(:,1:6,:);

% difference real, s+n and noise sims
ro(a).real=r(a).real-r(c).real; ro(a).sim=r(a).(simtype)-r(c).(simtype); ro(a).noi=r(a).noi-r(c).noi;
ro(b).real=r(a).real-r(b).real; ro(b).sim=r(a).(simtype)-r(b).(simtype); ro(b).noi=r(a).noi-r(b).noi;
ro(c).real=r(c).real-r(b).real; ro(c).sim=r(c).(simtype)-r(b).(simtype); ro(c).noi=r(c).noi-r(b).noi;

% copy ell over
ro(a).l=r(a).l; ro(b).l=r(b).l; ro(c).l=r(c).l;

% replace r with differenced one
r=ro;

return
