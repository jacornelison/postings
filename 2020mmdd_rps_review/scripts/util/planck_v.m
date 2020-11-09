function ret = planck_v(v,T,sm)
  
  % calculates the Planck spectra radiance as a function of frequency 
  %   and Temperature of the emitter
  % v is an array of frequencies in Hz
  % T is a single temperature in K
  % if sm = 0, the result is the generic Planck blackbody spectral radiance
  % if sm = 1, the result is the single moded expression. This is the default.
  
if(~exist('sm'))
  sm=[];
end
if(isempty(sm))
  sm=1;
end

[h,k,c]= constants() ;

% spectral radiance  (W m^-2 sr^-1 Hz^-1)
Bv = ( (2*h.*v.^3) ./ c^2 ) .* ( 1 ./ ( exp((v.*h)./(k.*T)) - 1 ) );

% Power per unit bandwidth absorbed by a detector is
% Qv =  A*Omega * Bv
% Antenna theorem states A*Omega = lambda^2 so
% Qv = lambda^2 * 2* h * v^3/c^2 *  1 ./ ( exp((v.*h)./(k.*T)) - 1 ) );
% Qv = h*v *  1 ./ ( exp((v.*h)./(k.*T)) - 1 ) );
% This latter expression is the single moded expression 
%for the power  per unit frequency absorbed by a single moded antenna, 1 polarization (that is 
% where the factor of 1/2 disappeared)

% single moded expression
Qv = (h.*v) .* ( 1 ./ ( exp((v.*h)./(k.*T)) - 1 ) );

if(sm)
  ret = Qv;
else
  ret = Bv;
end


