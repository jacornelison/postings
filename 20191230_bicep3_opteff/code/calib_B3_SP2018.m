% calib_B3_SP2018.m
%
% Sets calibration values for BICEP3 based on run04 at pole
%
% SNK 20151228

function calib=calib_B3_SP2018()


disp('Using BICEP3 (based on Run4) calibrations')

numCols = 32;

calib.R_SH = 0.0040; % Ohms (at 4K, different chip)

calib.R_WIRE = 151.0*ones(1,32); % ohms include in R_FB and R_BIAS

% Bias
calib.R_BIAS = 475.00*ones(1,32);
calib.BITS_BIAS = 15; % number of bits in TES bias DAC
calib.V_BIAS_MAX = 5.000*ones(1,32);

% FB
calib.R_FB = 2114*ones(1,32);
calib.V_FB_MAX = 0.972*ones(1,32);
calib.M_FB = 25*ones(1,32);

calib.BITS_FB = 14; % number of bits in FB DAC

% Bias ADC bins => bias current (Amps)
calib.BIAS_CAL = (calib.V_BIAS_MAX/(2^calib.BITS_BIAS)) ./ ...
  (calib.R_BIAS + calib.R_WIRE);
% SQ1 FB DAC bins => TES current (Amps)
calib.FB_CAL = (calib.V_FB_MAX/(2^calib.BITS_FB)) ./ ...
  (calib.R_FB + calib.R_WIRE) ./ calib.M_FB ;
