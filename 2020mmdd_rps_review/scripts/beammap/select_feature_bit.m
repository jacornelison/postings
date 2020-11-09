function [mask_slow, mask_fast] = select_feature_bit(data, fbit)
% [mask_slow, mask_fast] = select_feature_bit(data, fbit)
%
% Select periods of time series with feature bit set.
%
% [Arguments]
%   data       Data from arcfiles. Must include array.frame.features.
%   fbit       Value of the feature bit to select. Feature bit is binary, 
%              so setting bits f0+f1 corresponds to a value of 3.
% 
% [Returns]
%   mask_slow  Logical mask indicating when the feature bit is set for
%              slow registers.
%   mask_fast  Logical mask indicating when the feature bit is set for 
%              fast registers.

% 2014-02-06 CAB

% Slow register feature bit.
flag = data.array.frame.features;
mask_slow = (flag == fbit);

% Upsample to fast data rate.
sampratio = size(data.antenna0.time.utcfast, 1) / ...
    size(data.antenna0.time.utcslow, 1);
flag_fast = cvec(repmat(flag, 1, sampratio)');
mask_fast = (flag_fast == fbit);
