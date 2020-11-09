function [sampratio,samprate]=get_sampratio_rate(d)
% [sampratio,samprate]=get_sampratio_rate(d)
%
% in B2 partial loadcurves are recorded at a different rate - Walt has
% provided code which comes up with the correct numbers

% Number of fasts per slow (constant)
sampratio=length(d.t)/length(d.ts);

% Construct on-source cut (and fast-sample version)
% to select scans only: no bracketing cals.
c_onsrc=bitand(d.array.frame.features,2^1)>0;
c_onsrc_f=repmat(c_onsrc(:)',sampratio,1);
c_onsrc_f=c_onsrc_f(:);

% Note that the sample rate changes during the course of a scanset!
% The partial load curves are archived at the full rate of 100 Hz,
% but the el nods and scans are downsampled to 20 Hz as they are recorded.

% Calculate sample rate using *only* scans on source.
dtlist=diff(d.t(c_onsrc_f))*3600*24;
% Further throw out distant outliers
dtmedian=nanmedian(dtlist);
dtlist=dtlist(dtlist>=dtmedian/10 & dtlist<=dtmedian*3);
samprate=1/mean(dtlist);                         % fast in Hz

%disp(['using samprate=' num2str(samprate) ' Hz']);

return
