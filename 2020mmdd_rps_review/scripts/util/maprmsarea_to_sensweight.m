function [totalsens,surveyweight]=maprmsarea_to_sensweight(rms,area)
% [totalsens,surveyweight]=maprmsarea_to_sensweight(rms,area)
%
% We often quote rms Q/U map noise and effective area
% In the B2 instrument paper it defines "total sensitivity" as
% s=rms/sqrt(2*area) - which should be equivalent to NET/sqrt(t).
% The Keck paper goes on to define "survey weight" as 1/s^2
%
% e.g.
% Keck paper
% maprmsarea_to_sensweight(74,390) = [2.7,0.1424]

totalsens=rms/sqrt(2*area); % in units of temperature

surveyweight=1/totalsens^2; % in units of 1/T^2

return
