function ffbm_poweroutsider(bmopt)
%
% function ffbm_poweroutsider(bmopt)
%
% TSG 2017-09-29
%
% Take composite maps and, for each detector, calculate the 
% total integrated power outside radius r as a function of r.
%
% This should be used with large composites (8 deg).
%
% INPUTS (passed through struct bmopt):
%  
%  .expt = 'keck','bicep3'
%  .year = (Keck) 2012-2017
%            (B3) 2016,2017
%
% OPTIONAL INPUTS:
%
%   .compositedir:  directory from which to load composite maps
%                   (default beammaps/maps_composite)
%   .compositefile: file in compositedir to load
%                   (default ffbm_year_all_8deg_xpyp.mat)
%   .savedir:       directory where results are saved
%                   (default beammaps/integrated_power)
%   
% OUTPUTS:
%   
%   Saves per-det integrated power outside r vs r to specified output
%   directory.
%

tic;

% Parse input
expt = bmopt.expt;
year = bmopt.year;
if ~isfield(bmopt,'compositedir')
  bmopt.compositedir = 'beammaps/maps_composite';
  compositedir = bmopt.compositedir;
else
  compositedir = bmopt.compositedir;
end
if ~isfield(bmopt,'compositefile')
  bmopt.compositefile = ['ffbm_' num2str(year) '_all_8deg_xpyp'];
  compositefile = bmopt.compositefile;
else
  compositefile = bmopt.compositefile;
end
if ~isfield(bmopt,'savedir')
  bmopt.savedir = 'beammaps/integrated_power';
  savedir = bmopt.savedir;
else
  savedir = bmopt.savedir;
end

[p ind] = get_array_info([num2str(year) '0301']);
p = rmfield(p,'expt');

% Load composites
w = load([compositedir '/' compositefile]);
ad = w.ad;
map = w.map;

% Radius vector, degrees
R = ad.t_r*180/pi;
rmax = min(abs([ad.t_val_deg{1}(1),ad.t_val_deg{1}(end),...
                ad.t_val_deg{2}(1),ad.t_val_deg{2}(end)]));
r = linspace(0,rmax,round(ad.N_pix(1)/2) + 1);

% Store integrated power vectors here
ibp = NaN(length(map),length(r));

% Get integrated power for each det
for idet = 1:length(map)
  p_tot = nansum(map(idet).T(:));
  for ii = 1:length(r)
    if ~isempty(map(idet).T)
      z = map(idet).T(R >= r(ii));
      ibp(idet,ii) = nansum(z) / p_tot;
    end
  end
end

% Save
savename = [savedir '/ffbm_' num2str(year) '_all_ibp'];
save(savename,'r','ibp')
disp(['Saved ' savename ])
toc

return