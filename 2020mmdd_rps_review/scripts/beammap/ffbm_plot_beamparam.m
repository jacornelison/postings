function ffbm_plot_beamparam(compopt)
% ffbm_plot_beamparam(compopt)
%
% Load up beam parameter files and plot them
% Also print out the mean and uncertainties
%
% Steps to get here:
% 1. Make/save the cuts
%     ffbm_makecuts(compopt);
% 2. Get/save parameters with cuts applied
%    ffbm_makebeamparamfiles(compopt);
% 
% Note that this file uses all measurements in the beamparam files,
% whereas ffbm_plot_beamparam_tile just uses the final saved .csv file
%
% INPUTS (should all be sent in with compopt)
%
%  expt:             'keck','bicep3'
%  year:             keck: 2014/2015/2016
%                    bicep3:   2015
%
% OPTIONAL INPUTS
%
%   beamparamdir:    directory from which to load beamparam file
%                    (default beammaps/beamparams)
%   beamparamfile:   beamparam file to load
%                    (default beamparams_year)
%   onesig:          percentile range to get distribution widths
%                    (default [0.16 0.84], i.e. one sigma)
%   chflag:          0 = none (default)
%                    1 = load up and apply defaults
%   rxin:            rx to cut down get_array_info 
%                    (empty or [] for full array)
%   wid:             0 (use std)
%                    1 to use onesig above (default)
%   makepng:         1 to make nice .pngs
%                    0 (default) to just make .eps
%
% OUTPUTS
%
%    A bunch of .eps/png files are saved in the local directory 
  
expt = compopt.expt;
year = compopt.year;
if ~isfield(compopt,'beamparamdir')
  compopt.beamparamdir = 'beammaps/beamparams';
  beamparamdir = compopt.beamparamdir;
else
  beamparamdir = compopt.beamparamdir;
end
if ~isfield(compopt,'beamparamfile')
  compopt.beamparamfile = ['beamparams_' num2str(year)];
  beamparamfile = compopt.beamparamfile;
else
  beamparamfile = compopt.beamparamfile;
end
if ~isfield(compopt,'onesig')
  compopt.onesig = [0.16 0.84];
  onesig = compopt.onesig;
else
  onesig = compopt.onesig;
end
if ~isfield(compopt,'chflag')
  compopt.chflag = 0;
  chflag = compopt.chflag;
else
  chflag = compopt.chflag;
end
% Get array info
if chflag
  chflags = get_default_chflags(expt,year);
  [p ind] = get_array_info(sprintf('%d0201',year),...
      'obs','obs','obs','obs',chflags);
else
  [p ind] = get_array_info(sprintf('%d0201',year));
end

p = rmfield(p,'expt'); % for structcut to work

if ~isfield(compopt,'rxin')
  compopt.rxin = unique(p.rx);
  rxin = compopt.rxin;
else
  if isempty(compopt.rxin)
    compopt.rxin = unique(p.rx);
    rxin = compopt.rxin;
  else
    rxin = compopt.rxin;
  end
end
if ~isfield(compopt,'wid')
  compopt.wid = 1;
  wid = compopt.wid;
else
  wid = compopt.wid;
end
if ~isfield(compopt,'makepng')
  compopt.makepng = 0;
  makepng = compopt.makepng;
else
  makepng = compopt.makepng;
end

for ii = 1:length(rxin)
  
  rx = rxin(ii);
  % Get beam means, std, and widths over all runs
  % compopt.beamparamdir/beamparamfile determine which to use
  compopt.rxin = rx;
  cutind = (p.rx == rx);
  p_rx = structcut(p,cutind);
  ind_rx = make_ind(p_rx);
  % Since we set compopt.rxin, param_rx will be cut to match p_rx/ind_rx
  param_rx = ffbm_evalbeamparam(compopt);

  disp(' ')
  set(0,'defaultlinelinewidth',1.0);
  set(0,'DefaultAxesFontSize',14);

  % Make all the plots
  make_beamsummary_plot(ind_rx.e,param_rx.sig_med,...
      param_rx.sig_wid,'Sigma (deg)',rx,expt,year,makepng)
  make_beamsummary_plot(ind_rx.e,param_rx.e_p_med,...
      param_rx.e_p_wid,'Ellipticity plus',rx,expt,year,makepng)
  make_beamsummary_plot(ind_rx.e,param_rx.e_c_med,...
      param_rx.e_c_wid,'Ellipticity cross',rx,expt,year,makepng)
    
  make_beamsummary_plot(ind_rx.e,param_rx.dsig_med,...
      param_rx.dsig_wid,'Diff Sigma (deg)',rx,expt,year,makepng)
  make_beamsummary_plot(ind_rx.e,param_rx.aboffset_x,...
      param_rx.aboffset_wid_x,'Diff Pointing X',rx,expt,year,makepng)
  make_beamsummary_plot(ind_rx.e,param_rx.aboffset_y,...
      param_rx.aboffset_wid_y,'Diff Pointing Y',rx,expt,year,makepng)
  make_beamsummary_plot(ind_rx.e,param_rx.de_c_med,...
      param_rx.de_c_wid,'Diff Ellipticity cross',rx,expt,year,makepng)
  make_beamsummary_plot(ind_rx.e,param_rx.de_p_med,...
      param_rx.de_p_wid,'Diff Ellipticity plus',rx,expt,year,makepng)

  % Summarize statistics
  
  % Take stats mean/std over detectors
  index1 = ind_rx.rgl;
  stats.sig = nanmean(param_rx.sig_med(index1));
  stats.p = nanmean(param_rx.e_p_med(index1));
  stats.c = nanmean(param_rx.e_c_med(index1));
  
  disp(sprintf('### %s %d Rx%d parameters',expt,year,rx))
  if wid
    tmp = percentile(param_rx.sig_med(index1),onesig,1);
    stats.sig_scatter = (tmp(2) - tmp(1))/2;
    stats.sig_meas_unc = nanmedian(param_rx.sig_wid(index1));
    tmp = percentile(param_rx.e_p_med(index1),onesig,1);
    stats.p_scatter = (tmp(2) - tmp(1))/2;
    stats.p_meas_unc = nanmedian(param_rx.e_p_wid(index1));
    tmp = percentile(param_rx.e_c_med(index1),onesig,1);
    stats.c_scatter = (tmp(2) - tmp(1))/2;
    stats.c_meas_unc = nanmedian(param_rx.e_c_wid(index1));
  else
    stats.sig_scatter = nanstd(param_rx.sig_med(index1));
    stats.sig_meas_unc = nanmedian(param_rx.sig_std(index1));
    stats.p_scatter = nanstd(param_rx.e_p_med(index1));
    stats.p_meas_unc = nanmedian(param_rx.e_p_std(index1));
    stats.c_scatter = nanstd(param_rx.e_c_med(index1));
    stats.c_meas_unc = nanmedian(param_rx.e_c_std(index1));
  end
  
  % Differential parameters
  index2 = ind_rx.rgla;
  stats.dx = nanmean(param_rx.aboffset_x(index2));
  stats.dy = nanmean(param_rx.aboffset_y(index2));
  stats.dsig = nanmean(param_rx.dsig_med(index2));
  stats.dp = nanmean(param_rx.de_p_med(index2));
  stats.dc = nanmean(param_rx.de_c_med(index2));
  
  if wid
    tmp = percentile(param_rx.aboffset_x(index2),onesig,1);
    stats.dx_scatter = (tmp(2) - tmp(1))/2;
    stats.dx_meas_unc = nanmedian(param_rx.aboffset_wid_x(index2));
    tmp = percentile(param_rx.aboffset_y(index2),onesig,1);
    stats.dy_scatter = (tmp(2) - tmp(1))/2;
    stats.dy_meas_unc = nanmedian(param_rx.aboffset_wid_y(index2));
    tmp = percentile(param_rx.dsig_med(index2),onesig,1);
    stats.dsig_scatter = (tmp(2) - tmp(1))/2;
    stats.dsig_meas_unc = nanmedian(param_rx.dsig_wid(index2));
    tmp = percentile(param_rx.de_p_med(index2),onesig,1);
    stats.dp_scatter = (tmp(2) - tmp(1))/2;
    stats.dp_meas_unc = nanmedian(param_rx.de_p_wid(index2));
    tmp = percentile(param_rx.de_c_med(index2),onesig,1);
    stats.dc_scatter = (tmp(2) - tmp(1))/2;
    stats.dc_meas_unc = nanmedian(param_rx.de_c_wid(index2));
  else
    stats.dx_scatter = nanstd(param_rx.aboffset_x(index2));
    stats.dx_meas_unc = nanmedian(param_rx.aboffset_err_x(index2));
    stats.dy_scatter = nanstd(param_rx.aboffset_y(index2));
    stats.dy_meas_unc = nanmedian(param_rx.aboffset_err_y(index2));
    stats.dsig_scatter = nanstd(param_rx.dsig_med(index2));
    stats.dsig_meas_unc = nanmedian(param_rx.dsig_std(index2));
    stats.dp_scatter = nanstd(param_rx.de_p_med(index2));
    stats.dp_meas_unc = nanmedian(param_rx.de_p_std(index2));
    stats.dc_scatter = nanstd(param_rx.de_c_med(index2));
    stats.dc_meas_unc = nanmedian(param_rx.de_c_std(index2));
  end
  
  disp(['Beam width = ' num2str(stats.sig,'%0.3f') ...
        ' +/- ' num2str(stats.sig_scatter,'%0.3f') ...
        ' +/- ' num2str(stats.sig_meas_unc,'%0.3f')]);
  disp(['Plus ellip = ' num2str(stats.p,'%0.3f') ...
        ' +/- ' num2str(stats.p_scatter,'%0.3f') ...
        ' +/- ' num2str(stats.p_meas_unc,'%0.3f')]);
  disp(['Cross ellip = ' num2str(stats.c,'%0.3f') ...
        ' +/- ' num2str(stats.c_scatter,'%0.3f') ...
        ' +/- ' num2str(stats.c_meas_unc,'%0.3f')]);
  disp(['Diff X, dx = ' num2str(stats.dx*60,'%0.2f') ...
        ' +/- ' num2str(stats.dx_scatter*60,'%0.2f') ...
        ' +/- ' num2str(stats.dx_meas_unc*60,'%0.2f')]);
  disp(['Diff Y, dy = ' num2str(stats.dy*60,'%0.2f') ...
        ' +/- ' num2str(stats.dy_scatter*60,'%0.2f') ...
        ' +/- ' num2str(stats.dy_meas_unc*60,'%0.2f')]);
  disp(['Diff bw, dsig = ' num2str(stats.dsig,'%0.1d') ...
        ' +/- ' num2str(stats.dsig_scatter,'%0.3f') ...
        ' +/- ' num2str(stats.dsig_meas_unc,'%0.3f')]);
  disp(['Diff plus, dp = ' num2str(stats.dp,'%0.3f') ...
        ' +/- ' num2str(stats.dp_scatter,'%0.3f') ...
        ' +/- ' num2str(stats.dp_meas_unc,'%0.3f')]);
  disp(['Diff cross, dc = ' num2str(stats.dc,'%0.3f') ...
        ' +/- ' num2str(stats.dc_scatter,'%0.3f') ...
        ' +/- ' num2str(stats.dc_meas_unc,'%0.3f')]);
  
end

return % Main function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function make_beamsummary_plot(x,y,yerr,ylab,rx,expt,year,makepng)
 
if strfind(ylab,'Sigma') == 1
  switch expt
    case 'keck'
      switch year
        case 2015
          switch rx
            case {0,2}
              lim = [0.25 0.35];
            case {1,3}
              lim = [0.12 0.18];
            case 4
              lim = [0.19 0.23];
          end
        case 2016
          switch rx
            case {0,1,2,3}
              lim = [0.12 0.18];
            case 4
              lim = [0.19 0.23];
          end
        case 2017
          switch rx
            case {0,1,2,3}
              lim = [0.12 0.18];
            case 4
              lim = [0.10 0.16];
          end
        otherwise
          m = nanmedian(y);
          e = nanstd(y);
          lim = [m-10*e m+10*e];
      end
    case 'bicep3'
      m = nanmedian(y);
      e = nanstd(y);
      lim = [m-10*e m+10*e];
  end
  histlim = [lim(1):0.002:lim(2)];
elseif strfind(ylab,'Ellipticity') == 1
  lim = [-0.2 0.2];
  histlim = [lim(1):0.01:lim(2)];
elseif  strfind(ylab,'Diff Sigma') == 1
  lim = [-.025 .025];
  histlim = [lim(1):0.001:lim(2)];
elseif ~isempty(strfind(ylab,'Diff Pointing')) ||  ...
       ~isempty(strfind(ylab,'Diff Ellipticity'))
  lim = [-.1 .1];
  histlim = [lim(1):0.005:lim(2)];
end  

if strfind(ylab,'Diff')
  y(find(y == 0)) = NaN;
end

fh = figure(1); clf;
setwinsize(gcf,1200,400);
set(fh,'Visible','off')

subplot(121)
errorbar2(x,y,yerr,'X');
hold on;
if strcmp(expt,'keck') % Rx dividing lines
  plot([132 132],lim,'--k');
  plot([264 264],lim,'--k');
  plot([396 396],lim,'--k');
end
det_med = nanmedian(y);
plot([x(1) x(end)],[det_med det_med],'r--')
xlim([x(1)-1 x(end)+1]);
ylim(lim);
xlabel('gcp #');
ylabel(ylab);
grid();
tt = sprintf('%s %d Rx%d %s',expt,year,rx,ylab);
title(tt);

subplot(122)
hist(y,histlim)
xlabel(ylab);
ylabel('Number of hits');
xlim(lim);
grid();
tt = sprintf('%s %d Rx%d %s (Median = %2.3f)',expt,year,rx,ylab,det_med);
title(tt);
%suptitle(sprintf('%s %d',expt,year));
outname = strrep(ylab,' (deg)','');
outname = strrep(outname,' ','_');

savename = sprintf('%s_%d_%s_rx%d',expt,year,outname,rx);
print('-depsc2',savename);

if makepng
  system(['rm ' savename '.png']);
  system(['eps2png -B ' savename '.eps']);
end

return