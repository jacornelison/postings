function params = ffbm_evalbeamparam(compopt)
% function params = ffbm_evalbeamparam(compopt)
%
% Summarize beam parameters from a beam param file. 
% Outputs a struct with median individual/differential beam parameters, and
% for each the standard deviation or width of the n_beammaps measurements
%
% INPUTS (should all be sent in with compopt)
%   
%   expt:            'b2','keck','b3'
%   year:            2011/2012/'both' for bicep2   
%                    2012/2013/2014/2015 for keck
%                    2015                for bicep3
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
%
% OUTPUTS
%
%   params:          struct containing summary of beam parameters

% Parse compopt
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
if ~isfield(compopt,'rxin')
  compopt.rxin = [];
  rxin = compopt.rxin;
else
  rxin = compopt.rxin;
end

% Load beamparam file
filename = [beamparamdir '/' beamparamfile];
load(filename);
beampar = beamparam;

numFiles = size(beampar.sig,2);
  
% Get array info
if strcmp(expt,'b3')
  expt = 'bicep3';
end
if chflag
  flags = get_default_chflags(expt,year);
  [p ind] = get_array_info(sprintf('%d0201',year),...
      'obs','obs','obs','obs',flags);
else
  [p ind] = get_array_info(sprintf('%d0201',year));
end

% Cut down to an rx if requested
if ~isempty(rxin)
  cutind = (p.rx == rxin);
  p = rmfield(p,'expt');
  p = structcut(p,cutind);
  ind = make_ind(p);

  if size(beampar.sig,1) ~= size(ind.e,2)
    beampar.sig = beampar.sig(cutind,:);
    beampar.e_p = beampar.e_p(cutind,:);
    beampar.e_c = beampar.e_c(cutind,:);
    beampar.x = beampar.x(cutind,:);
    beampar.y = beampar.y(cutind,:);
  end
end

sig = beampar.sig;
e_p = beampar.e_p;
e_c = beampar.e_c;
x = beampar.x;
y = beampar.y;

% Process the beam parameters
params.sig_med = nanmedian(sig(ind.e,:),2);
params.sig_std = nanstd(sig(ind.e,:),[],2);
params.e_p_med = nanmedian(e_p(ind.e,:),2);
params.e_p_std = nanstd(e_p(ind.e,:),[],2);
params.e_c_med = nanmedian(e_c(ind.e,:),2);
params.e_c_std = nanstd(e_c(ind.e,:),[],2);
params.x_med = nanmedian(x(ind.e,:),2);
params.y_med = nanmedian(y(ind.e,:),2);

% Get distribution widths instead of stds
for ii = 1:length(ind.e)
  in = sig(ind.e(ii),:);
  % Make sure there's at least 2 measurements for an uncertainty
  if length(in(~isnan(in))) >= 2
    tmp = percentile(in(~isnan(in)),onesig);
    params.sig_wid(ind.e(ii),1) = (tmp(2) - tmp(1))/2;
  else
    params.sig_wid(ind.e(ii),1) = NaN;
  end
  in = e_p(ind.e(ii),:);
  if length(in(~isnan(in))) >= 2
    tmp = percentile(in(~isnan(in)),onesig);
    params.e_p_wid(ind.e(ii),1) = (tmp(2) - tmp(1))/2;
  else
    params.e_p_wid(ind.e(ii),1) = NaN;
  end
  in = e_c(ind.e(ii),:);
  if length(in(~isnan(in))) >= 2
    tmp = percentile(in(~isnan(in)),onesig);
    params.e_c_wid(ind.e(ii),1) = (tmp(2) - tmp(1))/2;
  else
    params.e_c_wid(ind.e(ii),1) = NaN;
  end
end

% Differential parameters

% A/B offset
dx(ind.e,numFiles) = NaN;
dy(ind.e,numFiles) = NaN;
dx(ind.a,:) = x(ind.a,:) - x(ind.b,:);
dy(ind.a,:) = y(ind.a,:) - y(ind.b,:);

params.aboffset_x = nanmedian(dx,2);
params.aboffset_y = nanmedian(dy,2);

params.aboffset_err_x(ind.e,1) = NaN;
params.aboffset_err_y(ind.e,1) = NaN;
params.aboffset_wid_x(ind.e,1) = NaN;
params.aboffset_wid_y(ind.e,1) = NaN;

params.aboffset_err_x(ind.a,1) = nanstd(dx(ind.a,:),[],2);
params.aboffset_err_x(ind.b,1) = params.aboffset_err_x(ind.a);
params.aboffset_err_y(ind.a,1) = nanstd(dy(ind.a,:),[],2);
params.aboffset_err_y(ind.b,1) = params.aboffset_err_y(ind.a);

% A/B offset: get distribution widths instead of stds
for ii = 1:length(ind.a)
  in = dx(ind.a(ii),:);
  if length(in(~isnan(in))) >= 2
    tmp = percentile(in(~isnan(in)),onesig);
    params.aboffset_wid_x(ind.a(ii),1) = (tmp(2) - tmp(1))/2;
    params.aboffset_wid_x(ind.b(ii),1) = (tmp(2) - tmp(1))/2;
  else
    params.aboffset_wid_x(ind.a(ii),1) = NaN;
    params.aboffset_wid_x(ind.b(ii),1) = NaN;
  end
  in = dy(ind.a(ii),:);
  if length(in(~isnan(in))) >= 2
    tmp = percentile(in(~isnan(in)),onesig);
    params.aboffset_wid_y(ind.a(ii),1) = (tmp(2) - tmp(1))/2;
    params.aboffset_wid_y(ind.b(ii),1) = (tmp(2) - tmp(1))/2;
  else
    params.aboffset_wid_y(ind.a(ii),1) = NaN;
    params.aboffset_wid_y(ind.b(ii),1) = NaN;
  end
end

% Differential beamwidth
dsig = NaN(length(ind.e),numFiles);
dsig(ind.a,:) = sig(ind.a,:) - sig(ind.b,:);
params.dsig_med = nanmedian(dsig,2);
params.dsig_std = nanstd(dsig,[],2);

params.dsig_wid(ind.e,1) = NaN;
params.de_c_wid(ind.e,1) = NaN;
params.de_p_wid(ind.e,1) = NaN;

% DSig: widths
for ii = 1:length(ind.a)
  in = dsig(ind.a(ii),:);
  if length(in(~isnan(in))) >= 2
    tmp = percentile(in(~isnan(in)),onesig);
    params.dsig_wid(ind.a(ii),1) = (tmp(2) - tmp(1))/2;
    params.dsig_wid(ind.b(ii),1) = (tmp(2) - tmp(1))/2;
  else
    params.dsig_wid(ind.a(ii),1) = NaN;
    params.dsig_wid(ind.b(ii),1) = NaN;
  end
end

% Differential ellipticity
de_p = NaN(length(ind.e),numFiles);
de_p(ind.a,:) = e_p(ind.a,:) - e_p(ind.b,:);
params.de_p_med = nanmedian(de_p,2);
params.de_p_std = nanstd(de_p,[],2);
de_c = NaN(length(ind.e),numFiles);
de_c(ind.a,:) = e_c(ind.a,:) - e_c(ind.b,:);
params.de_c_med = nanmedian(de_c,2);
params.de_c_std = nanstd(de_c,[],2);

% DE: widths
for ii = 1:length(ind.a)
  in = de_p(ind.a(ii),:);
  if length(in(~isnan(in))) >= 2
    tmp = percentile(in(~isnan(in)),onesig);
    params.de_p_wid(ind.a(ii),1) = (tmp(2) - tmp(1))/2;
    params.de_p_wid(ind.b(ii),1) = (tmp(2) - tmp(1))/2;
  else
    params.de_p_wid(ind.a(ii),1) = NaN;
    params.de_p_wid(ind.b(ii),1) = NaN;
  end
  in = de_c(ind.a(ii),:);
  if length(in(~isnan(in))) >= 2
    tmp = percentile(in(~isnan(in)),onesig);
    params.de_c_wid(ind.a(ii),1) = (tmp(2) - tmp(1))/2;
    params.de_c_wid(ind.b(ii),1) = (tmp(2) - tmp(1))/2;
  else
    params.de_c_wid(ind.a(ii),1) = NaN;
    params.de_c_wid(ind.b(ii),1) = NaN;
  end
end

return