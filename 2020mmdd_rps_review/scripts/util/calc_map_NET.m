function NET=calc_map_NET(map,m,coaddopt,time_pertag,plots,bin)
% function calc_map_NET(map,m,coaddopt,time_pertag,plots,bin)
%
% Written to analyze the NET (per detector and overall)
%
% Inputs:
%   m,map,coaddopt are the standard saved quantities
%   time_pertag (default = 2060 s)
%       This is the time per tag spent observing.
%       Used to turn the NET per detector into NET per map
%       The conversion is NET(map)=NET(det)/sqrt(ndet)=NET(det)/sqrt(integrated time/time)
%             = NET(det)/sqrt(pitime/(time_pertag*ntags))
%   plots (default = 1)
%       Will plot a histogram to screen if 1
%       Always check the histogram. This is a known issue in K2017 estimate.
%       http://bicep.rc.fas.harvard.edu/keck/analysis_logbook/analysis/20170809_K2017_prelim_timestream_NET/
%   bin (default = 3*95th percentile of the absolute values)
%       Specifies the edge of the histogram plotted. 
%       (It can be a single value for all maps, or it can be an array specifying the edge for each map.)
%       Also, points beyond that will not be used to calculate the NET since it's sensitive to outlier points
%
% Output:
%   NET - structure with Q,U,perchan_Q,perchan_U,and effective npix, 
%         each arrays of size(map,1)
%
% Example:
%   load('maps/1350/real_b_filtp3_weight3_gs_dp1102_jack21.mat')
%   map=make_map(ac,m,coaddopt); 
%   map=cal_coadd_maps(map,2900);
%   NET = calc_map_NET(map,m,coaddopt)
%

if ~exist('time_pertag','var') || isempty(time_pertag)
  % currently taken from a Keck pairmap
  time_pertags=2.06e+03;
end
if ~exist('plots','var')
  plots=1;
end
if ~exist('bin','var')
  % choice of binning for the histograms.  is also used to cut out crazy pixels
  bin = NaN;
end
if numel(map)>1 && numel(bin)==1
  bin = repmat(bin, size(map));
end

% save the pitime since it's removed when jackknifed
for ii=1:size(map,1)
    time(ii).Pitime=(map(ii,1).Pitime+map(ii,2).Pitime);
end

% for the effective pixels and turning it into a real number
ntags=length(coaddopt.tags);

% jackknife the map
map=jackknife_map(map);

% loop through and calculate the NET
for ii=1:size(map,1)

  % in calculations, only use the data with map.Qvar and map.Uvar are <sqrt(2)*median
  % otherwise the NETs are driven by the edge pixels.  
  % weighting by time should significantly down weight them, but it's not enough.
  % that is partly because the time is total polarization time, but the Q maps don't have the full
  % polarization coverage.  this is not as much an issue for >=2013 scan strategy
  qlimit=1.4*nanmedian(rvec(map(ii).Qvar));
  ulimit=1.4*nanmedian(rvec(map(ii).Uvar));
  Q=map(ii).Q(map(ii).Qvar<qlimit);
  U=map(ii).U(map(ii).Uvar<ulimit);
  usedQPitime=time(ii).Pitime(map(ii).Qvar<qlimit);
  usedUPitime=time(ii).Pitime(map(ii).Uvar<ulimit);

  % Calculate and plot the noise from the Q and U jackknifes- per chan
  dQ=Q.*sqrt(usedQPitime);
  dU=U.*sqrt(usedUPitime);
  if ~isfinite(bin(ii))
    bin(ii) = 3 * prctile(abs([cvec(dQ);cvec(dU)]), 95);
  end
  NET.perchan_Q(ii)=nanstd(dQ(abs(dQ)<bin(ii)));
  NET.perchan_U(ii)=nanstd(dU(abs(dU)<bin(ii)));

  if plots 
    qlegtext=['Q_{noise}=' num2str(NET.perchan_Q(ii),3)];
    ulegtext=['U_{noise}=' num2str(NET.perchan_U(ii),3)];
    bins = linspace(-bin(ii), bin(ii), 101);
    figure; plot_net_histogram({dQ,dU},bins,{qlegtext,'Q fit',ulegtext,'U fit'}, ...
      ['Q/U Scan Jackknife*sqrt(time) Histograms']);
  end

  % Save the effective number of pixels
  NET.effective_npix(ii)=nansum(rvec(time(ii).Pitime))./(ntags*time_pertags);

  % Convert the NET per channel into NET over the sub map
  NET.Q(ii)=NET.perchan_Q(ii)./sqrt(NET.effective_npix(ii)*2);
  NET.U(ii)=NET.perchan_U(ii)./sqrt(NET.effective_npix(ii)*2);

end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plot_net_histogram(map,bins,legtext,label)

colorder={'k','b','c'};
colorder_fit={'r','m','g'};

box on; hold on
for i=1:numel(map)
  d=rvec(map{i});
  [n,x]=hist(d,bins);
  stairs(x,n,colorder{i})
  [mu sig]=normfit(d);
  y=normpdf(x,mu,sig)*sum(n)*(bins(2)-bins(1));
  plot(x,y,colorder_fit{i})
end

title(label);
xlabel('\muK_{cmb}');
ylabel('pixel count');
axis tight
set(gca,'XGrid','on','YGrid','on');
l1=legend(legtext);
set(l1,'Box','off');
hold off

return
