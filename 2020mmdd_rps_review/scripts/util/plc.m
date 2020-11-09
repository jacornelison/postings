function [lc dat]=plc(d,p,lc,ch)
% [lc dat]=plc(d,p,lc,ch)
%
% Reduce partial load curve cals. The load curves are
% analyzed in two ways: using the entire load curve
% including the normal region, and using only the
% derivative at the bias point. The second method is
% less vulnerable to jumps and other discontinuities
% in the data, but it assumes constant P at the bias
% point.
%
% Outputs:
%   lc is the lc structure including tweaked start and end times
%   and the following fields from PLC analysis:
%
%   x(n,m,:) where n=1:2 for leading/trailing PLC
%            and m=1:nchan is the detector index.
%
%   Calibration / gain factor dP_J/dfb in W/fbu:
%    d(n,m,1)=dpj/dfb best estimate, W/fbu
%    d(n,m,2)=dpj/dfb from bias only, no load curve analysis
%    d(n,m,3)=dpj/dfb from derivatives near bias point assuming constant P
%    d(n,m,4)=dpj/dfb from full analysis including normal region
%
%   Detector resistances in Ohms:
%    r(n,m,1)=rdet at bias point from derivatives
%    r(n,m,2)=rdet at bias point from full analysis including normal region
%    r(n,m,3)=rnorm from full analysis including normal region
%
%   Chisq of idet relative to simple model
%    chisq(n,m,1)=chisq of idet for full analysis including normal region
%
%   Old-style output for backwards compatibility
%    g(n,m,1)=dpj/dfb best estimate, same as d(n,m,1)
%    g(n,m,2)=pj(end)
%    g(n,m,3)=rnorm
%    g(n,m,4)=rdet(end)
%
%   r_shunt value from calibration file
%    r_sh(m)=p.r_sh
%
%   dat(n,m) has fields pj, vdet, idet, rdet, and bias from full analysis
%   including normal region.
%

disp('partial loadcurve...')

% Get calibration values as provided by get_array_info
lc.r_sh=p.r_sh;

% standard biases for each column
% use these to determine when partial load curve has finished
std_bias=mode(double(d.mce0.tes.bias));

% Extract bias steps from user word
all_bias=calc_bias(d);

% Tweak the ends of plc a bit due to feature timing issue.
[sampratio,samprate]=get_sampratio_rate(d);
lc=tweak_plc(d,lc,sampratio);

% Fix the sizes of the outputs
n=length(lc.s); m=size(d.mce0.data.fb,2);
lc.g=NaN*zeros(n,m,4);
dat(n,m).idet=[]; dat(n,m).vdet=[]; dat(n,m).pj=[]; dat(n,m).rdet=[]; dat(n,m).bias=[];
lc.chisq=NaN*zeros(n,m,1);
lc.r=NaN*zeros(n,m,3);
lc.d=NaN*zeros(n,m,4);

% Number of chans per column -- different for BICEP3!
col_per_mce=size(d.mce0.tes.bias,2)/size(d.mce0.frame.status,2);
if floor(col_per_mce) ~= col_per_mce
  error(['Number of columns per MCE results in non-integer value.']);
end
colnum = 1 + p.mce_col + col_per_mce * p.mce;
colnum(~isfinite(colnum)) = 1;

% don't loop over channels with std_bias==0
chuse=ch(std_bias(colnum(ch))~=0);

for ii=1:length(lc.sf)

  sf=lc.sf(ii);
  ef=lc.ef(ii);

  for j=chuse

    % Get TES feedback, num flux jumps for channels in this MCE
    fb=d.mce0.data.fb(sf:ef,j);
    fj=double(d.mce0.data.numfj(sf:ef,j));
    bias=all_bias(sf:ef,1+p.mce(j));

    % Cut out portion between steps
    gapcut=(fb==0 & fj==0) | bias==0;
    dropcut=gapcut | isnan(fb) | isnan(bias);

    % Cut out portion below operating bias for this column
    biascut = (bias < std_bias(colnum(j)));

    % Cut out any dropped portion
    fb=fb(~dropcut & ~biascut);
    fj=fj(~dropcut & ~biascut);
    bias=bias(~dropcut & ~biascut);

    % First few samples are already at combine off,
    % but may ringing from large step in bias because
    % of MCE filter.
    fb=fb(3:end);
    fj=fj(3:end);
    bias=bias(3:end);

    % deglitching sometimes NaN's out so much data
    % reduction of the plc is not possible
    if isempty(fb) || mean(fb)==0
      continue
    end

    [tmpdat, dpjdfb, r, chisq]=fit_loadcurve(fb,bias,j,fj,p);
    lc.g(ii,j,1)=dpjdfb(1);
    lc.g(ii,j,2)=tmpdat.pj(end);
    lc.g(ii,j,3)=r(3);
    lc.g(ii,j,4)=r(2);
    lc.chisq(ii,j,1)=chisq;
    lc.r(ii,j,1:3)=r;
    lc.d(ii,j,1:4)=dpjdfb;
    dat(ii,j).pj=tmpdat.pj;
    dat(ii,j).vdet=tmpdat.vdet;
    dat(ii,j).idet=tmpdat.idet;
    dat(ii,j).rdet=tmpdat.rdet;
    dat(ii,j).bias=tmpdat.ibias;
    % catch
    %  error(['load curve analysis failed on channel ' num2str(j)]);
    %  g(ii,j,1:5)=NaN;
    %  dat(ii,j).pj=NaN;
    %  dat(ii,j).vdet=NaN;
    %  dat(ii,j).idet=NaN;
    %  dat(ii,j).rdet=NaN;
    %  dat(ii,j).bias=NaN;
    % end
  end
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TES bias steps are stored in user_word, which usually means something different.
% user_word > 20000 is a time stamp, not a valid bias.
% Assume bias steps down in uniform steps, get rid of anything that doesn't do this.
function bias=calc_bias(d)

bias=zeros(size(d.mce0.header.user_word));
for imce=(1:size(bias,2))-1
  % Get raw user word
  uw=double(d.mce0.header.user_word(:,imce+1));

  % For each uw value, find next distinct and nonzero value, i.e.
  % following bias step or userword that's not a bias step
  lastuw=uw(1);
  lastidx=1;
  next_bias=zeros(size(uw));
  while ~isempty(lastidx)
    idx=find(uw((lastidx+1):end)~=0 & uw((lastidx+1):end)~=lastuw,1,'first');
    if isempty(idx)
      break
    end
    idx=idx+lastidx;
    next_bias(lastidx:(idx-1))=uw(idx);
    lastuw=uw(idx);
    lastidx=idx;
  end

  % Get rid of values that are definitely not bias values
  uw(uw>20000)=0;

  % Find step size and lowest bias step
  bvals=unique(uw(isfinite(uw) & uw>0));
  stepsize=mode(diff(bvals));
  stepmin=min(bvals);

  % No bias?
  if isempty(stepmin)
    continue
  end

  % Good bias values are those that are one of the steps,
  % and have bias = next_bias+stepsize
  cgood=false(size(uw));
  cgood(uw==stepmin)=true;
  cgood(mod(uw-stepmin,stepsize)==0 & next_bias+stepsize==uw)=true;

  % Select good bias values
  bias(cgood,imce+1)=uw(cgood);
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fit a straight line to normal portion, remove offset, calculate load curve 
% parameters
function [dat, dpjdfbu, r, chisq]=fit_loadcurve(fb,bias,idx,fj,pcalib)
% INPUTS:
% fb_data = fb_data, correctly scaled for data_mode and filter gain
% tbias_data = applied detector bias data
% idx = channel index
% pcalib = calibration structure
%
% OUTPUTS:
% dat = structure with samples per tbias point:
%   pj = vector of detector joule power at each tbias point (W)
%   vdet = vector of detetor voltage at each tbias point (V)
%   idet = vector of detector current at each tbias point (A)
%   rdet = vector of detector resistance at each tbias point (Ohms)
% dpjdfbu = detector dP_J/dfb at bias point (W/fbu) (four values: best guess, from bias only, from derivs, from full)
% r = detector resistances (Ohms) (three values: at bias point from derivs; at bias point from full; rnorm)

% compute calibration factors from these circuit parameters
rsh = pcalib.r_sh(idx);
fbu_calfac = (pcalib.v_fb1_max(idx)/((2^pcalib.bits_fb1(idx))*(pcalib.r_fb1(idx)+pcalib.r_wire(idx))*pcalib.m_fb1(idx)));
bias_calfac= (pcalib.v_b_max(idx)/((2^pcalib.bits_bias(idx))*(pcalib.r_bias(idx))));

% apply calibration factors
input_current = -1 * fbu_calfac * fb;
tbias_current = bias_calfac * bias;
dat.ibias = tbias_current;

% dpjdfbu without any correction
dpjdfbu = NaN * [1 1 1 1];
if length(tbias_current)>=1
  dpjdfbu(1:2)=tbias_current(end) * rsh * fbu_calfac;
end

% initialize r list
r = NaN * [1 1 1];

% chisq not yet filled in
chisq = NaN;

% Analysis from derivs at bias point
if length(tbias_current)>=5
  nuse=3;
  ub=unique(tbias_current);
  if length(ub)>=nuse
    cc=tbias_current<=ub(min(nuse,length(ub)));
    p=polyfit(tbias_current(cc),input_current(cc),1);
    r(1) = rsh * (-1/p(1) - 1);  
    corrfac = (r(1) - rsh) / (r(1) + rsh);
    if ~(corrfac<0.1 || corrfac>1.0 || r(1)<0)
      dpjdfbu([1,3]) = dpjdfbu(2) * corrfac;
    end
  end
end

% Identify normal/bias branch crossover
% Don't accept a minimum that's within 20 samples of the start; this is just
% ringing / instability when the TES bias switches from normal to high
tmp_input_current=input_current;
tmp_input_current((1:length(input_current))<=20)=Inf;
[min_val, min_idx]=min(tmp_input_current);

% Select normal branch as 50% of samples from min_idx to max bias current
norm_start=1;
normfrac=.50;
norm_end=round(normfrac*min_idx);
if isempty(norm_end)
  norm_end=1;
end

% If not enough samples to do anything useful
if norm_end-norm_start<=5 || length(unique(tbias_current(norm_start:norm_end)))<2
  dat.pj=NaN*ones(max(1,length(bias)),1);
  dat.vdet=dat.pj;
  dat.idet=input_current;
  dat.rdet=dat.pj;
  return
end

% Subtract off y-intercept of linear fit to normal region
p=polyfit(tbias_current(norm_start:norm_end),input_current(norm_start:norm_end),1);
input_current=input_current-p(2);
rnorm=rsh*((1/p(1))-1);

% Calculate resistance, voltage, Joule power
dat.idet=input_current;
dat.rdet=rsh*(tbias_current./dat.idet-1);
dat.vdet=dat.rdet.*input_current;
dat.pj=dat.vdet.*input_current;

% Fill in r values
r(2) = dat.rdet(end);
r(3) = rnorm;

% dP_J/dI_s correction factor
corrfac = (r(2)-rsh)/(r(2)+rsh);

% Sanity checks
if ~(corrfac<0.1 || corrfac>1.0 || r(2)<0 || r(2)>rnorm)
  % Correccted dP_J/dI_s at bias point
  dpjdfbu([1,4]) = dpjdfbu(2) * corrfac;
end

% Calculate chisq relative to simple model
model_idet=max(dat.vdet./rnorm,dat.pj(end)./dat.vdet);
sig=nanstd(diff(dat.idet));
chisq=nanstd(dat.idet-model_idet)/sig;

return

