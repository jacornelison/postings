function d=timeconst_fitfunc_matmin(tc,func,bi,y,labtc)

tc

if(exist('labtc','var'))
  % we are doing a scaled dual tc fit against lab values
  labtc(1:2)=labtc(1:2)*tc;
  tc=labtc;
end

y=apply_tc(tc,y);

% pull out blip sections
yf=y(bi.f);
yb=y(bi.b);

% sum of square of vector of differences
d=sum((yf-yb).^2)

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y=apply_tc(tc,y)

% deconv timeconst from full scan group
fs.sf=1; fs.ef=length(y);
y=deconv_scans(y,fs,1,tc);
d.v=y; d=lowpassfilt(d,fs,{'v'}); y=d.v;

return
