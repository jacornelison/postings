function hc=dualexp_tfunc(tc,f)
% hc=dualexp_tfunc(tc,f)
%
% Generate dual exponential transfer func evalulated at freq f
% including effect of QUaD butterworth electronics filters

% analog butterworth electronics filters
[bb1,ab1]=butter(6,20,'s');
[bb2,ab2]=butter(2,30,'s');
% calc butterworth trans funcs for analog filters
hb1=freqs(bb1,ab1,f);
hb2=freqs(bb2,ab2,f);
hb=hb1.*hb2;

% exp filter in analog domain
be1=[0,1];
ae1=[tc(1),1];
% cal exp trans func for analog filter
he1=freqs(be1,ae1,f*2*pi); % 2pi required because tc in seconds
 
% if using 2nd tc
if(length(tc)>1)
  % and if 2nd tc has weight greater than zero for this det
  if(tc(3)>0)
    be2=[0,1];
    ae2=[tc(2),1];
    he2=freqs(be2,ae2,f*2*pi);
    
    hc=he1.*hb.*(1-tc(3))+he2.*hb.*tc(3);
  end
else
  hc=he1.*hb;
end

return
