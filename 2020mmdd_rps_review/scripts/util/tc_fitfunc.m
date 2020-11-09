function y=filtmodel_fitfunc(par,x,labtc)
% y=filtmodel_fitfunc(par,x,labtc)
%
% gauss blip plus filter fit function used in reduc_findtc

% apply exponential impulse response filter to source blip model

% convert from scale param to dual tc/weight form
if(exist('labtc','var'))
  labtc(1:2)=labtc(1:2)*par(4);
  par=[par(1:3),labtc'];
end

par(4:end)

% setup gaussian in x - this is obviously only a crude model when
% scanning a real source
y=gauss(par(1:3),x);

if(par(4)>2e-9)

  % do exactly the same thing as deconv_scans since we are trying to
  % get appropriate values to use there

  % setup butterworth electronics filter
  % this is bilinear transform of analog filter
  [bb,ab]=butter(6,20/50);
  sosb=tf2sos(bb,ab);

  % exp filter in analog domain
  ba1=[0,1];
  aa1=[par(4),1];
  
  % xform using bilinear to digital domain
  % - we must use blinear rather than impinvar as otherwise
  % the zero freq gain deviates from one for small tc
  [b1,a1]=bilinear(ba1,aa1,100,1);
  % put into second-order section form
  sos1=tf2sos(b1,a1);
  % combine with butterworth filter
  [b1,a1]=sos2tf([sos1;sosb]);

  if(length(par)>4)
    % do again for second tc
    ba2=[0,1];
    aa2=[par(5),1];
    [b2,a2]=bilinear(ba2,aa2,100,1);
    % put into second-order section form
    sos2=tf2sos(b2,a2);  
    % combine with butterworth filter
    [b2,a2]=sos2tf([sos2;sosb]);  
    % apply
    y1=filter(b1,a1,y);
    y2=filter(b2,a2,y);
    y=(1-par(6))*y1+par(6)*y2;
  else
    y=filter(b1,a1,y);
  end
    
else
  % impose huge penality for negative tc
  y=par(4)*1e6*ones(size(y));
end

% When applying butterworth electronic filter we are seeking to
% simulate an analog (continuous time) filter with digital (discrete
% time) operation. This is not exactly possible and several
% approximations exist. The "bilinear" transform used in the standard
% "butter" call is apparently more appropriate for continuous signals
% while the butter(,,'s') form followed by impinvar transform is
% better for pulses.

% It appears that initially I used impinvar version here to extract
% timeconstants in July 2005 and these are the values which were
% distributed to the collaboration.

% Some time later (Jan 2006?) I switched in this function to use
% bilinear version of butter on the basis that that is what is used in
% deconv_scans so would be better to use it here also.

% In Sept 2006 coming back to this to fit 2006 season time consts I
% found that using bilinear and fitting 2005 data I got results
% systematically about 5ms higher than before but didn't understand
% why - see http://find.uchicago.edu/~quest/logbook/20060911/

% Now coming back to this in Oct 2006 I got deconv data version to
% work and find that that produces timeconsts much more in line with
% my old July 2006 values. So switching back to impinvar in here I
% find that I get close to the "old" values 

% setup butterworth low pass (electronics) filter
% (this is bilinear transform of analog filter)
%[b,a]=butter(6,20/50);

% some notes:
% to combine two filters
% [b1,a1]=butter(6,20/50);
% [b2,a2]=butter(6,10/50,'high');
% sos1=tf2sos(b1,a1)
% sos2=tf2sos(b2,a2)
% [b,a]=sos2tf([sos1;sos2])
%
% [b,a]=butter(6,20/50);
% [bc,ac]=butter(6,20*2*pi,'s');
% [b2,a2]=bilinear(bc,ac,100,20);
% b==b2; a==a2;
%
% [b3,a3]=impinvar(bc,ac,100);

return

% look at diff in freq and imp response when exp a to d conversion is
% done using impinvar and bilinear

tc=0.02;
% setup exp in analog domain according to Jamie recipe
ba=[0,1];
aa=[tc,1];

% xform using impinvar
[b1,a1]=impinvar(ba,aa,100);
% xform using bilinear
[b2,a2]=bilinear(ba,aa,100,1);

% setup butterworth electronics filter
% this is bilinear transform of analog filter
[bb,ab]=butter(6,20/50);

% combine filters
sos1=tf2sos(b1,a1);
sos2=tf2sos(b2,a2);
sosb=tf2sos(bb,ab);
[b1,a1]=sos2tf([sos1;sosb]);
[b2,a2]=sos2tf([sos2;sosb]);

figure(1)
[h1,t1]=impz(b1,a1,[],100);
freqz(b1,a1,[],100)

figure(2)
[h2,t2]=impz(b2,a2,[],100);
freqz(b2,a2,[],100)

figure(3)
plot(t1,h1,'b.-');
hold on
plot(t2,h2,'r.-');
hold off


% look at diff in freq and imp response when exp a to d conversion is
% done using impinvar and bilinear

% setup exp in analog domain according to Jamie recipe
ba1=[0,0.5];
aa1=[0.02,1];

ba2=[0,0.5];
aa2=[0.5,1];

% xform
[b1,a1]=impinvar(ba1,aa1,100);
[b2,a2]=impinvar(ba2,aa2,100);

[h1,t1]=impz(b1,a1,[],100);
[h2,t2]=impz(b2,a2,[],100);

n=length(h1);
h2(1:n)=h2(1:n)+h1;

figure(2)
plot(t2,h2,'b.-');
