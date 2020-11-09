function xp=demod_data(x,y)
% xp=demod_data(x,y)
%
% demodulate the data in the columns of x using the squarewave demod
% template y

disp('demod_data');

% make y logical
y=y>mean(y);

% find the on states
fon=find_blk(y==1);

% find the off states
fof=find_blk(y==0);

% make sure off comes first
if(fon.s(1)<fof.s(1))
  fon=structcut(fon,2:length(fof.s));
end

% make sure on/off same length (and avoid end)
l=min([length(fon.s),length(fof.s)]);
fon=structcut(fon,1:l-1);
fof=structcut(fof,1:l-1);

% construct quadrature waveforms

% number of points in each pair
l=(fon.e-fon.s)+(fof.e-fof.s)+2;
% number of points in quarter phase
qp=l/4;
% take fractional part
f=rem(qp,1);
% When few samples per cycle quarter phase is poorly defined - to get
% best approx on average to quadrature waveform round down/up randomly
% with probability of fractional part of number of cycles per quarter phase.
% If frac part less than uniform random number add one
qp=floor(qp)+double(rand(size(qp))<f);

% add the quarter phase to each start/end pointer
fonq.s=fon.s+qp; fonq.e=fon.e+qp;
fofq.s=fof.s+qp; fofq.e=fof.e+qp;

% make on one shorter than off so last on has braketing offs
l=length(fon.s);
fon=structcut(fon,1:l-1);
fonq=structcut(fonq,1:l-1);

xp=zeros(size(x));

% take quadrature sum of diff of on and braketing off states
for i=1:(length(fon.s)-1)
  
  % out of phase part
  indon=fonq.s(i):fonq.e(i);
  indof=[fofq.s(i):fofq.e(i),fofq.s(i+1):fofq.e(i+1)];
  indal=[indon,indof];
  mvc=mean(x(indof,:),1)-mean(x(indon,:),1);

  % in phase part
  indon=fon.s(i):fon.e(i);
  indof=[fof.s(i):fof.e(i),fof.s(i+1):fof.e(i+1)];
  indal=[indon,indof];  
  mvs=mean(x(indof,:),1)-mean(x(indon,:),1);
  
  % store quadrature diff value in in phase on and off slots
  xp(indal,:)=repmat(pyth(mvs,mvc),[length(indal),1]);
end

return
