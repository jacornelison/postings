function sc=find_scans(x,f,sampratio,rmcorrupt)
% sc=find_scans(x,f,sampratio,rmcorrupt)
%
% find start and end of scans defined by offset x
% where start of scan has f set
%
% optional arg rmcorrupt tells to shrink down scan periods to avoid
% regions at start end which have been corrupted in
% deconv/lowpass operations

if(~exist('f','var'))
  f=[];
end

if(~exist('sampratio','var'))
  sampratio=[];
end

if(~exist('rmcorrupt','var'))
  rmcorrupt='no_end_rm';
end

if(isempty(f))
  f=ones(size(x));
end

% don't define a default sampratio - asking for trouble

sc.s=[]; sc.e=[];

% note that start/end points have been empirically tuned to hit true
% scan start/end as measured with d.antenna0.pmac.fast_az_pos etc
n=1; i=0;
while(i<length(x))
  i=i+1;
  l=0;
  if(x(i)~=0&f(i)>0) % found start - parse scan
    sc.s(n,1)=i+3;
    while(x(i)>=l&i<length(x)) % look for turnaround 
      l=x(i);
      i=i+1;
    end
    sc.e(n,1)=i+1;
    sc.inc(n,1)=1; % denis changed the name to inc
    n=n+1;
    sc.s(n,1)=i+3;
    %while(x(i)~=0&i<length(x)) % look for end
    while(x(i)<=l&x(i)~=0&i<length(x)) % look for end or turnaround
      l=x(i);
      i=i+1;
    end
    sc.e(n,1)=i+3;
    sc.inc(n,1)=0;
    n=n+1;
  end
end

% if last scan end at end run assume interrupted and kick out
if ~strcmp(rmcorrupt,'keepend')
  ind=sc.e<length(x);
  sc=structcut(sc,ind);
else
  ind=sc.s<length(x);
  sc=structcut(sc,ind);
end

% remove extra at start/end
if(strcmp(rmcorrupt,'end_rm'))
  sc.s=sc.s+2;
  sc.e=sc.e-2;
end

% Introduce sf,ef - start,end fast sample
% (In BICEP 1st fast sample is contemporaneous with 1st slow sample
% according to antenna0.time.utcfast/slow)
sc.sf=(sc.s-1)*sampratio+1;
sc.ef=(sc.e-1)*sampratio+1;

if strcmp(rmcorrupt,'keepend') & sc.e(end)>length(x)
  sc.ef(end)=length(x)*sampratio;
end

return
