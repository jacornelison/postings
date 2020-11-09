function fs=tweak_field_scans(d,fs,turnfrac)
% fs=tweak_field_scans(d,fs,turnfrac)
%
% tweak fast start/end indices so that point to const velocity scan
% section and also so exact length such that lead trail fields can be
% point by point subtracted
%
% Modified Aug 6 2007 to regress azoff to find best "global" start/end
% and hopefully do a little better in the presence of jerky motion
%
% kludge paramater turnfrac includes part of turnarounds - Cynthia set
% to 1.13 in BICEP1 analysis

if ~exist('turnfrac','var') || isempty(turnfrac)
  turnfrac=1;
end

% assume samprate same for all scans - take from the first - in Keck
% data this line is highly unstable as there are almost equal numbers
% of two different integer values - skips can effect which is selected
% - there was a bug in the second loop which implictly meant that the
% last samprate of the first loop was used for all rep of the second
% loop - this showed up when Chris stopped the second loop from being
% active in scans containing skips
samprate=1/nanmedian(diff(d.t(fs.sf(1):fs.ef(1)))*3600*24);

for i=1:length(fs.s)
  
  sf=fs.sf(i); ef=fs.ef(i);
  azoff=d.azoff(sf:ef);

  n=round(samprate*turnfrac*fs.throw(i)/fs.rate(i));
    
  % look for point closest to -fs.throw/2
  [m,ind]=min(abs(azoff+turnfrac*fs.throw(i)/2));
  
  if(fs.inc(i)==1) 
    % for forward scan this is the beginning
    fs.sf(i)=sf+ind-1;
    fs.ef(i)=sf+ind-1+n;
  else
    % for backward scan this is the end
    fs.ef(i)=sf+ind-1;
    fs.sf(i)=sf+ind-1-n;
  end
end

% refine further by regressing versus linear motion - this
% refinement moves the sf/ef by only a few samples.
% (can't go straight to this algorithm as to start with fs.sf:fs.ef
% includes bits of accel/decel on the ends of the scans which mess up
% the fits)
if(1)
  for i=1:length(fs.s)
    
    sf=fs.sf(i); ef=fs.ef(i);
    azoff=d.azoff(sf:ef);
    
    n=round(samprate*turnfrac*fs.throw(i)/fs.rate(i));
    
    % skipped samples causes regression to fail or be unreliable; don't regress these
    % data, they must be cut in mapping
    dazoff=abs(azoff(2:end)-azoff(1:end-1));
    lothresh=median(dazoff)/2;
    hithresh=median(dazoff)*2;
    if(~any(dazoff<lothresh | dazoff>hithresh));

      m=ef-sf+1;
      X=[ones(m,1) [1:m]'];
      b=regress(azoff,X);
      % find point corresponding to -fs.throw/2
      ind=round((-fs.throw(i)*turnfrac/2-b(1))/b(2));
      
      if(fs.inc(i)==1) 
	% for forward scan this is the beginning
	fs.sf(i)=sf+ind-1;
	fs.ef(i)=sf+ind-1+n;
      else
	% for backward scan this is the end
	fs.ef(i)=sf+ind-1;
	fs.sf(i)=sf+ind-1-n;
      end
    
    end
  end
end

return
