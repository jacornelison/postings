function p=do_hwp_rot(mjd,p)
% p=do_hwp_rot(mjd,p)
%
% Rotate the detector polarization sensitivity angle (chi) according
% to the waveplate rotation angle as given in external file
%
% mjd is the time value as in d.t - the mjd day with fraction in units
% of days.
%
% e.g. p=do_hwp_rot(d.t(1),p);

% If only one receiver (i.e. we are B2) do nothing 
if(length(unique(p.rx))==1)
  return
end

% translate mjd to utc
utc=mjd+datenum(2005,5,17)-53507;

for i=sort(unique(p.rx))'
  
  % get the angle data for this rx
  hwpf=sprintf('aux_data/hwp/angle_rx%1d.csv',i);
  hwp=ParameterRead(hwpf);
  
  % look for line where requested date>start and <end
  starttime=datenum(hwp.starttime,'yyyymmdd HH:MM:SS');
  endtime=datenum(hwp.endtime,'yyyymmdd HH:MM:SS');
  
  % find the relevant line in file
  after=utc>starttime;
  before=utc<endtime;
  hit=and(after,before);
  if(sum(hit)>1)
    error(['Two lines apply to same date in file ',hwpf]);
  end
  
  % find the relevant entries in p structure
  j=p.rx==i;
  
  % modify the pol sensitivity angle for this rx
  rotang=hwp.angle(hit);
  p.chi(j)=p.chi(j)-rotang;
  
  disp(sprintf('chi for rx%d rotated by %f',i,rotang));
  
end

return
