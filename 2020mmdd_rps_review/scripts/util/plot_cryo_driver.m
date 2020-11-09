exst=dir('cryoplots/2*.png');

if isempty(exst)
  numdate='20120212';
else
  dates=[{exst(:).name}];
  numdate=dates{end}(1:(end-4));
end

last = str2double(numdate);
lastnum=datenum(num2str(last),'yyyymmdd');
today=datenum(date,'dd-mmm-yyyy');

daystart=lastnum:7:today;

if isempty(daystart)
  last = max(str2double(numdate));
  lastnum=datenum(num2str(last),'yyyymmdd');
  daystart=lastnum;
end

dayend=daystart+6;

for ii=1:length(daystart)
  ds=str2double(datestr(daystart(ii),'yyyymmdd'));
  de=str2double(datestr(dayend(ii),'yyyymmdd'));
  plot_cryo(ds,de)
end

% Make cumulative plot
plot_cryo(20120213,de(end))

make_cryo_html(de(end))

system('chmod g+w cryoplots/*.html')
system('chmod g+w cryoplots/*.png')
