function d=lowpassfilt(d,sc,el)
% d=lowpassfilt(d,sc,el)
%
% low pass filter timestream data
% el is cell array of strings containing d struct fields to be low passed

disp('lowpassfilt...')

% prepare Butterworth low pass filter
[b,a]=butter(3,6/50);

% apply to each scan
for i=1:length(sc.sf)
  s=sc.sf(i); e=sc.ef(i);
  for j=1:length(el)
    eval(sprintf('d.%s(s:e,:)=filtfilt(b,a,d.%s(s:e,:));',el{j},el{j}));
  end
end

return
