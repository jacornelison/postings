function track_job_resources(jobidentifier,username)
% function track_job_resources(jobidentifier,username)
% make a histogram of the time and the mem used by finished jobs found in sacct
% that have 'jobidentifier' in the jobname
% 
% username: if not handed in automatically uses your username
%
% Examples:
% >> track_job_resources('1351')
% >> track_job_resources('1351','srh,sstokes,pryke')
  
ff = '~/job_recources.txt';

if ~exist('username','var')
  % fetch the information with the acct command:
  username = whoami();
end

system_safe(['sacct --format="jobid,jobname%60,state,start,elapsed,TimeLimit,MaxVMSize,ReqMem,exitcode" -u ',username,' | grep -A 1 ',jobidentifier,' | grep '' batch '' > ',ff])


% extract used mem and time
m = get_mem(ff);
t = get_temp(ff);

% plot it
figure(1); clf;
set(gcf,'Position',[0 0 840 520])
subplot(2,2,1)
stairhist(m,linspace(0,max(m)*1.05,50))
ylabel('Entries');xlabel('Mem used [GB]');
xlim([0,max(m)*1.05]);

subplot(2,2,2)
stairhist(t,linspace(0,max(t)*1.05,50))
ylabel('Entries');xlabel('Runtime [min]');
xlim([0,max(t)*1.05]);

subplot(2,2,3)
plot(m,t,'.')
xlabel('Mem used [GB]');ylabel('Runtime [min]');

% get rid of the txt file used for the sacct dump:
system_safe(['rm -f ',ff]);

return

function ms = get_mem(ff)
[s,d] = system_safe(['cat ',ff,' | awk ''{print substr($0,120,20)}'' ']);
d=strread(d,'%s');
ms = zeros(size(d));
for ii=1:length(d)
  m = d{ii}(1:end-1);
  u = d{ii}(end);
  f = 1024;
  switch u
    case 'M'
      f=f^2;
    case 'G'
      f=f^3;
  end
  ms(ii) = str2num(m)*f/1024^3;
end
return

function ts = get_temp(ff)
[s,d] = system_safe(['cat ',ff,' | awk ''{print substr($0,108,8)}'' ']);
d=strread(d,'%s');
ts = zeros(size(d));
for ii=1:length(d)
  h = str2num(d{ii}(1:2)); m = str2num(d{ii}(4:5)); s = str2num(d{ii}(7:8));
  ts(ii) = h*60 + m + s/60;
end
return
  
  
