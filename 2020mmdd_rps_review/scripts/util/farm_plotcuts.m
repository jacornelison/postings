function farm_plotcuts(tags)
% farm_plotcuts(tags,ntagsperjob)
% farms reduc_plotcuts per phase.  
% specify job details:

username=whoami();
GROUP=['/BICEP2_cuts_' username];
queue='short_serial';
system(['bgadd ' GROUP]);
system(['bgmod -L 20 ' GROUP]);
args=' -R rusage[mem=6000]';

for i=1:length(tags) 
  farmit('farmfiles/',['reduc_plotcuts({''' tags{i} '''},''calc'');'],...
    'group',GROUP,'queue',queue,'args',args);
end
