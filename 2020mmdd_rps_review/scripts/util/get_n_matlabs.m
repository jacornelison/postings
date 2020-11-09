function nm=get_n_matlabs()
% function nm=get_n_matlabs()
% return the number of matlab session running
% for your on the current host.
  username=whoami();
  [status,nm] = system_safe(['top -n 1 -b -u ',username,' | grep MATLAB | wc -l']);
  nm=str2num(nm);
return
