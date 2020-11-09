function [simdata b k] = rps_make_chopped_data(chopref,amplitude,offset,tau,nsamp,MCE)
% [simdata b k] = rps_make_chopped_data(chopref,amplitude,offset,tau,nsamp,MCE)

if ~exist('offset','var')
   offset = 0; 
end

if ~exist('tau','var')
    fprintf('No tau input. Defaulting to 1ms\n')
   tau = 0.001; 
end

if ~exist('nsamp','var')
   nsamp = 100; 
end

if ~exist('MCE','var')
   MCE = true; 
end

simdata = chopref.*amplitude+offset;

%Apply detector response
if tau~=0
% tt=(0:(nsamp-1))'/150;
% b=exp(-tt/tau);
%b = [b(end:-1:1);b];
b = gausswin(2*nsamp,(2*nsamp-1)/2/(tau*150/2)^(1/2));
b=b/sum(b);
simdata = filter2(b,simdata,'full');
simdata = simdata((nsamp):(size(chopref,1)+(nsamp-1)));
else
    b = [];
end

if MCE
%Apply mce/gcp response
load('rps_gcpkernel')
k = rps_gcpkernel.kernel(end:-1:1);
simdata = filter2(k,simdata,'full');
simdata = simdata(1:size(chopref,1));
else
    k = [];
end


