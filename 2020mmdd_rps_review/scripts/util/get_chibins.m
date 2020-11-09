function chibins=get_chibins(bins)
% chibins=get_chibins(bins)
%
% return the sets of bandpower bins to be used when calculating chisq
% values
% dB Sept 1 2008, removed the 4 bins 3 to 6 after sim004.
  
if(~exist('bins','var'))
  bins=[];
end

if isempty(bins)
  bins=6;
end

chibins{1}=2:bins; % TT
chibins{2}=2:bins; % TE
chibins{3}=2:bins; % EE
chibins{4}=2:bins; % BB
chibins{5}=2:bins; % TB
chibins{6}=2:bins; % EB
  
  %disp(sprintf('you are using %0d ell bins',bins))

  chibins{7}=chibins{2};
  chibins{8}=chibins{5};
  chibins{9}=chibins{6};

  
return
