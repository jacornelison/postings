function hplotb(be,n,opt)
% hplotb(be,n,opt)
%
% Plot a histogram from given bin edges and contents
% as a line "paw style"
%
% opt = 'L' makes y axis log

if(~exist('opt','var'))
  opt=' ';
end

% Set zeros to a small value for log(y) plots
n(n==0)=1e-99;

% Double up all but first and last be entries
be2=be(2:end-1);
be2=[be2;be2];
be2=[be(1),be2(:)',be(end)];

n2=[n;n];
n2=n2(:)';

%plot(be2,n2,opt(opt~='L'));
semilogx(be2,n2,opt(opt~='L'));
xlim([be2(1),be2(end)]);

m=max(n);
if(any(opt=='L'))
  ylim([3e-1,max(n)*2]);
  set(gca,'YScale','log');
else
  ylim([0,max(n)*1.1]);
end
