function h=hplot(bc,n,opt,lim)
% hplot(bc,n,opt,lim)
%
% Plot a histogram as a line "paw style"
%
% opt = 'L' makes y axis log
% lim = lower/upper limits to plot a subset of the histogram

if(~exist('opt','var'))
  opt=' ';
end

% Bin width (assume that bins are equal width)
bw=bc(2)-bc(1);

if(~exist('lim','var'))
  lim=[bc(1)-bw/2,bc(end)+bw/2];
end

% Use limits to remove unwanted bins
ind=bc>lim(1)&bc<lim(2);
bc=bc(ind);
n=n(ind);

% calc bin edges from centers
be=bc(1)+bw/2:bw:bc(end);
be2=[be;be];
be2=[bc(1)-bw/2,be2(:)',bc(end)+bw/2];

% Set zeros to a small value for log(y) plots
n(n==0)=1e-99;

n2=[n;n];
n2=n2(:)';

optp=opt(opt~='L'&opt~='S');

if(~any(opt=='S'))
  h=plot(be2,n2,optp);
  xlim([be2(1),be2(end)]);
  
  m=max(n);
  if(any(opt=='L'))
    m=max(n)*2;
    ylim([0.7,m]);
    set(gca,'YScale','log');
  else
    ylim([0,max(n)*1.1]);
  end
else
  hold on
  h=plot(be2,n2,optp);
  hold off
end
