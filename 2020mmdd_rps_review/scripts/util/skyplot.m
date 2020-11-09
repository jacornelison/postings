function skyplot(option)
% skyplot()
%
% Set image display astronomy style
%
% Optional cell array of string option arguments
% defaults to cross hair on
% Available options:
% xhair - green cross hairs are overplotted at zero
% sym - color axis is forced symettric about zero
%
% e.g.: skyplot({'xhair','sym'})

if(~exist('option','var'))
  option{1}='xhair';
end

axis xy; axis equal; axis tight; set(gca,'XDir','reverse');

if(any(strcmp(option,'xhair')))
  hold on
  x=xlim; y=ylim;
  line([0,0],y,'Color','g','LineStyle',':');
  line(x,[0,0],'Color','g','LineStyle',':');
  hold off
end

if(any(strcmp(option,'sym')))
  x=max(abs(caxis));
  caxis([-x,x]);
end
