function snow=plot_lfbsnow(tags,showfig,dosave)
%plot_lfbsnow Summary of this function goes here
%   Detailed explanation goes here

if(~exist('tags','var'))
    tags=[];
end
if(isempty(tags))
    disp('Error: no tags to process');
    return;
end
if(~iscell(tags))
    tags={tags};
end
if(~exist('showfig','var'))
    showfig=[];
end
if(isempty(showfig))
    showfig=1;
end
if(~exist('dosave','var'))
    dosave=[];
end
% if(isempty(dosave))
%     dosave=0;
% end


if(showfig==0)
    set(0,'DefaultFigureVisible','off');
end

snow=lfbsnow(tags);

numpks=0;
if(length(snow.acc)>=3)
  [pks locs]=findpeaks(snow.acc,'minpeakdistance',3,'minpeakheight',1700);
  if(~isempty(locs))
    numpks=length(locs);
  end
else
  snow.t=0;
  snow.acc=0;
end
snwx=1705*ones(1,length(snow.t));

if(numpks>5)
  [sval sidx]=sort(pks,'descend');
  [locs lidx]=sort(locs(sidx(1:5)));
  pks=sval(lidx);
  numpks=5;
end

clf
setwinsize(gcf,1400,600);

subplot(2,5,1:5)
plot(snow.t,snow.acc,'k-*');
hold on;
plot(snow.t,snwx,'r-');
ylim([1000 6000]);
datetick('x','HH:MM','keepticks');
grid on;
ylabel('snow acc (arb)');
title(['Snow Accumulation on Lower Forebaffle from ' datestr(snow.t(1),'yyyy-mmm-dd HH:MM') ' to ' datestr(snow.t(end),'yyyy-mmm-dd HH:MM')]);

if(numpks>0)
  plot(snow.t(locs),pks,'rs');
   
  for ii=1:numpks
    subplot(2,5,5+ii);
    imagesc(snow.jpgs{locs(ii)}); colormap('Gray');  caxis([0 200]); axis square; axis off;
    title(datestr(snow.t(locs(ii)),'yyyy-mmm-dd HH:MM'));
  end
end

if(~isempty(dosave))
  if(ischar(dosave))
    mkpng([dosave '.png']);
  else
    mkpng(['reducplots/' tags{1}(1:6) '/snow_' tags{1}(1:8) '.png']);
    setpermissions(['reducplots/' tags{1}(1:6) 'snow_' tags{1}(1:8) '.png']);
  end
end

return

