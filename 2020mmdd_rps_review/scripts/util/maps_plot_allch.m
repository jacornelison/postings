function maps_plot_allch(map,p,ind,tag)
% maps_plot_allch(map,p,ind,tag)
%
% Plot full set of ind channel maps which may be a/b or sum/diff

setwinsize(gcf,1200,900); clf
for i=ind
  x=(i-1)/2;
    
  subplot(8,8,x*2+1)
  m=map.rdoff.map(:,:,i);
  imagesc(map.rdoff.x_tic,map.rdoff.y_tic,m);
  axis xy; axis tight; title(sprintf('%s %.1e',p.channel_name{i},rms(m(:))));
  set(gca,'XDir','reverse');
  
  subplot(8,8,x*2+2)
  m=map.rdoff.map(:,:,i+1);
  imagesc(map.rdoff.x_tic,map.rdoff.y_tic,m);
  axis xy; axis tight; title(sprintf('%s %.1e',p.channel_name{i+1},rms(m(:))));
  set(gca,'XDir','reverse');
end

gtitle(sprintf('%s src=%s dk=%.0f',tag,map.src,map.dk));

if(numel(strfind(tag,'through'))>0)
  mkpng(sprintf('reducplots/all/maps_%s_dk%.0f.png',map.src,map.dk));
end
if(numel(strfind(tag,'+'))==0)
  mkpng(sprintf('reducplots/%s/maps_%s_dk%.0f.png',tag,map.src,map.dk));
end
  
return
