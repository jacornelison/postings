function plot_legend(data,text_x,text_y,delta,textFontSize)
%plot_legend(data,text_x,text_y,delta,textFontSize)
% on a plot print the following:
% valid, nanmean, nanmedian, nanstd
text(text_x,text_y,['Valid: ',num2str(sum(~isnan(data)))],'units','normalized','Fontsize',textFontSize);
text(text_x,text_y-1*delta,['Mean: ' num2str(nanmean(data),'%.3f') ' '],'units','normalized','Fontsize',textFontSize);
text(text_x,text_y-2*delta,['Median: ' num2str(nanmedian(data),'%.3f') ' '],'units','normalized','Fontsize',textFontSize);
text(text_x,text_y-3*delta,['Std: ' num2str(nanstd(data),'%.3f') ' '],'units','normalized','Fontsize',textFontSize);
end