function dummy=plot_mce_masks(dead_lists,my_title)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Plot the MCE masks from .cfg files in dead_lists.  You will need to point
% to the files manually because they are not normally in the pipeline.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Grab the files with '.cfg' in their name.
listing=dir(dead_lists);
listing=struct2cell(listing);
listing=listing(1,:)';
my_cfg={};
for jj=1:length(listing)
    if ~isempty(strfind(listing{jj},'.cfg'))
        my_cfg=[my_cfg listing(jj)];
    end
end

% Loop through each of the .cfg files and get their masks.
my_mask=false(528,length(my_cfg));
for jj=1:length(my_cfg)
    my_mask_temp=get_mask_from_dead_cfg(strcat(dead_lists,my_cfg{jj}));
    % Go from 41*32 array to 528*1 ordering.
    my_mask(:,jj)=my_mask_temp(:);
end

% Truncate the extraneous '.cfg' ending for the labels on the plot.
my_cfg=regexprep(my_cfg,'.cfg','');

% Plot black on white.
imagesc(~my_mask(1:528,:)')
colormap(gray)
title(my_title)
set(gca,'TickLength',[0 0])
%set(gca,'XGrid','on')
grid on
set(gca,'XTick',1:33:528)
set(gca,'XTickLabel',(0:33:527)')
set(gca,'YTick',1:length(my_cfg))
set(gca,'YTickLabel',my_cfg)

dummy=true;