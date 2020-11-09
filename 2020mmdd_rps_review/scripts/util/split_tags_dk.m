function [tags_1]=split_tags_dk(tags,dk)
  
  %given a list of tags, splits the list into a sublist with only the dk angles spacified.
 % dk angles must be specified in same manner as last part of tag ( ie 315, 000, 135, 180)  
  
 n_tags=length(tags);
 k=1;
 
 for i=1:length(tags)
 
   if strfind(tags{i}(end-2:end),dk)
     idx(k)=i;
     k=k+1;
   end
 
 end
 
tags_1=tags(idx);
l=length(tags_1);
disp(sprintf('you have reduced your list to %d tags at dk angle %s',l,dk));

 
  return
  
  
  
  