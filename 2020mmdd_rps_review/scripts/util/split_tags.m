function [tags_1,tags_2,thresh]=split_tags(tags,cutstats)
% given a tag list and a cut statistics, this returns the cut threashold that will 
%  yield 2 evenly filled taglists
% length of tags and aof cutstats must be the same
  
  
  n_tags=length(tags);
  thresh_start=median(cutstats);
  
  condition=10;
  inc=1;
  thresh=thresh_start;
  
  while condition  ~= 0
    
    q=find(cutstats < thresh);
    condition=(length(q) - n_tags/2);
    
    if condition > 0
      thresh=0.99*thresh;
    elseif condition < 0
      thresh=1.01*thresh;
    else
      %do nothing
    end
    
    
  end
  
  
  tags_1=tags(find(cutstats < thresh));
  tags_2=tags(find(cutstats > thresh));
  
  return
  
  
  
  