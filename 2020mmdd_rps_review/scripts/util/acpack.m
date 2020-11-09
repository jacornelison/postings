function ac=acpack(ac)
% ac=acpack(ac)
%
% Convert fields of ac into lists of pixval and a list of pixel
% indices. This results in factor 2 smaller compressed files on disk.

if ~isfield(ac,'pixind')
  
  fn=fieldnames(ac);
  
  % record the map size to use when reinstating
  ac(1).matsize=size(ac(1).wsum);
  
  for i=1:size(ac,1)
    for j=1:size(ac,2)
      % find the indecies of the non-zero elements for this pair
      ac(i,j).pixind=find(getfield(ac(i,j),'wsum')~=0);
      for k=1:length(fn)
        % get this field
        x=getfield(ac(i,j),fn{k});
        % store a vector list of the non-zero elements - make this full
        % rather than sparse
        ac(i,j)=setfield(ac(i,j),fn{k},full(x(ac(i,j).pixind)));
      end
    end
  end
  
end

return
