function ac=acunpack(ac)
% ac=acunpack(ac)
%
% Convert ac back into regular sparse matricies

% if ac appears to be packed unpack it
if(isfield(ac,'pixind'))
  fn=fieldnames(ac);
  fn=setdiff(fn,{'matsize','pixind'});
  
  for i=1:size(ac,1)
    for j=1:size(ac,2)
      pixind=getfield(ac(i,j),'pixind');
      for k=1:length(fn)
	x=sparse(zeros(ac(1).matsize));
	x(pixind)=getfield(ac(i,j),fn{k});
	ac(i,j)=setfield(ac(i,j),fn{k},x);
      end
    end
  end
  
  ac=rmfield(ac,{'pixind','matsize'});
end

return

