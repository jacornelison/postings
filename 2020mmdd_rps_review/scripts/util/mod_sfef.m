function varargout=mod_sfef(n,varargin)
% varargout=mod_sfef(n,varargin)
%
% If have pointer structures gen by find_blk or find_scans and then
% downsample or use meanof_n_fastreg on the d structure the .sf/ef
% pointers become invalid. This function makes them good again.
%
% e.g.
% d=meanofn_fastreg(d,20,length(d.tf));
% [rc,cr,en]=mod_sfef(20,rc,cr,en);

for i=1:length(varargin)
  x=varargin{i};
  if(~isfield(x,'sf'))
    names=fieldnames(x);
    for i=1:length(names)
      eval(sprintf('x.%s=mod_sfef(n,x.%s);',names{i},names{i}));
    end
  else
    x.sf=(x.sf-1)/n+1;
    x.ef=(x.ef-1)/n+1;
  end
  varargout{i}=x;
end

return
