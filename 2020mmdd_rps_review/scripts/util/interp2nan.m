function zi=interp2nan(x,y,z,xi,yi,varargin)
% interp2nan(x,y,z,xi,yi)
%
% interp2 fails when xi/yi are NaN which is silly as it should just
% return NaN for the interp val

ind=~isnan(xi)&~isnan(yi);

zi=NaN*zeros(size(xi));
zi(ind)=interp2(x,y,z,xi(ind),yi(ind),varargin{:});

return
