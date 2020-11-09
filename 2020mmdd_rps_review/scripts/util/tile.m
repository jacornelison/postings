function out=tile(axis,x,y,tile,weight)
% out=tile(axis,x,y,tile,weight)
%
% Stack up a bunch of usually smaller arrays at given coords
% to build a "mosaic"
%
% axis is tic val of out array - tile resolution assumed same
% x,y are the location vectors to place center of each tile
% tile is 2D or 3D with 3rd dim same length as x,y
% for even tile size center assumed in N/2+1 element
% for odd tile size center assumed to be in N/2+0.5 element
% weight is optional multiplier for each tile
%
% At the moment tile is added centered on nearest pixel
% location in out array. Could interp tile to sub pixel
% center location...

if(exist('weight','var')~=1)
  weight=[];
end

if(isempty(weight))
  weight=ones(size(x));
end

% size of tile and out
N=length(axis);
M=size(tile,1);

% Make output array
out=zeros(N);

for i=1:length(x)
  % Find indecies in out of center of tile
  %[vx,m]=min(abs(axis-x(i)));
  %[vy,l]=min(abs(axis-y(i)));    
  
  m=1+(x(i)-axis(1))/(axis(2)-axis(1));
  l=1+(y(i)-axis(1))/(axis(2)-axis(1));
  
  % Find the start and end indecies of regions in out to place tile
  l1=l-M/2; l2=l+M/2-1;
  m1=m-M/2; m2=m+M/2-1;
  
  % Round up makes it work for odd size tiles also
  l1=ceil(l1); l2=ceil(l2);
  m1=ceil(m1); m2=ceil(m2);
  
  % If tile falls right outside reject
  if(l1<N & l2>1 & m1<N & m2>1)
  
    % Test if tile falls partially outside correct as necessary
    l3=1; l4=M;
    m3=1; m4=M;
    if(l1<1)
      l3=-l1+2; l1=1;
    end;
    if(l2>N)
      l4=l4-(l2-N); l2=N;
    end;
    if(m1<1)
      m3=-m1+2; m1=1;
    end;
    if(m2>N)
      m4=m4-(m2-N); m2=N;
    end;
    
    % Enter tile with appropriate weight
    if(size(tile,3)==1)
      out(l1:l2,m1:m2)=out(l1:l2,m1:m2)+tile(l3:l4,m3:m4)*weight(i);
    else
      out(l1:l2,m1:m2)=out(l1:l2,m1:m2)+tile(l3:l4,m3:m4,i)*weight(i);
    end
  end
end

return

