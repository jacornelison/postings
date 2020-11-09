function sp=subplot_grid(n,m,k,dox,doy,parentfig)
% subplot_grid(n,m,x,dox,doy,parentfig)
%
% version of subplot which stacks the plots with no spacing in x direction, y
% direction, or both
%
% e.g. subplot_grid(3,2,1);
%      plot(1:10);
%      subplot_grid2(3,2,1);
%
% dox = controls spacing in x (default = false)
% doy = controls spacing in y (default = true)


if ~exist('dox','var') || isempty(dox)
  dox=false;
end

if ~exist('doy','var') || isempty(doy)
  doy=true;
end



% Run this to get plot x-spacing
%for n=2:10                            
%  subplot(n,1,n);p1=get(gca,'Position');
%  subplot(n,1,n-1);p2=get(gca,'Position');
%  w(n)=p2(2)-(p1(2)+p1(4));
%end


switch n
  case 1
    st=0;
  case 2
    st=.1327;
  case 3
    st=.0839;
  case 4
    st=.0613;
  case 5
    st=.0483;
  case 6
    st=.0399;
  case 7
    st=.0340;
  case 8
    st=.0296;
  case 9
    st=.0262;
  case 10
    st=.0235;
  case 11
    st=.0173;
  otherwise
    error('1<=n<=10');
end

% Run this to get plot x-spacing
%for n=2:10                            
%  subplot(1,n,1);p1=get(gca,'Position');
%  subplot(1,n,2);p2=get(gca,'Position');
%  w(n)=p2(1)-(p1(1)+p1(3));
%end



x=(m-1)/m;
 
switch m
  case 1
    w=0;
  case 2
    w=.1057;
  case 3
    w=.0674;
  case 4
    w=.0495;
  case 5
    w=.0391;
  case 6
    w=.0323;
  case 7
    w=.0275;
  case 8
    w=.0240;
  case 9
    w=.0212;
  case 10
    w=.0191;
  otherwise
    error('1<=m<=10');
end

if ~exist('parentfig','var')
  parentfig = [];
end

if ~isempty(parentfig)
  sp=subplot(n,m,k,'align','Parent',parentfig); box on
else
  sp=subplot(n,m,k,'align');
end

p=get(gca,'Position');

if dox
  xpos = mod(max(k)-1,m)+1;
  p(3) = p(3) + w*(m-1)/m - 0.003;
  p(1) = p(1) - (xpos-1)*w*(1-(m-1)/m);
end

if doy
  ypos = n-floor((max(k)-1)/m);
  p(4) = p(4) + st*(n-1)/n - 0.003;
  p(2) = p(2) - (ypos-1)*st*(1-(n-1)/n);
end

%% leave some space for axis labels if not gridding in that dimenison
%if ~doy
%  p(4)=p(4)-0.1/n;
%end
%
%if ~dox
%  p(3)=p(3)-0.1/m;
%end

set(gca,'Position',p);

return
