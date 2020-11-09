function xaxis_time(h)
% xaxis_time(h)
%
% Convert xaxis of plot handle h to time

xl=xlim(h);
xr=xl(2)-xl(1);

if(xr>100)
  tv=datenum(2010,1:36,1);
  ds=12;
end
if(xr<100&xr>15)
  [xg,yg]=meshgrid([1:36],[1:7:28]);
  tv=datenum(2010,xg,yg);
  ds=6;
end
if(xr<15&xr>1)
  tv=datenum(2010,1,1)+[1:1000];
  ds=6;
end
if(xr<1)
  tv=datenum(2010,1,1)+[1:(1/24):1000];
  ds=15;
end
ind=tv<xl(2)&tv>xl(1);

set(h,'XTick',tv(ind));
datetick(h,'x',ds,'keepticks','keeplimits') 

grid

%set(h,'XMinorTick','on');

return
