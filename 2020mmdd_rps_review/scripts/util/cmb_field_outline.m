function [x,y]=cmb_field_outline(type,fi)
% [x,y]=cmb_field_outline(type)
%
% Get x/y data delineating various sets of CMB fields
%
% dasi  = DASI fields 50% beam point
% dasio = DASI fields 10% beam point
% cbi   = CBI mosaic fields
% boom2 = elipse used in 2nd analysis
%
% For dasi cases extra arg fi defines field rows 'ABCDE'
% - defaults to 'ABCD'
%
% e.g: [x,y]=cmb_field_outline('dasio'); plotm(x,y,'k');

if(strcmp(type,'dasi')|strcmp(type,'dasio'))
  % DASI field coords
  o=ones(1,8);
  dec=[-67*o;-61*o;-55*o;-49*o;-43*o];
  ra=[21.5 22.5 23.5 0.5 1.5 2.5 3.5 4.5; ... % B
      22 23 0 1 2 3 4 5; ...                  % A
      22.5 23.5 0.5 1.5 2.5 3.5 4.5 5.5; ...  % C
      23 0 1 2 3 4 5 6; ...                   % D
      23.5 0.5 1.5 2.5 3.5 4.5 5.5 6.5];      % E

  % DASI primary beam points
  % 2.85 is 10% pow radius at 31 GHz
  % 3.60 is 1% pow radius at 31 GHz
  
  % 1.66 is 50% pow radius for agregate beam 26-36 GHz
  % 2.88 is 10% pow radius for agregate beam 26-36 GHz
  % 3.79 is 1% pow radius for agregate beam 26-36 GHz

  if(~exist('fi'))
    fi='ABCD';
  end
  
  ind=[];
  if(any(fi=='B'))
    ind(end+1)=1;
  end
  if(any(fi=='A'))
    ind(end+1)=2;
  end
  if(any(fi=='C'))
    ind(end+1)=3;
  end
  if(any(fi=='D'))
    ind(end+1)=4;
  end
  if(any(fi=='E'))
    ind(end+1)=5;
  end

  ra=ra(ind,:); dec=dec(ind,:);
  ra=ra(:); dec=dec(:);  
end

switch type
  case 'dasi'
    r=1.66*ones(size(ra));  
    [x,y]=scircle1(dec,ra*15,r);
  case 'dasio'
    r=2.88*ones(size(ra));  
    [x,y]=scircle1(dec,ra*15,r);
  case 'cbi'
    [x(:,1),y(:,1)]=rect([-4.5,-2.5],[40.5,44.5]);
    [x(:,2),y(:,2)]=rect([-4.5,-2.5],[130.5,134.5]);
    [x(:,3),y(:,3)]=rect([-4.5,-2.5],[220.5,224.5]);
    [x(:,4),y(:,4)]=rect([-4.5,-2.5],[310.5,314.5]);
  case 'boomo'
    % Numbers digitized by Kovac from some no longer
    % existing Coble plot:
    % http://astro.uchicago.edu/home/web/coble/b98_sky_cov.ps
    y=[58 115 120 140 142 150 160 161 179 173 162 140 142 110 100 ...
       40 38 38 43 40 42 38 38 45 46 44 51 58];
    x=[-29 -29 -31 -31 -42 -48 -48 -52 -50 -60 -64 -64 -60 -61 ...
       -62 -62 -61 -60 -56 -55 -52 -50 -48 -45 -42 -41 -31 -29];
  case 'boom1'
    % Numbers digitized by Kovac off the plot
    y=[70,70,96.56,96.5645,96.9529,97.3435,97.7366,98.1325,98.5314,98.9337,...
       99.3397,99.749,100.1643,100.5836,101.0082,101.4385,101.8751,102.3183,...
       102.7689,103.2273,103.6943,104.1707,104.6570,105.1543,105.6635,...
       106.1854,106.4,70];
    x=[-55,-35,-35,-35.0978,-35.9827,-36.8683,-37.7544,-38.6411,-39.5284,...
       -40.4160,-41.3041,-42.1926,-43.0815,-43.9706,-44.8600,-45.7497,...
       -46.6395,-47.5294,-48.4193,-49.3093,-50.1993,-51.0891,-51.9787,...     
       -52.8682,-53.7573,-54.6460,-55,-55];
  case 'boom2'
    [x,y]=ellipse(-46,85,12,20);
  case 'quad05'
    [x,y]=rect_rh([-43,-50],[5*15,6*15]);
  case 'bicep06'
    [x,y]=rect_rh([-45,-70],[-30,30]);
end
  
return

%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x,y]=rect(xc,yc)
% Given opposite corners return a rectangle on the sphere
x=[xc(1),xc(2),xc(2),xc(1),xc(1)];
y=[yc(1),yc(1),yc(2),yc(2),yc(1)];

[x,y]=track('gc',x,y);
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x,y]=rect_rh(xc,yc)
% Use rhumb lines instead of great circle
x=[xc(1),xc(2),xc(2),xc(1),xc(1)];
y=[yc(1),yc(1),yc(2),yc(2),yc(1)];

[x,y]=track('rh',x,y);
return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x,y]=ellipse(xc,yc,a,b)
% Draw an ellipse on the sphere

az=linspace(0,2*pi,360);
x=a*cos(az); y=b*sin(az);
r=sqrt(x.^2+y.^2);
az=atan2(y,x);

az=az*180/pi;
x=xc*ones(size(az)); y=yc*ones(size(az));
[x,y]=reckon('gc',x,y,r,az);
return
