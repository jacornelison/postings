function azelView = coordinate_system_draw(azelView,rthetaPixel,elazBoresight,dkdrum,ealphaEllipse,chiPol,rthetaPrime,doFig)
%  function azelView = coordinate_system_draw(azelView,rthetaPixel,elazBoresight,dkdrum,ealphaEllipse,doFig)
%  defaults:
%  azelView=[19.5,28];     the view point of the 3D figure
%  rthetaPixel=[30,160];   the angular distance/direction between the pixel and the boresight
%  elazBoresight=[45,95];  telescope pointing
%  dkdrum=[10,0];          dk and drum angle (drum not implemented)
%  ealphaEllipse=[0.9,36]; beam ellipse parameters, not c,p convention, alpha counted as chi
%  chiPol=30;              polarization direction of the channel
%  rthetaPrime=[30,270];   the angular distance/direction between the P and the P'
%  doFig=1;                make a figure panel, or (0) use the existing one
%  
%  for ealphaEllipse, rthetaPrime / chiPol one can input a third /second value, to switch off the 
%  display of that property, for instance rthetaPrime=[30,250,0], to switch off the prime
%  coordinate system

  if ~exist('azelView','var') | isempty(azelView) azelView=[19.5,28]; end
  if ~exist('rthetaPixel','var') | isempty(rthetaPixel) rthetaPixel=[30,160]; end
  if ~exist('elazBoresight','var') | isempty(elazBoresight) elazBoresight=[45,95]; end
  if ~exist('dkdrum','var') | isempty(dkdrum) dkdrum=[10,0]; end
  if ~exist('ealphaEllipse','var') | isempty(ealphaEllipse) ealphaEllipse=[0.9,36,1]; end
  if ~exist('chiPol','var') | isempty(chiPol) chiPol=[150,0,1]; end
  if ~exist('rthetaPrime','var') | isempty(rthetaPrime) rthetaPrime=[30,250,1]; end
  if ~exist('doFig','var') doFig=1; end
  
  if doFig
    clf
    set(gcf,'Position',[1922 186 1108 755]);
  end
  
  % that is the gray sphere in the back
  [X,Y,Z] = sphere(100);
  hold off;
  surf(X,Y,Z,'EdgeColor','none','LineStyle','none');
  hold on;
  
  % the lines have to lie on a radius thas is little wider than the sphere:
  sr=1.002;
  
  % these indicate Az=0,90 and el=0:
%    drawCircle(90,0,90,sr,'k--');
  drawCircle(90,0,90,sr,'color','k','linestyle','--');
  drawCircle(0,0,90,sr,'color','k','linestyle','--');
  drawCircle(0,90,90,sr,'color','k','linestyle','--');
  annotatePoint(30,0,1,'A=0',0,10,'k');
  annotatePoint(-10,90,1,'A=90',0,10,'k');
  
  % copies for convenience:
  K = dkdrum(1); % dk rotation angle
  r = rthetaPixel(1);
  theta = rthetaPixel(2);
  bs = elazBoresight; %the boresight (E, A)
    
  % few special points:
  annotatePoint(90,0,1,'Z');
  annotatePoint(0,0,1,'N');
  % the boresight
  annotatePoint(bs(1),bs(2),1,'B');
  % that is a circle at the pixel distance r around the boresight:
  drawCircle(bs(1),bs(2),r,sr,'color','k','linestyle','-');
  
  % indicate az and el of boresight
  drawArc(0,0,bs(2),90,sr,'color','b','linestyle','-'); % A
  drawArc(bs(1),bs(2),bs(1),180,sr,'color','b','linestyle','-'); % E
  annotatePoint(0,bs(2)/2,1,'A',0,15,'b');
  annotatePoint(bs(1)/2,bs(2),1,'E',0,15,'b');
  
  % dk rotation lines
  % drawArc starts a great circle at the first two points
  % and the goes the length of the third argument. The angle towards
  % which the arc is directed is given by the 4th argument, counted
  % counterclockwise when looked at the figure.
  drawArc(bs(1),bs(2),360,90,sr,'color','b','linestyle','-.'); % dk=0
  drawArc(bs(1),bs(2),r,90+K,sr,'color','b','linestyle',':'); % K
  % annotate the angle of the dk rotation.
  % draw Angle uses the first two arguments as center of a ring 
  % segment with radius of the 3rd argument. The it draws counterclockwise
  % between the direction towards the 4th and 5th argument, using the later to
  % place a (clumpsy) arrow head.
  drawAngle(bs(1),bs(2),r/1.5,90,90+K,sr,'b','K');
  
  % r and r extended, at the end of r is the Pixel:
  [p(1),p(2),rm(1),rm(2)]=drawArc(bs(1),bs(2),r,90+K-theta,sr,'color','r','linestyle','-');
  drawArc(bs(1),bs(2),r*2,90+K-theta,sr,'color','r','linestyle',':');
  annotatePoint(p(1),p(2),1,'P',1,15);
  annotatePoint(rm(1),rm(2),0.95,'r',0,15,'r');  
  drawAngle(bs(1),bs(2),r/3,90+K,90+K-theta,sr,'r','\theta');
    
  % x and y
  % the circle corresponding to r has a different angle at P
  % when counted from the zenith, calculate the angle towards
  % r at P with the sine formula for spherical trigonometry,
  % or matlab azimuth. The 180- make it count in the same direction
  % as used here.
  rAngle = 180-azimuth(p(1),-p(2),bs(1),-bs(2));

  % if we draw...
  % drawArc(p(1),p(2),60,180-rAngle,sr,'k--')
  % ...this will be on top of rextended
  % the new coordinate system is setup with x counted -theta from r
  [xe(1),xe(2)]=drawArc(p(1),p(2),20,rAngle+theta,sr,'color','k','linestyle','-'); % x
  drawArc(p(1),p(2),360,rAngle+theta,sr,'color','k','linestyle',':'); % x
  [ye(1),ye(2)]=drawArc(p(1),p(2),20,rAngle+theta-90,sr,'color','k','linestyle','-'); % y
  drawArc(p(1),p(2),360,rAngle+theta-90,sr,'color','k','linestyle',':');
  annotatePoint(xe(1),xe(2),0.95,'x''',0,15);
  annotatePoint(ye(1),ye(2),0.95,'y''',0,15);
  % and the angle between r and x:
  drawAngle(p(1),p(2),r/3,rAngle,rAngle+theta,sr,'r','-\theta',0.25);
  
  % draw an ellipse:
  if checkDraw(ealphaEllipse)
    drawEllipse(p(1),p(2),15*sqrt(1+ealphaEllipse(1)),sr,ealphaEllipse(1),90-rAngle+ealphaEllipse(2),'color','k','linestyle','--');
  end
  
  % draw polarization direction chi
  if checkDraw(chiPol)
    drawArc(p(1), p(2), 30, rAngle-chiPol(1)+180, sr,'color','m','linestyle','--');
    [c(1), c(2), cm(1), cm(2)] = drawArc(p(1), p(2), 30, rAngle-chiPol(1), sr,'color','m','linestyle','--');
    %annotatePoint(c(1), c(2), 1, '\chi', 1, 15, 'm');
    drawAngle(p(1), p(2), 20, rAngle, rAngle-chiPol(1), sr, 'm', '\chi');
  end
  
  % draw the stuff belonging to P'
  % first r' and P'
  if checkDraw(rthetaPrime)
    ppcol=[0,0.7,0];
    [p_prime(1),p_prime(2),ppm(1),ppm(2)]=drawArc(p(1),p(2),rthetaPrime(1),rAngle+theta-rthetaPrime(2),sr,'color',ppcol,'linestyle','-');
    drawArc(p(1),p(2),max(rthetaPrime(1)*1.5,25),rAngle+theta-rthetaPrime(2),sr,'color',ppcol,'linestyle',':');
    annotatePoint(p_prime(1),p_prime(2),1,'P''',1,15,ppcol);
%      keyboard
    drawAngle(p(1),p(2),17,rAngle+theta,rAngle+theta-rthetaPrime(2),sr,ppcol,'\theta''')
    
    if checkDraw(chiPol)
      % now where chi' is identified:
      drawAngle(p(1), p(2),25, rAngle+theta-rthetaPrime(2), rAngle - chiPol(1), sr, ppcol, '\chi''');
      p_prime_angle = 180-azimuth(p_prime(1),-p_prime(2),p(1),-p(2));
      % and chi' around P':
      chi_prime = theta + chiPol(1) - rthetaPrime(2);
      [cp(1), cp(2), cpm(1), cpm(2)] = drawArc(p_prime(1), p_prime(2), 30, p_prime_angle-chi_prime, sr,'color','m','linestyle',':');
      drawAngle(p_prime(1), p_prime(2), 10, p_prime_angle, p_prime_angle-chi_prime, sr, ppcol, '\chi''');    
    end
  end    
      
  xlim([-1,1]);
  ylim([-1,1]);
  zlim([-1,1]);
  xlabel('X');
  ylabel('Y');
  axis equal;
  set(gca,'visible','off');
  view(azelView(1),azelView(2));
  rotate3d on;
  colormap gray;
return

function doDraw = checkDraw(toCheck)
% use the third element to decide if this should be drawn or not
  doDraw=1;
  if numel(toCheck)==3 && ~toCheck(3)
    doDraw=0;
  end
return

function drawEllipse(el,az,alpha,r,e,rot,varargin)
  [cel,caz] = circOnSphere(el,az,alpha,0,360,e,rot);  
  r = ones(size(cel))*r;
  [x,y,z] = sph2xyz(cel,caz,r);
  plot3(x,y,z,'Linewidth',2,varargin{:})
return

function [elEnd,azEnd,elMid,azMid] = drawArc(el,az,alpha,kappa,r,varargin)
  % drawArc starts a great circle at the first two points
  % and the goes the length of the third argument. The angle towards
  % which the arc is directed is given by the 4th argument, counted
  % counterclockwise when looked at the figure.
  [cel,caz] = arcOnSphere(el,az,alpha,kappa);
  elEnd = cel(end);
  azEnd = caz(end);
  elMid = cel(round(numel(cel)/2));
  azMid = caz(round(numel(cel)/2));
  r = ones(size(cel))*r;
  [x,y,z] = sph2xyz(cel,caz,r);
  plot3(x,y,z,'Linewidth',2,varargin{:})
return

function [cel,caz] = arcOnSphere(el,az,alpha,kappa,numberOfPoints)
  if ~exist('numberOfPoints','var'), numberOfPoints=100; end

  %source point (theta,phi)
  [theta,phi] = elaz2thetaphi(el,az);
  
  [n(1),n(2),n(3)] = tpr2xyz(theta,phi,1);
  n=n';
  
  Or1t = 90-theta;
  Or1p = phi+180;
  [Or1(1),Or1(2),Or1(3)]=tpr2xyz(Or1t,Or1p,1);
  Or1 = Or1';
  
  Or2 = cross(n,Or1);
  Or1 = Or1/norm(Or1);
  Or2 = Or2/norm(Or2);
  
  theta = deg2rad(theta);
  phi = deg2rad(phi);
  alpha = deg2rad(alpha);
  kappa = deg2rad(kappa);
  
  i = 1;
  for alp = linspace(0,alpha,numberOfPoints)
    K = n*cos(alp) + Or1*sin(alp)*cos(kappa) + Or2*sin(alp)*sin(kappa);
    K = K/norm(K);
    ctheta(i) = acos(K(3));
    cphi(i) = atan2(K(2),K(1));
    i=i+1;
  end
  ctheta = rad2deg(ctheta);
  cphi = rad2deg(cphi);
  [cel,caz] = elaz2thetaphi(ctheta,cphi);
return


function [theta,phi] = elaz2thetaphi(el,az)
  theta = 90-el;
  phi = -az;
return

function drawAngle(el,az,alpha,kappaS,kappaE,r,col,t,annR)
  % annotate the angle of the dk rotation.
  % draw Angle uses the first two arguments as center of a ring 
  % segment with radius of the 3rd argument. Then it draws counterclockwise
  % between the direction towards the 4th and 5th argument, using the later to
  % place a (clumpsy) arrow head.
  if ~exist('annR','var') annR=0.5; end
  [cel,caz] = circOnSphere(el,az,alpha,kappaS,kappaE);  
  r = ones(size(cel))*r;
  [x,y,z] = sph2xyz(cel,caz,r);
  plot3(x,y,z,'Linewidth',2,'color',col)
  plot3(x(end-3),y(end-3),z(end-3),'^','MarkerSize',6,'MarkerFaceColor',col,'color',col)
  elMid = cel(round(numel(cel)*annR));
  azMid = caz(round(numel(cel)*annR));
  annotatePoint(elMid,azMid,0.95,t,0,15,col)
return

function drawCircle(el,az,alpha,r,varargin)
  [cel,caz] = circOnSphere(el,az,alpha);  
  r = ones(size(cel))*r;
  [x,y,z] = sph2xyz(cel,caz,r);
  plot3(x,y,z,'Linewidth',1,varargin{:})
return

function [x,y,z] = sph2xyz(el,az,r)
  [theta,phi] = elaz2thetaphi(el,az);
  [x,y,z] = tpr2xyz(theta,phi,r);
return

function [x,y,z] = tpr2xyz(theta,phi,r)
  x = r.*sind(theta).*cosd(phi);
  y = r.*sind(theta).*sind(phi);
  z = r.*cosd(theta);
return


function annotatePoint(el,az,r,t,doMark,Fontsize,col)
  if ~exist('doMark','var') doMark=1; end
  if ~exist('Fontsize','var') Fontsize=15; end
  if ~exist('col','var') col='k'; end
  [x,y,z] = sph2xyz(el,az,r*1.1);
  text(x,y,z,t,'Fontsize',Fontsize,'HorizontalAlignment','center','color',col)
  if doMark
    [x,y,z] = sph2xyz(el,az,r);
    plot3(x,y,z,'ko','MarkerSize',5,'MarkerFaceColor','k','color',col)
  end
return

function [cel,caz] = circOnSphere(el,az,alpha,kappaS,kappaE,e,rot,numberOfPoints)
  if ~exist('kappaS','var'), kappaS=0; end
  if ~exist('kappaE','var'), kappaE=360; end
  if ~exist('e','var'), e=0; end
  if ~exist('rot','var'), rot=0; end
  if ~exist('numberOfPoints','var'), numberOfPoints=100; end
  
  %source point (theta,phi)
  [theta,phi] = elaz2thetaphi(el,az);
  
  [n(1),n(2),n(3)] = tpr2xyz(theta,phi,1);
  n=n';
  
  Or1t = 90-theta;
  Or1p = phi+180;
  [Or1(1),Or1(2),Or1(3)]=tpr2xyz(Or1t,Or1p,1);
  Or1 = Or1';
  
  Or2 = cross(n,Or1);
  Or1 = Or1/norm(Or1);
  Or2 = Or2/norm(Or2);
  
  theta = deg2rad(theta);
  phi = deg2rad(phi);
  alpha = deg2rad(alpha);
  kappaS=deg2rad(kappaS);
  kappaE=deg2rad(kappaE);
  rot = deg2rad(rot);
  % alphamin = alpha*sqrt(1-e^2);
  alphamin = alpha*sqrt((1-e)/(1+e));
    
  i = 1;
  for kappa = linspace(kappaS,kappaE,numberOfPoints)
    calpha = alpha*alphamin/sqrt((alpha*cos(kappa+rot))^2+(alphamin*sin(kappa+rot))^2);
    K = n*cos(calpha) + Or1*sin(calpha)*cos(kappa) + Or2*sin(calpha)*sin(kappa);
    K = K/norm(K);
    ctheta(i) = acos(K(3));
    cphi(i) = atan2(K(2),K(1));
    i=i+1;
  end
  ctheta = rad2deg(ctheta);
  cphi = rad2deg(cphi);
  [cel,caz] = elaz2thetaphi(ctheta,cphi);
return


