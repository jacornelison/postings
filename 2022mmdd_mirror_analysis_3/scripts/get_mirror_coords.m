function [x_mirr, y_mirr] = get_mirror_coords(dk,x,y,drumangle,mount,mirror,PLOT)
%
%
% Point at which the central ray for a given detector hits the mirror.

if ~exist('PLOT','var')
    PLOT = false;
end

if length(dk)==1 & length(x)~=1
    dk = repmat(dk,length(x),1);
    dk = reshape(dk,size(x,1),size(x,2));
end
az = dk*0;
el = dk*0;
nsamp = length(az);

% Calculate aperture position, boresight pointing, and boresight
% orientation.
[pos, pnt, ort] = kbmp_mount(az, el, dk, mount, drumangle);

% Right-handed coordinates:
%   pnt = B3, ort = B1, yhat = B3 cross B1
ort2 = cross(pnt, ort, 2);

% Calculate mirror position and orientation.
%[mpos, a, mtilt] = kbmp_mount(dk*0, dk*0, dk*0, mount, 0);

[mpos, mnorm] = kbmp_mirror(az, el, mount, mirror);

% Mirror X/Y is defined as we're looking at the face of the mirror.
% +X is going toward the right of the mirror and
% +Y is going toward the top of the mirror.
% The mirror coordinates with no tilt or roll thus will have
% same +X convention and Y is defined in the right-handed sense (Norm x X).
mort = repmat([1, 0, 0], nsamp, 1);

% Apply mirror roll.
euler = zeros(nsamp, 3);
euler(:,2) = mirror.roll;
euler(:,1) = 90;
mort = rotate_3d(mort, euler);

% Apply mirror tilt.
euler = zeros(nsamp, 3);
euler(:,2) = mirror.tilt;
mort = rotate_3d(mort, euler);

% Move mount to encoder az/el position.
euler = zeros(nsamp, 3);
euler(:,2) = 90 - el;
euler(:,1) = -1 * az;
mort = rotate_3d(mort, euler);

mort2 = cross(mnorm,mort,2);

%
% We now have a plane defined by the mirror norm.
[x_mirr, y_mirr, h] = deal(NaN(size(x)));
[port, port2,ppnt,Q] = deal(NaN(size(mnorm)));
[theta,r] = cart2pol(x,y);
for i = 1:length(r)
    % Convert pointing to cartesian
    [port(i,:), port2(i,:), ppnt(i,:)] = rtheta_rotation(r(i), theta(i)*180/pi, ...
        ort(i,:), ort2(i,:), pnt(i,:));
    
    
    % Find the intersection point of the plane in mirror coords
        
    Q(i,:) = (pos(i,:)-mpos(i,:));
    denom = dot(ppnt(i,:),mnorm(i,:));
    h(i) = -dot(mnorm(i,:),Q(i,:))/denom;
    
    x_mirr(i) = dot(mort(i,:),-1*(Q(i,:)+ppnt(i,:)*h(i)));
    y_mirr(i) = dot(mort2(i,:),-1*(Q(i,:)+ppnt(i,:)*h(i)));
    
end


%%
if PLOT
    %keyboard()
    figure(1)
    clf; hold on;
    
    uscale = 0.2;
    % mirror object
    mnormrot = (rodmatrix(90,mort(1,:))*mnorm(1,:)')';
    quiv3(mpos,mnormrot,'showarrowhead','off','Color','k')
    %quiver(mpos(1,1),mpos(1,3),mnormrot(1),mnormrot(3),'showarrowhead','off','Color','k')
    mnormrot = (rodmatrix(-90,mort(1,:))*mnorm(1,:)')';
    %quiver(mpos(1,1),mpos(1,3),mnormrot(1),mnormrot(3),'showarrowhead','off','Color','k')
    quiv3(mpos,mnormrot,'showarrowhead','off','Color','k')
    
    % Mirror Normal
    %quiver(mpos(1,1),mpos(1,3),mnorm(1,1)*uscale,mnorm(1,3)*uscale,'Color','r')
    %quiver(mpos(1,1),mpos(1,3),mort2(1,1)*uscale,mort2(1,3)*uscale,'Color','g')
    quiv3(mpos,mnorm*uscale,'Color','r')
    quiv3(mpos,mort2*uscale,'Color','g')
    quiv3(mpos,mort*uscale,'Color','b')
    
    % aperture
%     quiver(pos(1,1),pos(1,3),pnt(1,1)*uscale,pnt(1,3)*uscale)
%     quiver(pos(1,1),pos(1,3),ort(1,1)*uscale,ort(1,3)*uscale)
%     quiver(pos(1,1),pos(1,3),ort2(1,1)*uscale,ort2(1,3)*uscale)
    quiv3(pos,pnt*uscale,'Color','m')
    quiv3(pos,ort*uscale,'Color','m')
    quiv3(pos,ort2*uscale,'Color','m')

    
    for i = 1:length(r)
        %quiver(pos(i,1),pos(i,3),ppnt(i,1)*2,ppnt(i,3)*2,'showarrowhead','off')
        quiv3(pos(i,:),ppnt(i,:)*2,'showarrowhead','off')
    end
    
    for i = 1:length(r)
        R = pos(i,:)+ppnt(i,:)*h(i);
        plot3(R(1),R(2),R(3),'rx')
        
    end
    
    xlabel('Cart X (m)')
    ylabel('Cart Y (m)')
    zlabel('Cart Z (m)')
    grid on
    view(0,0)
    %axis image
end

function quiv3(V1,V2,varargin)
quiver3(V1(1,1),V1(1,2),V1(1,3),V2(1,1),V2(1,2),V2(1,3),varargin{:})

